import os
import numpy as np
import pandas as pd
from scipy import interpolate

def calc_events(N: pd.Series, S: pd.Series) -> pd.Series:
    """
    The Kaplan-Meier estimator of survival function is given by
    S(t) = product_{i:t_i<=t}(1 - D_i / N_i)
    Where S is event-free probability, D is number of events, and N is number of observations
    Assume S(t+10) = S(t) * ((N(t+10) - D(t+10)) / N(t+10))
    solve for D
    """
    D = N * (1 - S / S.shift(1))
    return D

def calc_hazard_rate(D: pd.Series, N: pd.Series) -> pd.Series:
    # H is noted as events per patient-year
    H = D / (N * (10 / 12)) # time_interval = 10 months
    return H

def calc_variance(S: pd.Series, D: pd.Series, N: pd.Series, num: int) -> list:
    """
    Greenwood's formula:
    Var_S(t) = S(t)^2 * sum_{i:t_i<=t}(D_i / (N_i * (N_i - D_i)))
    """
    Var_S = [np.nan] # set variance equal to NaN when t = 0
    val = 0
    for i in range(1, num, 1):
        S_i = S.loc[i]
        D_i = D.loc[i]
        N_i = N.loc[i]
        val += D_i / (N_i * (N_i - D_i))
        Var_S_i = S_i ** 2 * val
        Var_S.append(Var_S_i)
    return Var_S

def get_S_draws(S: pd.Series, Var_S: list, num: int, ndraws: int, seed=123) -> pd.DataFrame:
    """
    Sample n draws of S(t) from Normal(mean=S(t), sd=sqrt(var_S(t))) for each t
    """
    random.seed(seed)
    percentiles = np.random.uniform(low=0, high=1, size=ndraws)
    df = pd.DataFrame({'t_0': [1]*ndraws})
    for i in range(1, num, 1):
        S_i = []
        dist = norm(loc=S.loc[i], scale=np.sqrt(Var_S[i]))
        for p in percentiles:
            S_i.append(dist.ppf(p))
        df[f't_{i*10}'] = np.clip(S_i, 0, df[f't_{(i-1)*10}']) # 0 <= S(t) <= S(t-1)  
    return df.transpose().reset_index(drop=True)

def get_H_draws(N: pd.Series, S_data: pd.DataFrame, ndraws: int) -> pd.DataFrame:
    df = pd.DataFrame()
    for draw in range(ndraws):
        S = S_data[draw]
        D = calc_events(N, S)
        H = calc_hazard_rate(D, N)
        df[f'draw_{draw}'] = H
    return df

def interp(t_data: pd.Series, hazard_data: pd.DataFrame, ndraws: int, duration=60, step=1/30) -> pd.DataFrame:
    t_target = np.arange(0, duration, step)
    t_step_start = np.arange(0, len(t_target), 1)
    t_step_end = t_step_start + 1
    df = pd.DataFrame({'time_since_entrance_start': t_step_start, 'time_since_entrance_end': t_step_end})
    for i in range(ndraws):
        h_data = hazard_data[f'draw_{i}']
        h = interpolate.interp1d(t_data[1:], h_data[1:], kind='nearest', fill_value='extrapolate') # piece-wise constant
        val = h(t_target)
        df[f'draw_{i}'] = val
    return df

def get_results(input_dir: str, survival_outcome_type: str, line: str, ndraws: int) -> pd.DataFrame:
    if survival_outcome_type == 'overall_survival':
        s_name = 'Survival probability, S(t)'
    else:
        s_name = 'Progression-free probability, P(t)'
    data = pd.read_excel(os.path.join(input_dir, f'{survival_outcome_type}_by_time.xlsx'),
                         sheet_name = line,
                         engine = 'openpyxl')
    t = data['Time in months, t']
    num = len(t)
    N = data['Number at risk, N(t)']
    S = data[s_name]
    D = calc_events(N, S)
    Var_S = calc_variance(S, D, N, num)
    S_draws = get_S_draws(S, Var_S, num, ndraws)
    H_draws = get_H_draws(N, S_draws, ndraws)
    res = interp(t, H_draws, ndraws)
    return res

def select_draws(input_dir: str, line: str, ndraws: int, idx_cols: list):
    """
    output 1000 draws where incidence > mortality for each draw
    """
    assert ndraws > 1000, 'Insert ndraws greater than 1000'
    mortality = get_results(input_dir, 'overall_survival', line, ndraws).set_index(indx_cols)
    incidence = get_results(input_dir, 'progression_free_survival', line, ndraws).set_index(idx_cols)
    unwanted_draws = []
    draws = list(range(ndraws))
    for draw in draws:
        if (incidence[f'draw_{draw}'] - mortality[f'draw_{draw}'] < 0).any():
            unwanted_draws.append(draw)
    if unwanted_draws:
        wanted_draws = np.setdiff1d(draws, unwanted_draws)
        selected_draws = np.random.choice(wanted_draws, size=1000, replace=False)
        selected_draws = [f'draw_{i}' for i in selected_draws]
        mortality = mortality.loc[:, selected_draws]
        incidence =  incidence.loc[:, selected_draws]
        reordered_draws = [f'draw_{i}' for i in range(1000)]
        mortality.columns = reordered_draws
        incidence.columns = reordered_draws
    return mortality.reset_index(), incidence.reset_index()


if __name__ == '__main__':
    idx_cols = ['time_since_entrance_start', 'time_since_entrance_end']
    input_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/treatment_landscape/braunlin_et_al_2020'
    output_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/cause_model_input'
    for line in ['First-line', 'Second-line', 'Third-line', 'Fourth-line', 'Fifth-line']:
        df_mortality = get_results(input_dir, 'overall_survival', line, 1000)
        df_incidence = get_results(input_dir, 'progression_free_survival', line, 1000)
        if (df_incidence.set_index(idx_cols) - df_mortality.set_index(idx_cols) < 0).any():
            df_mortality, df_incidence = select_draws(input_dir, line, 2000, idx_cols)        
        df_mortality.to_csv(os.path.join(output_dir, f'mortality {line}.csv'), index=False)
        df_incidence.to_csv(os.path.join(output_dir, f'incidence {line}.csv'), index=False)
