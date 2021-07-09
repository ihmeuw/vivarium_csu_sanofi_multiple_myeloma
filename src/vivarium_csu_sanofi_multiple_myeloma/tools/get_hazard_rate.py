import os
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.stats import norm

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

def get_S_draws(S: pd.Series, Var_S: list, num: int, ndraws, seed) -> pd.DataFrame:
    """
    Sample n draws of S(t) from Normal(mean=S(t), sd=sqrt(var_S(t))) for each t
    """
    df = pd.DataFrame({'t_0': [1]*ndraws})
    for i in range(1, num, 1):
        S_i = []
        np.random.seed(seed) 
        for draw in list(range(0,ndraws)):
            S_i.append(norm(loc=S.loc[i], scale=np.sqrt(Var_S[i])).ppf(np.random.uniform(0,1)))
        df[f't_{i*10}'] = np.clip(S_i, 0, df[f't_{(i-1)*10}']) # 0 <= S(t) <= S(t-1)  
    return df.transpose().reset_index(drop=True)

def get_H_draws(N: pd.Series, S_data: pd.DataFrame, ndraws) -> pd.DataFrame:
    df = pd.DataFrame()
    for draw in range(ndraws):
        S = S_data[draw]
        D = calc_events(N, S)
        H = calc_hazard_rate(D, N)
        df[f'draw_{draw}'] = H
    return df

def interp(ndraws, t_data: pd.Series, hazard_data: pd.DataFrame, duration=60, step=1/30) -> pd.DataFrame:
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

def get_results(input_dir: str, survival_outcome_type: str, line: str, ndraws: int, seed: int) -> pd.DataFrame:
    if survival_outcome_type == 'overall_survival':
        s_name = 'Survival probability, S(t)'
    else:
        s_name = 'Progression-free probability, P(t)'
    data = pd.read_excel(f'{input_dir}/{survival_outcome_type}_by_time.xlsx',
                         sheet_name = line,
                         engine = 'openpyxl')
    t = data['Time in months, t']
    num = len(t)
    N = data['Number at risk, N(t)']
    S = data[s_name]
    D = calc_events(N, S)
    Var_S = calc_variance(S, D, N, num)
    S_draws = get_S_draws(S, Var_S, num, ndraws, seed)
    H_draws = get_H_draws(N, S_draws, ndraws)
    res = interp(ndraws, t, H_draws)
    return res

def test_for_illogical_draws(pfs_results, os_results, line):
    data = os_results.set_index(['time_since_entrance_start','time_since_entrance_end']) < pfs_results.set_index(['time_since_entrance_start','time_since_entrance_end'])
    data = data.stack().reset_index().rename(columns={'level_2':'draw',0:'val'})
    data = data.loc[data.val==False]
    bad_draws = data.draw.unique()
    assert len(bad_draws) == 0, f"Error: OS hazard rate greater than PFS hazard rate for {line} draws: {bad_draws}"
    
if __name__ == '__main__':
    input_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/treatment_landscape/braunlin_et_al_2020'
    output_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/cause_model_input'
    treatment_lines = ['First-line', 'Second-line', 'Third-line', 'Fourth-line', 'Fifth-line']
    seeds = [9,7,5,2,8]
    for i in list(range(0,len(treatment_lines))):
        df_mortality = get_results(input_dir, 'overall_survival', treatment_lines[i], 1000, seeds[i])
        df_incidence = get_results(input_dir, 'progression_free_survival', treatment_lines[i], 1000, seeds[i])
        #test_for_illogical_draws(df_incidence, df_mortality, treatment_lines[i])
        df_mortality.to_csv(os.path.join(output_dir, f'mortality {treatment_lines[i]}.csv'), index=False)
        df_incidence.to_csv(os.path.join(output_dir, f'incidence {treatment_lines[i]}.csv'), index=False)
