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
    # R is noted as events per person-year
    H = D / (N * (10 / 12)) # time_interval = 10 months
    return H

def calc_variance(S: pd.Series, D: pd.Series, N: pd.Series, num: int) -> list:
    """
    Greenwood's formula:
    Var_S(t) = S(t)^2 * product_{i:t_i<=t}(D_i / (N_i * (N_i - D_i)))
    """
    Var_S = [np.nan] # set variance equal to NaN when t = 0
    Values = []
    for i in range(1, num, 1):
        S_i = S.loc[i]
        D_i = D.loc[i]
        N_i = N.loc[i]
        val = D_i / (N_i * (N_i - D_i))
        Values.append(val)
        Var_S_i = S_i ** 2 * np.product(Values)
        Var_S.append(Var_S_i)
    return Var_S

def get_S_draws(S: pd.Series, Var_S: list, num: int) -> pd.DataFrame:
    """
    Sample 1000 draws of S(t) from Normal(mean=S(t), sd=sqrt(var_S(t))) for each t
    """
    df = pd.DataFrame({'t_0': [1]*1000})
    for i in range(1, num, 1):
        S_i = np.random.normal(loc=S.loc[i], scale=np.sqrt(Var_S[i]), size=1000)
        df[f't_{i*10}'] = S_i
    return df.transpose().reset_index(drop=True)

def get_H_draws(N: pd.Series, S_data: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame()
    for draw in range(1000):
        S = S_data[draw]
        D = calc_events(N, S)
        H = calc_hazard_rate(D, N)
        df[f'draw_{draw}'] = H
    return df

def interp(measure: str, t_data: pd.Series, hazard_data: pd.DataFrame, duration=60, step=28/30, k=3) -> pd.DataFrame:
    t_target = np.arange(0, duration, step)
    t_step_start = np.arange(0, len(t_target), 1)
    t_step_end = t_step_start + 1
    df = pd.DataFrame({'t_start': t_step_start, 't_end': t_step_end})
    for i in range(1000):
        h_data = hazard_data[f'draw_{i}']
        h = interpolate.UnivariateSpline(t_data[1:], h_data[1:], k=k) # cubic
        val = h(t_target)
        df[f'{measure}_draw_{i}'] = val
    return df

def get_results(input_dir: str, survival_outcome_type: str, line: str) -> pd.DataFrame:
    if survival_outcome_type == 'overall_survival':
        s_name = 'Survival probability, S(t)'
        measure = 'mortality'
    else:
        s_name = 'Progression-free probability, P(t)'
        measure = 'incidence'
    data = pd.read_excel(os.path.join(input_dir, f'/{survival_outcome_type}_by_time.xlsx'),
                         sheet_name = line,
                         engine = 'openpyxl')
    t = data['Time in months, t']
    num = len(t)
    N = data['Number at risk, N(t)']
    S = data[s_name]
    D = calc_events(N, S)
    Var_S = calc_variance(S, D, N, num)
    S_draws = get_S_draws(S, Var_S, num)
    H_draws = get_H_draws(N, S_draws)
    res = interp(measure, t, H_draws)
    return res


if __name__ == '__main__':
	input_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/treatment_landscape/braunlin_et_al_2020'
	output_dir = '/home/j/Project/simulation_science/multiple_myeloma/data/cause_model_input'
	for line in ['First-line', 'Second-line', 'Third-line', 'Fourth-line', 'Fifth-line']:
		df_mortality = get_results('overall_survival', line)
		df_mortality.to_csv(os.path.join(output_dir, f'/mortality {line}.csv'), index=False)
		df_incidence = get_results('progression_free_survival', line)
		df_incidence.to_csv(os.path.join(output_dir, f'/incidence {line}.csv'), index=False)
