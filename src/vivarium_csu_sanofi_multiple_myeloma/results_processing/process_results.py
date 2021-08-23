from pathlib import Path
from typing import NamedTuple, List

import pandas as pd
import yaml

from vivarium_csu_sanofi_multiple_myeloma.constants import results
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import RISKS


SCENARIO_COLUMN = 'scenario'
HAZARD_RATE_SOURCE_COLUMN = 'hazard_rate_source'
RACE_IMPACT_SCENARIO_COLUMN = 'race_impact_scenario'
GROUPBY_COLUMNS = [
    results.INPUT_DRAW_COLUMN,
    SCENARIO_COLUMN,
    HAZARD_RATE_SOURCE_COLUMN,
    RACE_IMPACT_SCENARIO_COLUMN
]
OUTPUT_COLUMN_SORT_ORDER = [
    'age_group',
    'sex',
    'year',
    'risk',
    'cause',
    'measure',
    'input_draw'
]


def make_measure_data(data):
    measure_data = MeasureData(
        population=get_population_data(data),
        person_time=get_measure_data(data, 'person_time', stratified_by_treatment=True, stratified_by_risks=True),
        ylls=get_by_cause_measure_data(data, 'ylls', stratified_by_treatment=True, stratified_by_risks=True),
        ylds=get_by_cause_measure_data(data, 'ylds'),
        deaths=get_by_cause_measure_data(data, 'deaths', stratified_by_treatment=True, stratified_by_risks=True),
        state_person_time=get_state_person_time_measure_data(data, 'state_person_time', stratified_by_treatment=True,
                                                             stratified_by_risks=True),
        transition_count=get_transition_count_measure_data(data, 'transition_count', stratified_by_treatment=True,
                                                           stratified_by_risks=True),
        treatment_count=get_treatment_count_measure_data(data, 'treatment_count'),
        survival=get_survival_measure_data(data),
        registry=get_registry_measure_data(data),
    )
    return measure_data


class MeasureData(NamedTuple):
    population: pd.DataFrame
    person_time: pd.DataFrame
    ylls: pd.DataFrame
    ylds: pd.DataFrame
    deaths: pd.DataFrame
    state_person_time: pd.DataFrame
    transition_count: pd.DataFrame
    treatment_count: pd.DataFrame
    survival: pd.DataFrame
    registry: pd.DataFrame

    def dump(self, output_dir: Path):
        for key, df in self._asdict().items():
            # FIXME: HDFs fail to build due to `OverflowError: value too large to convert to int`, but not required
            # df.to_hdf(output_dir / f'{key}.hdf', key=key)
            df.to_csv(output_dir / f'{key}.csv')


def read_data(path: Path, single_run: bool) -> (pd.DataFrame, List[str]):
    data = pd.read_hdf(path)
    # noinspection PyUnresolvedReferences
    data = (data
            .drop(columns=data.columns.intersection(results.THROWAWAY_COLUMNS))
            .reset_index(drop=True)
            .rename(columns={results.OUTPUT_SCENARIO_COLUMN: SCENARIO_COLUMN})
            .rename(columns={results.HAZARD_RATE_SOURCE_COLUMN: HAZARD_RATE_SOURCE_COLUMN})
            .rename(columns={results.RACE_IMPACT_SCENARIO_COLUMN: RACE_IMPACT_SCENARIO_COLUMN})
            )
    if single_run:
        data[results.INPUT_DRAW_COLUMN] = 0
        data[results.RANDOM_SEED_COLUMN] = 0
        data[SCENARIO_COLUMN] = 'baseline'
        data[HAZARD_RATE_SOURCE_COLUMN] = 'population'
        data[RACE_IMPACT_SCENARIO_COLUMN] = 'commpass_registry'
        keyspace = {results.INPUT_DRAW_COLUMN: [0],
                    results.RANDOM_SEED_COLUMN: [0],
                    results.OUTPUT_SCENARIO_COLUMN: ['baseline'],
                    results.HAZARD_RATE_SOURCE_COLUMN: ['population'],
                    results.RACE_IMPACT_SCENARIO_COLUMN: ['commpass_registry']}
    else:
        data[results.INPUT_DRAW_COLUMN] = data[results.INPUT_DRAW_COLUMN].astype(int)
        data[results.RANDOM_SEED_COLUMN] = data[results.RANDOM_SEED_COLUMN].astype(int)
        with (path.parent / 'keyspace.yaml').open() as f:
            keyspace = yaml.full_load(f)
    return data, keyspace


def filter_out_incomplete(data, keyspace):
    output = []
    for draw in keyspace[results.INPUT_DRAW_COLUMN]:
        # For each draw, gather all random seeds completed for all scenarios.
        random_seeds = set(keyspace[results.RANDOM_SEED_COLUMN])
        draw_data = data.loc[data[results.INPUT_DRAW_COLUMN] == draw]
        for scenario in keyspace[results.OUTPUT_SCENARIO_COLUMN]:
            seeds_in_data = draw_data.loc[data[SCENARIO_COLUMN] == scenario,
                                          results.RANDOM_SEED_COLUMN].unique()
            random_seeds = random_seeds.intersection(seeds_in_data)
        draw_data = draw_data.loc[draw_data[results.RANDOM_SEED_COLUMN].isin(random_seeds)]
        output.append(draw_data)
    return pd.concat(output, ignore_index=True).reset_index(drop=True)


def aggregate_over_seed(data):
    non_count_columns = []
    for non_count_template in results.NON_COUNT_TEMPLATES:
        non_count_columns += results.RESULT_COLUMNS(non_count_template)
    count_columns = [c for c in data.columns if c not in non_count_columns + GROUPBY_COLUMNS]

    # non_count_data = data[non_count_columns + GROUPBY_COLUMNS].groupby(GROUPBY_COLUMNS).mean()
    count_data = data[count_columns + GROUPBY_COLUMNS].groupby(GROUPBY_COLUMNS).sum()
    return pd.concat([
        count_data,
        # non_count_data
    ], axis=1).reset_index()


def pivot_data(data):
    return (data
            .set_index(GROUPBY_COLUMNS)
            .stack()
            .reset_index()
            .rename(columns={f'level_{len(GROUPBY_COLUMNS)}': 'process', 0: 'value'}))


def sort_data(data):
    sort_order = [c for c in OUTPUT_COLUMN_SORT_ORDER if c in data.columns]
    other_cols = [c for c in data.columns if c not in sort_order]
    data = data[sort_order + other_cols].sort_values(sort_order)
    return data.reset_index(drop=True)


def split_processing_column(data, stratified_by_treatment=False, stratified_by_risks=False):
    # TODO find a better way to do this:
    #       FutureWarning: Columnar iteration over characters will be deprecated in future releases.
    if stratified_by_risks:
        data['process'], data['renal_function_at_diagnosis'] = data.process.str.split(
            '_renal_function_at_diagnosis_').str
        data['process'], data['race_and_cytogenetic_risk_at_diagnosis'] = data.process.str.split(
            '_race_and_cytogenetic_risk_at_diagnosis_').str
    if stratified_by_treatment:
        data['process'], data['retreated'] = data.process.str.split('_retreated_').str
        data['process'], data['treatment'] = data.process.str.split('_treatment_state_').str
    data['process'], data['age'] = data.process.str.split('_in_age_group_').str
    data['process'], data['sex'] = data.process.str.split('_among_').str
    data['year'] = data.process.str.split('_in_').str[-1]
    data['measure'] = data.process.str.split('_in_').str[:-1].apply(lambda x: '_in_'.join(x))
    return data.drop(columns='process')


def new_split_processing_column(data, stratified_by_treatment=False, stratified_by_risks=False):
    out = {'measure': [], 'year': [], 'sex': [], 'age': []}
    if stratified_by_treatment:
        out['treatment'] = []
        out['retreated'] = []
    if stratified_by_risks:
        out['renal_function_at_diagnosis'] = []
        out['race_and_cytogenetic_risk_at_diagnosis'] = []
    out['value'] = []

    for k, v in data.iterrows():
        pass

    # return data.drop(columns='process')


def get_population_data(data):
    total_pop = pivot_data(data[[results.TOTAL_POPULATION_COLUMN]
                                + results.RESULT_COLUMNS('population')
                                + GROUPBY_COLUMNS])
    total_pop = total_pop.rename(columns={'process': 'measure'})
    return sort_data(total_pop)


def get_measure_data(data, measure, stratified_by_treatment=False, stratified_by_risks=False):
    data = pivot_data(data[results.RESULT_COLUMNS(measure) + GROUPBY_COLUMNS])
    data = split_processing_column(data, stratified_by_treatment, stratified_by_risks)
    return sort_data(data)


def get_by_cause_measure_data(data, measure, stratified_by_treatment=False, stratified_by_risks=False):
    data = get_measure_data(data, measure, stratified_by_treatment, stratified_by_risks)
    data['measure'], data['cause'] = data.measure.str.split('_due_to_').str
    return sort_data(data)


def get_state_person_time_measure_data(data, measure, stratified_by_treatment, stratified_by_risks):
    data = get_measure_data(data, measure, stratified_by_treatment, stratified_by_risks)
    data['measure'], data['cause'] = 'state_person_time', data.measure.str.split('_person_time').str[0]
    return sort_data(data)


def get_transition_count_measure_data(data, measure, stratified_by_treatment, stratified_by_risks):
    # Oops, edge case.
    data = data.drop(columns=[c for c in data.columns if 'event_count' in c and '2026' in c])
    data = get_measure_data(data, measure, stratified_by_treatment, stratified_by_risks)
    return sort_data(data)


def get_treatment_count_measure_data(data, measure):
    data = pivot_data(data[results.RESULT_COLUMNS(measure) + GROUPBY_COLUMNS])
    data['process'], data['year'] = data.process.str.split('_year_').str
    data['process'], data['treatment'] = data.process.str.split('_treatment_').str
    data['process'], data['treatment_line'] = data.process.str.split('line_').str
    data = data.drop(columns='process')
    return sort_data(data)


def get_survival_measure_data(data):
    data = pivot_data(data[results.RESULT_COLUMNS('survival_alive') + results.RESULT_COLUMNS('survival_other') + GROUPBY_COLUMNS])
    for s in reversed(RISKS):
        data['process'], data[s] = data.process.str.split('_' + s + '_').str
    data['measure'], data['process'] = data.process.str.split('_period_').str
    data['period'], data['treatment_line'] = data.process.str.split('_line_').str

    data = data.drop(columns='process')
    return sort_data(data)


def get_registry_measure_data(data):
    data = pivot_data(data[results.RESULT_COLUMNS('registry_status') + GROUPBY_COLUMNS])
    data['process'], data['renal_function_at_diagnosis'] = data.process.str.split('_renal_function_at_diagnosis_').str
    data['process'], data['race_and_cytogenetic_risk_at_diagnosis'] = data.process.str.split('_race_and_cytogenetic_risk_at_diagnosis_').str
    data['process'], data['age'] = data.process.str.split('_in_age_group_').str
    data['process'], data['sex'] = data.process.str.split('_among_').str
    data['year'] = data.process.str.split('_in_').str[-1]
    data['measure'] = data.process.str.split('_in_').str[:-1].apply(lambda x: '_in_'.join(x))
    data = data.drop(columns='process')
    return sort_data(data)
