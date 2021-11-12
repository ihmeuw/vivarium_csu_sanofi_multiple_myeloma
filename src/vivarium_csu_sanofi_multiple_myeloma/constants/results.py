import itertools

import pandas as pd

from vivarium_csu_sanofi_multiple_myeloma.constants import models
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import RISKS, RISK_LEVEL_MAP, UNDIAGNOSED

#################################
# Results columns and variables #
#################################

TOTAL_POPULATION_COLUMN = 'total_population'
TOTAL_YLDS_COLUMN = 'years_lived_with_disability'
TOTAL_YLLS_COLUMN = 'years_of_life_lost'

# Columns from parallel runs
INPUT_DRAW_COLUMN = 'input_draw'
RANDOM_SEED_COLUMN = 'random_seed'
OUTPUT_SCENARIO_COLUMN = 'mm_scenarios.mm_treatment_scenario'
HAZARD_RATE_SOURCE_COLUMN = 'mm_scenarios.hazard_rate_source'
RACE_IMPACT_SCENARIO_COLUMN = 'mm_scenarios.race_impact_scenario'

STANDARD_COLUMNS = {
    'total_population': TOTAL_POPULATION_COLUMN,
    'total_ylls': TOTAL_YLLS_COLUMN,
    'total_ylds': TOTAL_YLDS_COLUMN,
}

THROWAWAY_COLUMNS = [f'{state}_event_count' for state in models.STATES]

TOTAL_POPULATION_COLUMN_TEMPLATE = 'total_population_{POP_STATE}'
PERSON_TIME_COLUMN_TEMPLATE = 'person_time_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{RETREATED}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS_EXT}'
DEATH_COLUMN_TEMPLATE = 'death_due_to_{CAUSE_OF_DEATH}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{RETREATED}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS_EXT}'
YLLS_COLUMN_TEMPLATE = 'ylls_due_to_{CAUSE_OF_DEATH}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{RETREATED}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS_EXT}'
YLDS_COLUMN_TEMPLATE = 'ylds_due_to_{CAUSE_OF_DISABILITY}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}'
STATE_PERSON_TIME_COLUMN_TEMPLATE = '{STATE}_person_time_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{RETREATED}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS_EXT}'
TRANSITION_COUNT_COLUMN_TEMPLATE = '{TRANSITION}_event_count_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{RETREATED}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS_EXT}'
TREATMENT_COUNT_COLUMN_TEMPLATE = 'line_{TREATMENT_LINE}_treatment_{TREATMENT}_year_{YEAR}'
SURVIVAL_ALIVE_TEMPLATE = 'alive_at_period_{LEFT_PERIOD}_line_{TREATMENT_LINE}' + '_' + '_'.join([f'{s}_{{{s.upper()}}}' for s in RISKS])
SURVIVAL_OTHER_TEMPLATE = '{SURVIVAL_METRIC}_period_{RIGHT_PERIOD}_line_{TREATMENT_LINE}' + '_' + '_'.join([f'{s}_{{{s.upper()}}}' for s in RISKS])
REGISTRY_TEMPLATE = 'registry_status_{REGISTRY_STATUS}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_race_and_cytogenetic_risk_at_diagnosis_{RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS}_renal_function_at_diagnosis_{RENAL_FUNCTION_AT_DIAGNOSIS}'

COLUMN_TEMPLATES = {
    'population': TOTAL_POPULATION_COLUMN_TEMPLATE,
    'person_time': PERSON_TIME_COLUMN_TEMPLATE,
    'deaths': DEATH_COLUMN_TEMPLATE,
    'ylls': YLLS_COLUMN_TEMPLATE,
    'ylds': YLDS_COLUMN_TEMPLATE,
    'state_person_time': STATE_PERSON_TIME_COLUMN_TEMPLATE,
    'transition_count': TRANSITION_COUNT_COLUMN_TEMPLATE,
    'treatment_count': TREATMENT_COUNT_COLUMN_TEMPLATE,
    'survival_alive': SURVIVAL_ALIVE_TEMPLATE,
    'survival_other': SURVIVAL_OTHER_TEMPLATE,
    'registry_status': REGISTRY_TEMPLATE
}

NON_COUNT_TEMPLATES = [
]

POP_STATES = ('living', 'dead', 'tracked', 'untracked')
SEXES = ('male', 'female')
YEARS = tuple(range(2011, 2026))
AGE_GROUPS = (
    '15_to_19',
    '20_to_24',
    '25_to_29',
    '30_to_34',
    '35_to_39',
    '40_to_44',
    '45_to_49',
    '50_to_54',
    '55_to_59',
    '60_to_64',
    '65_to_69',
    '70_to_74',
    '75_to_79',
    '80_to_84',
    '85_to_89',
    '90_to_94',
    '95_plus',
)
CAUSES_OF_DEATH = (
    'other_causes',
    models.SUSCEPTIBLE_STATE_NAME,
    models.MULTIPLE_MYELOMA_1_STATE_NAME,
    models.MULTIPLE_MYELOMA_2_STATE_NAME,
    models.MULTIPLE_MYELOMA_3_STATE_NAME,
    models.MULTIPLE_MYELOMA_4_STATE_NAME,
    models.MULTIPLE_MYELOMA_5_STATE_NAME,
)
CAUSES_OF_DISABILITY = (
    models.MULTIPLE_MYELOMA_1_STATE_NAME,
    models.MULTIPLE_MYELOMA_2_STATE_NAME,
    models.MULTIPLE_MYELOMA_3_STATE_NAME,
    models.MULTIPLE_MYELOMA_4_STATE_NAME,
    models.MULTIPLE_MYELOMA_5_STATE_NAME,
)

# CAUTION: This needs to change with the time step size.
ALL_PERIODS = [i * 28 for i in range(61)] + [28*1000]


TEMPLATE_FIELD_MAP = {
    'POP_STATE': POP_STATES,
    'YEAR': YEARS,
    'SEX': SEXES,
    'AGE_GROUP': AGE_GROUPS,
    'CAUSE_OF_DEATH': CAUSES_OF_DEATH,
    'CAUSE_OF_DISABILITY': CAUSES_OF_DISABILITY,
    'STATE': models.STATES,
    'TRANSITION': models.TRANSITIONS,
    'TREATMENT': models.TREATMENTS,
    'TREATMENT_LINE': list(range(1, 6)),
    'RETREATED': ('True', 'False'),
    'SURVIVAL_METRIC': ('progressed_by', 'died_by', 'sim_end_on'),
    # CAUTION: This needs to change with the time step size.
    'LEFT_PERIOD': ALL_PERIODS[:-1],
    'RIGHT_PERIOD': ALL_PERIODS[1:],
    'SEX_AT_DIAGNOSIS': RISK_LEVEL_MAP[RISKS.sex_at_diagnosis],
    'AGE_AT_DIAGNOSIS': RISK_LEVEL_MAP[RISKS.age_at_diagnosis],
    'RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS': RISK_LEVEL_MAP[
        RISKS.race_and_cytogenetic_risk_at_diagnosis],
    'RENAL_FUNCTION_AT_DIAGNOSIS': RISK_LEVEL_MAP[RISKS.renal_function_at_diagnosis],
    'RACE_AND_CYTOGENETIC_RISK_AT_DIAGNOSIS_EXT': RISK_LEVEL_MAP[
                                                      RISKS.race_and_cytogenetic_risk_at_diagnosis] + [UNDIAGNOSED],
    'RENAL_FUNCTION_AT_DIAGNOSIS_EXT': RISK_LEVEL_MAP[RISKS.renal_function_at_diagnosis] + [UNDIAGNOSED],
    'REGISTRY_STATUS': ['newly_eligible', 'newly_enrolled', 'enrolled']
}


def RESULT_COLUMNS(kind='all'):
    if kind not in COLUMN_TEMPLATES and kind != 'all':
        raise ValueError(f'Unknown result column type {kind}')
    columns = []
    if kind == 'all':
        for k in COLUMN_TEMPLATES:
            columns += RESULT_COLUMNS(k)
        columns = list(STANDARD_COLUMNS.values()) + columns
    else:
        template = COLUMN_TEMPLATES[kind]
        filtered_field_map = {field: values
                              for field, values in TEMPLATE_FIELD_MAP.items() if f'{{{field}}}' in template}
        fields, value_groups = filtered_field_map.keys(), itertools.product(*filtered_field_map.values())
        for value_group in value_groups:
            columns.append(template.format(**{field: value for field, value in zip(fields, value_group)}))
    return columns


def RESULTS_MAP(kind):
    if kind not in COLUMN_TEMPLATES:
        raise ValueError(f'Unknown result column type {kind}')
    columns = []
    template = COLUMN_TEMPLATES[kind]
    filtered_field_map = {field: values
                          for field, values in TEMPLATE_FIELD_MAP.items() if f'{{{field}}}' in template}
    fields, value_groups = list(filtered_field_map.keys()), list(itertools.product(*filtered_field_map.values()))
    for value_group in value_groups:
        columns.append(template.format(**{field: value for field, value in zip(fields, value_group)}))
    df = pd.DataFrame(value_groups, columns=map(lambda x: x.lower(), fields))
    df['key'] = columns
    df['measure'] = kind  # per researcher feedback, this column is useful, even when it's identical for all rows
    return df.set_index('key').sort_index()
