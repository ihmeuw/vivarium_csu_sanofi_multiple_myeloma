import itertools

from vivarium_csu_sanofi_multiple_myeloma.constants import models

#################################
# Results columns and variables #
#################################

TOTAL_POPULATION_COLUMN = 'total_population'
TOTAL_YLDS_COLUMN = 'years_lived_with_disability'
TOTAL_YLLS_COLUMN = 'years_of_life_lost'

# Columns from parallel runs
INPUT_DRAW_COLUMN = 'input_draw'
RANDOM_SEED_COLUMN = 'random_seed'
OUTPUT_SCENARIO_COLUMN = 'mm_treatment_scenario'

STANDARD_COLUMNS = {
    'total_population': TOTAL_POPULATION_COLUMN,
    'total_ylls': TOTAL_YLLS_COLUMN,
    'total_ylds': TOTAL_YLDS_COLUMN,
}

THROWAWAY_COLUMNS = [f'{state}_event_count' for state in models.STATES]

TOTAL_POPULATION_COLUMN_TEMPLATE = 'total_population_{POP_STATE}'
PERSON_TIME_COLUMN_TEMPLATE = 'person_time_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{BOOL}'
DEATH_COLUMN_TEMPLATE = 'death_due_to_{CAUSE_OF_DEATH}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{BOOL}'
YLLS_COLUMN_TEMPLATE = 'ylls_due_to_{CAUSE_OF_DEATH}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{BOOL}'
YLDS_COLUMN_TEMPLATE = 'ylds_due_to_{CAUSE_OF_DISABILITY}_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}'
STATE_PERSON_TIME_COLUMN_TEMPLATE = '{STATE}_person_time_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{BOOL}'
TRANSITION_COUNT_COLUMN_TEMPLATE = '{TRANSITION}_event_count_in_{YEAR}_among_{SEX}_in_age_group_{AGE_GROUP}_treatment_state_{TREATMENT}_retreated_{BOOL}'
TREATMENT_COUNT_COLUMN_TEMPLATE = 'line_{TREATMENT_LINE}_treatment_{TREATMENT}_year_{YEAR}'
SURVIVAL_ALIVE_TEMPLATE = 'alive_at_period_{LEFT_PERIOD}_line_{TREATMENT_LINE}'
SURVIVAL_OTHER_TEMPLATE = '{SURVIVAL_METRIC}_period_{RIGHT_PERIOD}_line_{TREATMENT_LINE}'

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
}

NON_COUNT_TEMPLATES = [
]

POP_STATES = ('living', 'dead', 'tracked', 'untracked')
SEXES = ('male', 'female')
#XXX
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
# TODO - add causes of death
CAUSES_OF_DEATH = (
    'other_causes',
    models.SUSCEPTIBLE_STATE_NAME,
    models.MULTIPLE_MYELOMA_1_STATE_NAME,
    models.MULTIPLE_MYELOMA_2_STATE_NAME,
    models.MULTIPLE_MYELOMA_3_STATE_NAME,
    models.MULTIPLE_MYELOMA_4_STATE_NAME,
    models.MULTIPLE_MYELOMA_5_STATE_NAME,
)
# TODO - add causes of disability
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
    'BOOL': ('True', 'False'),
    'SURVIVAL_METRIC': ('progressed_by', 'died_by', 'sim_end_on'),
    # CAUTION: This needs to change with the time step size.
    'LEFT_PERIOD': ALL_PERIODS[:-1],
    'RIGHT_PERIOD': ALL_PERIODS[1:],
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
