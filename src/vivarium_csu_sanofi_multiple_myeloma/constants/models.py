from vivarium_csu_sanofi_multiple_myeloma.constants import data_keys


class TransitionString(str):

    def __new__(cls, value):
        # noinspection PyArgumentList
        obj = str.__new__(cls, value.lower())
        obj.from_state, obj.to_state = value.split('_TO_')
        return obj


###########################
# Disease Model variables #
###########################

MULTIPLE_MYELOMA_MODEL_NAME = data_keys.MULTIPLE_MYELOMA.name
SUSCEPTIBLE_STATE_NAME = f'susceptible_to_{MULTIPLE_MYELOMA_MODEL_NAME}'
MULTIPLE_MYELOMA_1_STATE_NAME = MULTIPLE_MYELOMA_MODEL_NAME + "_1"
MULTIPLE_MYELOMA_2_STATE_NAME = MULTIPLE_MYELOMA_MODEL_NAME + "_2"
MULTIPLE_MYELOMA_3_STATE_NAME = MULTIPLE_MYELOMA_MODEL_NAME + "_3"
MULTIPLE_MYELOMA_4_STATE_NAME = MULTIPLE_MYELOMA_MODEL_NAME + "_4"
MULTIPLE_MYELOMA_5_STATE_NAME = MULTIPLE_MYELOMA_MODEL_NAME + "_5"
MULTIPLE_MYELOMA_MODEL_STATES = (SUSCEPTIBLE_STATE_NAME, MULTIPLE_MYELOMA_1_STATE_NAME, MULTIPLE_MYELOMA_2_STATE_NAME,
                                 MULTIPLE_MYELOMA_3_STATE_NAME, MULTIPLE_MYELOMA_4_STATE_NAME,
                                 MULTIPLE_MYELOMA_5_STATE_NAME)
MULTIPLE_MYELOMA_MODEL_TRANSITIONS = (
    TransitionString(f'{SUSCEPTIBLE_STATE_NAME}_TO_{MULTIPLE_MYELOMA_1_STATE_NAME}'),
    TransitionString(f'{MULTIPLE_MYELOMA_1_STATE_NAME}_TO_{MULTIPLE_MYELOMA_2_STATE_NAME}'),
    TransitionString(f'{MULTIPLE_MYELOMA_2_STATE_NAME}_TO_{MULTIPLE_MYELOMA_3_STATE_NAME}'),
    TransitionString(f'{MULTIPLE_MYELOMA_3_STATE_NAME}_TO_{MULTIPLE_MYELOMA_4_STATE_NAME}'),
    TransitionString(f'{MULTIPLE_MYELOMA_4_STATE_NAME}_TO_{MULTIPLE_MYELOMA_5_STATE_NAME}'),
)

STATE_MACHINE_MAP = {
    MULTIPLE_MYELOMA_MODEL_NAME: {
        'states': MULTIPLE_MYELOMA_MODEL_STATES,
        'transitions': MULTIPLE_MYELOMA_MODEL_TRANSITIONS,
    },
}


STATES = tuple(state for model in STATE_MACHINE_MAP.values() for state in model['states'])
TRANSITIONS = tuple(state for model in STATE_MACHINE_MAP.values() for state in model['transitions'])

