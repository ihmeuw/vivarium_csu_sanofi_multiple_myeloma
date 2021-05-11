from typing import NamedTuple

####################
# Project metadata #
####################

PROJECT_NAME = 'vivarium_csu_sanofi_multiple_myeloma'
CLUSTER_PROJECT = 'proj_csu'

CLUSTER_QUEUE = 'all.q'
MAKE_ARTIFACT_MEM = '10G'
MAKE_ARTIFACT_CPU = '1'
MAKE_ARTIFACT_RUNTIME = '3:00:00'
MAKE_ARTIFACT_SLEEP = 10

LOCATIONS = [
    'United States of America'
]


class __Scenarios(NamedTuple):
    baseline: str = 'baseline'
    # TODO - add scenarios here


SCENARIOS = __Scenarios()