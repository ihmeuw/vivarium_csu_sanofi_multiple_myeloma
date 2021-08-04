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
    baseline: str
    alternative: str


SCENARIOS = __Scenarios(*__Scenarios._fields)


class __HazardRateSources(NamedTuple):
    population: str
    clinical_trial: str


HAZARD_RATE_SOURCES = __HazardRateSources(*__HazardRateSources._fields)


class __RaceImpactScenario(NamedTuple):
    commpass_registry: str
    no_impact: str


RACE_IMPACT_SCENARIO = __RaceImpactScenario(*__RaceImpactScenario._fields)
