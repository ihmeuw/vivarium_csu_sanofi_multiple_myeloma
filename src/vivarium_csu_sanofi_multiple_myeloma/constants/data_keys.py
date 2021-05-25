from typing import NamedTuple

from vivarium_public_health.utilities import TargetString


#############
# Data Keys #
#############

METADATA_LOCATIONS = 'metadata.locations'


class __Population(NamedTuple):
    LOCATION: str = 'population.location'
    STRUCTURE: str = 'population.structure'
    AGE_BINS: str = 'population.age_bins'
    DEMOGRAPHY: str = 'population.demographic_dimensions'
    TMRLE: str = 'population.theoretical_minimum_risk_life_expectancy'
    ACMR: str = 'cause.all_causes.cause_specific_mortality_rate'

    @property
    def name(self):
        return 'population'

    @property
    def log_name(self):
        return 'population'


POPULATION = __Population()


# TODO - sample key group used to identify keys in model
# For more information see the tutorial:
# https://vivarium-inputs.readthedocs.io/en/latest/tutorials/pulling_data.html#entity-measure-data
class __MultipleMyeloma(NamedTuple):

    # Keys that will be loaded into the artifact. must have a colon type declaration
    PREVALENCE: TargetString = TargetString('cause.multiple_myeloma.prevalence')
    INCIDENCE_RATE: TargetString = TargetString('cause.multiple_myeloma.incidence_rate')
    DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma.disability_weight')
    EMR: TargetString = TargetString('cause.multiple_myeloma.excess_mortality_rate')
    CSMR: TargetString = TargetString('cause.multiple_myeloma.cause_specific_mortality_rate')
    RESTRICTIONS: TargetString = TargetString('cause.multiple_myeloma.restrictions')

    # Useful keys not for the artifact - distinguished by not using the colon type declaration
    # RAW_DISEASE_PREVALENCE = TargetString('sequela.raw_disease.prevalence')
    # RAW_DISEASE_INCIDENCE_RATE = TargetString('sequela.raw_disease.incidence_rate')

    @property
    def name(self):
        return 'multiple_myeloma'

    @property
    def log_name(self):
        return 'multiple myeloma'


MULTIPLE_MYELOMA = __MultipleMyeloma()

MAKE_ARTIFACT_KEY_GROUPS = [
    POPULATION,
    MULTIPLE_MYELOMA
]
