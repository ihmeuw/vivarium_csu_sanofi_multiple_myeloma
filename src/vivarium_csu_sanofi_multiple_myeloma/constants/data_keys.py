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


class __MultipleMyeloma(NamedTuple):

    # Keys that will be loaded into the artifact. must have a colon type declaration
    PREVALENCE: TargetString = TargetString('cause.multiple_myeloma.prevalence')
    INCIDENCE_RATE: TargetString = TargetString('cause.multiple_myeloma.incidence_rate')
    DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma.disability_weight')
    CSMR: TargetString = TargetString('cause.multiple_myeloma.cause_specific_mortality_rate')
    RESTRICTIONS: TargetString = TargetString('cause.multiple_myeloma.restrictions')

    LINE_1_PREVALENCE: TargetString = TargetString('cause.multiple_myeloma_1.prevalence')
    LINE_1_DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma_1.disability_weight')
    LINE_1_EMR: TargetString = TargetString('cause.multiple_myeloma_1.excess_mortality_rate')

    LINE_2_PREVALENCE: TargetString = TargetString('cause.multiple_myeloma_2.prevalence')
    LINE_2_DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma_2.disability_weight')
    LINE_2_EMR: TargetString = TargetString('cause.multiple_myeloma_2.excess_mortality_rate')

    LINE_3_PREVALENCE: TargetString = TargetString('cause.multiple_myeloma_3.prevalence')
    LINE_3_DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma_3.disability_weight')
    LINE_3_EMR: TargetString = TargetString('cause.multiple_myeloma_3.excess_mortality_rate')

    LINE_4_PREVALENCE: TargetString = TargetString('cause.multiple_myeloma_4.prevalence')
    LINE_4_DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma_4.disability_weight')
    LINE_4_EMR: TargetString = TargetString('cause.multiple_myeloma_4.excess_mortality_rate')

    LINE_5_PREVALENCE: TargetString = TargetString('cause.multiple_myeloma_5.prevalence')
    LINE_5_DISABILITY_WEIGHT: TargetString = TargetString('cause.multiple_myeloma_5.disability_weight')
    LINE_5_EMR: TargetString = TargetString('cause.multiple_myeloma_5.excess_mortality_rate')

    GBD_CSMR: TargetString = TargetString('cause.multiple_myeloma_gbd.cause_specific_mortality_rate')
    SUSCEPTIBLE_EMR: TargetString = TargetString('cause.susceptible_to_multiple_myeloma.excess_mortality_rate')

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
