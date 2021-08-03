import itertools
from typing import NamedTuple

from vivarium_csu_sanofi_multiple_myeloma.constants import models


# Population-based Hazard Ratio PFS Distributions
PFS_HR = {
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.932, 0.647, 1.365),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.932, 0.647, 1.365),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.residual, False): (1.002, 0.989, 1.018),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.283, 0.878, 1.718),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, True): (1.632, 0.905, 2.733),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.146, 1.000, 1.318),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.333, 0.995, 1.702),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, False): (0.962, 0.920, 1.000),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.405, 0.924, 2.020),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, True): (1.883, 0.974, 3.100),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.133, 0.977, 1.296),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.345, 0.993, 1.747),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.residual, False): (0.930, 0.852, 1.001),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.residual, True): (0.930, 0.852, 1.001),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.736, 0.394, 1.265),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, True): (0.878, 0.653, 1.583),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.098, 0.877, 1.327),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.275, 0.981, 1.843),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, False): (0.955, 0.822, 1.081),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, True): (0.955, 0.822, 1.081),
}

# Population-based Hazard Ratio OS Distributions
OS_HR = {
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.971, 0.627, 1.488),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.971, 0.627, 1.488),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.residual, False): (1.001, 0.986, 1.011),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.517, 0.939, 2.349),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, True): (2.085, 0.946, 3.634),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.225, 1.035, 1.443),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.502, 1.051, 1.944),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, False): (0.941, 0.887, 0.987),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.453, 0.896, 2.407),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, True): (2.008, 0.975, 3.790),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.265, 1.078, 1.457),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.685, 1.231, 2.152),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.residual, False): (0.865, 0.773, 0.951),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.residual, True): (0.865, 0.773, 0.951),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.627, 0.948, 2.628),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, True): (2.333, 1.031, 4.074),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.217, 0.976, 1.467),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.620, 1.008, 2.210),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, False): (0.834, 0.702, 0.969),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, True): (0.834, 0.702, 0.969),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.592, 0.103, 1.947),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.isatuximab, True): (0.914, 0.493, 2.643),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.217, 0.976, 1.467),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.427, 0.834, 2.410),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.residual, False): (0.952, 0.744, 1.145),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.residual, True): (0.952, 0.744, 1.145),
}

# Randomized Clinical Trial Progression Free Survival Hazard Ratios
RCT_PFS_HR = {
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.506, 0.402, 0.620),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.506, 0.402, 0.620),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.residual, False): (1.015, 1.011, 1.018),
    # 2nd-4th lines
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.814, 0.593, 1.056),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, True): (0.927, 0.714, 1.077),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.949, 0.581, 1.681),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, True): (0.987, 0.892, 1.207),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, False): (1.023, 0.627, 1.272),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, True): (1.023, 0.627, 1.272),
}

# Fill in 3rd and 4th lines using the same values as in the 2nd line
LATER_LINES = [models.MULTIPLE_MYELOMA_3_STATE_NAME, models.MULTIPLE_MYELOMA_4_STATE_NAME]
ACTUAL_TREATMENTS = [t for t in models.TREATMENTS if t != models.TREATMENTS.not_treated]
for line, treatment, retreated in itertools.product(LATER_LINES, ACTUAL_TREATMENTS, [True, False]):
    RCT_PFS_HR[(line, treatment, retreated)] = RCT_PFS_HR[(models.MULTIPLE_MYELOMA_2_STATE_NAME, treatment, retreated)]

# Randomized Clinical Trial Overall Survival Hazard Ratios
RCT_OS_HR = {
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.760, 0.645, 0.895),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.760, 0.645, 0.895),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.residual, False): (1.007, 1.003, 1.010),
    # 2nd-5th lines
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.031, 0.960, 1.105),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, True): (1.056, 0.928, 1.181),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.031, 0.960, 1.105),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.056, 0.928, 1.181),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, False): (0.984, 0.929, 1.020),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, True): (0.984, 0.929, 1.020),
}
# Fill in 3rd-5th lines using the same values as in the 2nd line
LATER_LINES += [models.MULTIPLE_MYELOMA_5_STATE_NAME]
for line, treatment, retreated in itertools.product(LATER_LINES, ACTUAL_TREATMENTS, [True, False]):
    RCT_OS_HR[(line, treatment, retreated)] = RCT_OS_HR[(models.MULTIPLE_MYELOMA_2_STATE_NAME, treatment, retreated)]


PROBABILITY_RETREAT = 0.15
REGISTRY_ENROLL_PROBABILITY = 0.05


class __Risks(NamedTuple):
    sex_at_diagnosis: str
    age_at_diagnosis: str
    race_and_cytogenetic_risk_at_diagnosis: str
    renal_function_at_diagnosis: str


RISKS = __Risks(*__Risks._fields)


class __RiskExposureLevels(NamedTuple):
    Male: str
    Female: str
    over_65: str
    under_65: str
    high_cytogenetic_risk_and_black: str
    high_cytogenetic_risk_and_non_black: str
    low_cytogenetic_risk_and_black: str
    low_cytogenetic_risk_and_non_black: str
    renal_impaired: str
    renal_unimpaired: str


RISK_EXPOSURE_LEVELS = __RiskExposureLevels(*__RiskExposureLevels._fields)

RISK_LEVEL_MAP = {
    RISKS.sex_at_diagnosis: [RISK_EXPOSURE_LEVELS.Male, RISK_EXPOSURE_LEVELS.Female],
    RISKS.age_at_diagnosis: [RISK_EXPOSURE_LEVELS.over_65, RISK_EXPOSURE_LEVELS.under_65],
    RISKS.race_and_cytogenetic_risk_at_diagnosis: [RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black,
                                                   RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black,
                                                   RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black,
                                                   RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black],
    RISKS.renal_function_at_diagnosis: [RISK_EXPOSURE_LEVELS.renal_impaired, RISK_EXPOSURE_LEVELS.renal_unimpaired]
}

UNDIAGNOSED = 'undiagnosed'

# Risk Exposure Distributions for Risk Effects Calculation
SEX_RISK_EXPOSURE = 0.539  # exposed: male
AGE_RISK_EXPOSURE = 0.647  # exposed: 65+
RACE_RISK_EXPOSURE = 0.177  # exposed: black
CYTOGENETIC_RISK_EXPOSURE = 0.872  # exposed: high
RENAL_RISK_EXPOSURE = 0.081  # exposed: impaired

RACE_AND_CYTO_EXPOSURES = {
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: RACE_RISK_EXPOSURE * CYTOGENETIC_RISK_EXPOSURE,
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1 - RACE_RISK_EXPOSURE) * CYTOGENETIC_RISK_EXPOSURE,
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: RACE_RISK_EXPOSURE * (1 - CYTOGENETIC_RISK_EXPOSURE),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (1 - RACE_RISK_EXPOSURE) * (1 - CYTOGENETIC_RISK_EXPOSURE)
}

RENAL_RISK_EXPOSURE = {
    RISK_EXPOSURE_LEVELS.renal_impaired: RENAL_RISK_EXPOSURE,
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (1 - RENAL_RISK_EXPOSURE)
}

# Hazard ratios from the CoMMpass registry (2a)
RISK_PFS_HR_2A = {
    RISK_EXPOSURE_LEVELS.Male: (1.12, 1.02, 1.21),
    RISK_EXPOSURE_LEVELS.Female: (0.86, 0.76, 0.97),
    RISK_EXPOSURE_LEVELS.over_65: (1.17, 1.11, 1.23),
    RISK_EXPOSURE_LEVELS.under_65: (0.69, 0.59, 0.8),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: (1.31, 1.06, 1.59),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1.31, 1.06, 1.59),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: (1.31, 1.06, 1.59),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.75, 0.62, 0.89),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.37, 1.14, 1.61),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.97, 0.95, 0.99)
}

RISK_OS_HR_2A = {
    RISK_EXPOSURE_LEVELS.Male: (1.26, 1.11, 1.38),
    RISK_EXPOSURE_LEVELS.Female: (0.7, 0.56, 0.87),
    RISK_EXPOSURE_LEVELS.over_65: (1.24, 1.16, 1.3),
    RISK_EXPOSURE_LEVELS.under_65: (0.57, 0.44, 0.71),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: (1.52, 1.15, 1.91),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (0.94, 0.85, 1.03),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: (1.52, 1.15, 1.91),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.53, 0.35, 0.79),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.79, 1.36, 2.3),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.93, 0.89, 0.97)
}

# Assumption of no impact of race on multiple myeloma survival outcomes independent of age (2b)
# i.e., race_impact_scenario == no_impact
RISK_PFS_HR_2B = {
    RISK_EXPOSURE_LEVELS.Male: (1.12, 1.02, 1.21),
    RISK_EXPOSURE_LEVELS.Female: (0.86, 0.76, 0.97),
    RISK_EXPOSURE_LEVELS.over_65: (1.17, 1.11, 1.23),
    RISK_EXPOSURE_LEVELS.under_65: (0.69, 0.59, 0.8),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: (1.06, 1.03, 1.08),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1.06, 1.03, 1.08),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: (0.63, 0.49, 0.79),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.63, 0.49, 0.79),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.37, 1.14, 1.61),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.97, 0.95, 0.99)
}

RISK_OS_HR_2B = {
    RISK_EXPOSURE_LEVELS.Male: (1.26, 1.11, 1.38),
    RISK_EXPOSURE_LEVELS.Female: (0.7, 0.56, 0.87),
    RISK_EXPOSURE_LEVELS.over_65: (1.24, 1.16, 1.3),
    RISK_EXPOSURE_LEVELS.under_65: (0.57, 0.44, 0.71),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: (1.05, 1.02, 1.07),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1.05, 1.02, 1.07),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: (0.66, 0.52, 0.84),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.66, 0.52, 0.84),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.79, 1.36, 2.3),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.93, 0.89, 0.97)
}
