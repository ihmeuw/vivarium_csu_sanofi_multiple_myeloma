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

RISK_PFS_HR = {
    RISK_EXPOSURE_LEVELS.Male: (1.117379506, 1.023631793, 1.206248002),
    RISK_EXPOSURE_LEVELS.Female: (0.862760188, 0.758855373, 0.97236977),
    RISK_EXPOSURE_LEVELS.over_65: (1.168929731, 1.109344393, 1.225450674),
    RISK_EXPOSURE_LEVELS.under_65: (0.690375253, 0.58678021, 0.799586905),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: 1.,
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1.029416953, 1.007427712, 1.051070122),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: 1.,
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.799597005, 0.652084796, 0.949398709),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.366074674, 1.140038711, 1.60581245),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.967734441, 0.946604126, 0.987657089)
}

RISK_OS_HR = {
    RISK_EXPOSURE_LEVELS.Male: (1.255109255, 1.112976053, 1.375112298),
    RISK_EXPOSURE_LEVELS.Female: (0.701726923, 0.561419677, 0.867908693),
    RISK_EXPOSURE_LEVELS.over_65: (1.236326008, 1.157185865, 1.302944895),
    RISK_EXPOSURE_LEVELS.under_65: (0.566847231, 0.44474406, 0.711900129),
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_black: 1.,
    RISK_EXPOSURE_LEVELS.high_cytogenetic_risk_and_non_black: (1.059433242, 1.019805071, 1.089014929),
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_black: 1.,
    RISK_EXPOSURE_LEVELS.low_cytogenetic_risk_and_non_black: (0.595111041, 0.393585797, 0.865077955),
    RISK_EXPOSURE_LEVELS.renal_impaired: (1.788742852, 1.35961846, 2.299354949),
    RISK_EXPOSURE_LEVELS.renal_unimpaired: (0.930480771, 0.885475788, 0.968303487)
}
