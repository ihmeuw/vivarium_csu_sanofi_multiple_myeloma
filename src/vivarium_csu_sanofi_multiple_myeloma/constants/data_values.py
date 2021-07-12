from vivarium_csu_sanofi_multiple_myeloma.constants import models


# Population-based Hazard Ratio PFS Distributions
PFS_HR = {
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.932, 0.647, 1.365),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.daratumumab, False): (0.932, 0.647, 1.365),
    (models.MULTIPLE_MYELOMA_1_STATE_NAME, models.TREATMENTS.residual, False): (1.002, 0.989, 1.018),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.283, 0.878, 1.178),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.isatuximab, True): (1.632, 0.905, 2.733),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.146, 1.000, 1.318),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.333, 0.995, 1.702),
    (models.MULTIPLE_MYELOMA_2_STATE_NAME, models.TREATMENTS.residual, False): (0.962, 0.920, 1.000),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.405, 0.924, 2.020),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.isatuximab, True): (1.883, 0.974, 3.100),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.133, 0.977, 1.296),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.345, 0.993, 1.747),
    (models.MULTIPLE_MYELOMA_3_STATE_NAME, models.TREATMENTS.residual, False): (0.930, 0.852, 1.001),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.736, 0.394, 1.265),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, True): (0.878, 0.653, 1.583),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.098, 0.877, 1.327),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.275, 0.981, 1.843),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, False): (0.955, 0.822, 1.081),
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
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, False): (1.627, 0.948, 2.628),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.isatuximab, True): (2.333, 1.031, 4.074),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.217, 0.976, 1.467),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.620, 1.008, 2.210),
    (models.MULTIPLE_MYELOMA_4_STATE_NAME, models.TREATMENTS.residual, False): (0.834, 0.702, 0.969),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.isatuximab, False): (0.592, 0.103, 1.947),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.isatuximab, True): (0.914, 0.493, 2.643),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.daratumumab, False): (1.217, 0.976, 1.467),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.daratumumab, True): (1.427, 0.834, 2.410),
    (models.MULTIPLE_MYELOMA_5_STATE_NAME, models.TREATMENTS.residual, False): (0.952, 0.744, 1.145),
}

PROBABILITY_RETREAT = 0.15
