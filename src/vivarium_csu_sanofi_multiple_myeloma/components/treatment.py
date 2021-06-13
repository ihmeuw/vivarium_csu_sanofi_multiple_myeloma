from typing import NamedTuple, TYPE_CHECKING

import pandas as pd

from vivarium_csu_sanofi_multiple_myeloma.constants import models
from vivarium_csu_sanofi_multiple_myeloma.constants.metadata import SCENARIOS

if TYPE_CHECKING:
    from vivarium.framework.engine import Builder


class __Treatments(NamedTuple):
    isatuxamib: str
    daratumamab: str
    residual: str


TREATMENTS = __Treatments(*__Treatments._fields)

TREATMENT_LINES = pd.Index([
    models.MULTIPLE_MYELOMA_1_STATE_NAME,
    models.MULTIPLE_MYELOMA_2_STATE_NAME,
    models.MULTIPLE_MYELOMA_3_STATE_NAME,
    models.MULTIPLE_MYELOMA_4_STATE_NAME,
    models.MULTIPLE_MYELOMA_5_STATE_NAME,
], name=models.MULTIPLE_MYELOMA_MODEL_NAME)

def make_treatment_coverage(year, scenario):
    coverages = {
        (2021, SCENARIOS.baseline): (
            [0.000, 0.008, 0.013, 0.015, 0.009],
            [0.016, 0.169, 0.256, 0.297, 0.171],
        ),
        (2025, SCENARIOS.baseline): (
            [0.000, 0.100, 0.100, 0.100, 0.100],
            [0.029, 0.306, 0.463, 0.537, 0.309],
        ),
        (2025, SCENARIOS.alternative): (
            [0.100, 0.100, 0.100, 0.100, 0.100],
            [0.029, 0.306, 0.463, 0.537, 0.309],
        )
    }
    coverages[(2021, SCENARIOS.alternative)] = coverages[(2021, SCENARIOS.baseline)]

    coverage_data = coverages[(year, scenario)]
    coverage = pd.DataFrame({
        TREATMENTS.isatuxamib: coverage_data[0],
        TREATMENTS.daratumamab: coverage_data[1],
    }, index=TREATMENT_LINES)
    coverage[TREATMENTS.residual] = 1 - coverage.sum(axis=1)
    return coverage.reset_index()


class MMTreatment:

    configuration_defaults = {
        'mm_treatment_scenario': SCENARIOS.baseline,
    }

    def setup(self, builder: 'Builder') -> None:




