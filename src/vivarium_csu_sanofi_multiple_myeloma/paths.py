from pathlib import Path

import vivarium_csu_sanofi_multiple_myeloma
from vivarium_csu_sanofi_multiple_myeloma.constants import metadata

BASE_DIR = Path(vivarium_csu_sanofi_multiple_myeloma.__file__).resolve().parent

ARTIFACT_ROOT = Path(f"/share/costeffectiveness/artifacts/{metadata.PROJECT_NAME}/")
MODEL_SPEC_DIR = BASE_DIR / 'model_specifications'
RESULTS_ROOT = Path(f'/share/costeffectiveness/results/{metadata.PROJECT_NAME}/')

CSV_RAW_DATA_ROOT = BASE_DIR / 'data' / 'raw_data'
MORTALITY_FIRST_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality First-line.csv'
