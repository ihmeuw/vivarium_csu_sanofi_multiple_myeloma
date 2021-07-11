from pathlib import Path

import vivarium_csu_sanofi_multiple_myeloma
from vivarium_csu_sanofi_multiple_myeloma.constants import metadata

BASE_DIR = Path(vivarium_csu_sanofi_multiple_myeloma.__file__).resolve().parent

ARTIFACT_ROOT = Path(f"/share/costeffectiveness/artifacts/{metadata.PROJECT_NAME}/")
MODEL_SPEC_DIR = BASE_DIR / 'model_specifications'
RESULTS_ROOT = Path(f'/share/costeffectiveness/results/{metadata.PROJECT_NAME}/')

CSV_RAW_DATA_ROOT = BASE_DIR / 'data' / 'raw_data'
MORTALITY_FIRST_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality First-line.csv'
MORTALITY_SECOND_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality Second-line.csv'
MORTALITY_THIRD_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality Third-line.csv'
MORTALITY_FOURTH_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality Fourth-line.csv'
MORTALITY_FIFTH_LINE_PATH = CSV_RAW_DATA_ROOT / 'mortality Fifth-line.csv'

INCIDENCE_FIRST_LINE_PATH = CSV_RAW_DATA_ROOT / 'incidence First-line.csv'
INCIDENCE_SECOND_LINE_PATH = CSV_RAW_DATA_ROOT / 'incidence Second-line.csv'
INCIDENCE_THIRD_LINE_PATH = CSV_RAW_DATA_ROOT / 'incidence Third-line.csv'
INCIDENCE_FOURTH_LINE_PATH = CSV_RAW_DATA_ROOT / 'incidence Fourth-line.csv'
