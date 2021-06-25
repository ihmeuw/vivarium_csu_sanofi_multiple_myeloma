from typing import NamedTuple

from vivarium_csu_sanofi_multiple_myeloma.utilities import LogNormalHazardRate


# Progression hazard ratios
class __ProgressionHazardRate(NamedTuple):
    FIRST_LINE_ISA: LogNormalHazardRate = LogNormalHazardRate("first_line_isa", 0.429, 0.495) # lower = 0.368, upper = 0.495
    FIRST_LINE_DARA: LogNormalHazardRate = LogNormalHazardRate("first_line_dara", 0.429, 0.495) # lower = 0.368, upper = 0.495
    FIRST_LINE_RESIDUAL: LogNormalHazardRate = LogNormalHazardRate("first_line_residual", 1.00581, 1.0064)  # lower = 1.0051, upper = 1.0064

    LATER_LINE_ISA: LogNormalHazardRate = LogNormalHazardRate("later_line_isa", 0.530,  0.803)  # lower = 0.356, upper = 0.803
    LATER_LINE_ISA_RETREAT: LogNormalHazardRate = LogNormalHazardRate("later_line_isa_retreat", 0.765, 0.902)  # lower = 0.678, upper = 0.902
    LATER_LINE_DARA: LogNormalHazardRate = LogNormalHazardRate("later_line_dara", 0.217, 0.231)  # lower = 0.203, upper = 0.231
    LATER_LINE_DARA_RETREAT: LogNormalHazardRate = LogNormalHazardRate("later_line_dara_retreat", 0.609, 0.616) # lower = 0.601, upper = 0.616
    LATER_LINE_RESIDUAL: LogNormalHazardRate = LogNormalHazardRate("later_line_residual", 1.331, 1.337) # lower = 1.324, upper = 1.337

    @property
    def name(self):
        return 'progression_hazard_rate'

    @property
    def log_name(self):
        return 'progression hazard rate'


PROGRESSION_HAZARD_RATE = __ProgressionHazardRate()


# Mortality hazard ratios
class __MortalityHazardRate(NamedTuple):

    FIRST_LINE_ISA: LogNormalHazardRate = LogNormalHazardRate("first_line_isa", 0.760, 0.895)  # LOWER = 0.645, UPPER = 0.895
    FIRST_LINE_DARA: LogNormalHazardRate = LogNormalHazardRate("first_line_dara", 0.760, 0.895)  # LOWER = 0.645, UPPER = 0.895
    FIRST_LINE_RESIDUAL: LogNormalHazardRate = LogNormalHazardRate("first_line_residual", 1.0024, 1.0036)  # LOWER = 1.0011, UPPER = 1.0036

    LATER_LINE_ISA: LogNormalHazardRate = LogNormalHazardRate("later_line_isa", 1.116, 1.185)  # LOWER = 1.044, UPPER = 1.185
    LATER_LINE_ISA_RETREAT: LogNormalHazardRate = LogNormalHazardRate("later_line_isa_retreat", 1.232, 1.370)  # LOWER = 1.088, UPPER = 1.370
    LATER_LINE_DARA: LogNormalHazardRate = LogNormalHazardRate("later_line_dara", 0.572, 0.594) # LOWER = 0.551, UPPER = 0.594
    LATER_LINE_DARA_RETREAT: LogNormalHazardRate = LogNormalHazardRate("later_line_dara_retreat", 0.786, 0.797) # LOWER = 0.776, UPPER = 0.797
    LATER_LINE_RESIDUAL: LogNormalHazardRate = LogNormalHazardRate("later_line_residual", 1.181, 1.190)  # LOWER = 1.171, UPPER = 1.190

    @property
    def name(self):
        return 'mortality_hazard_rate'

    @property
    def log_name(self):
        return 'mortality hazard rate'


MORTALITY_HAZARD_RATE = __MortalityHazardRate()

PROBABILITY_RETREAT = 0.15