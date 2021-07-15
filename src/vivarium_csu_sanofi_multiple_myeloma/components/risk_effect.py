import itertools
from collections import Counter

from typing import Callable, Dict, Iterable, List, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd


from vivarium.framework.engine import Builder
from vivarium.framework.event import Event
from vivarium.framework.randomness import get_hash
from vivarium.framework.population import SimulantData

from vivarium_public_health.metrics.utilities import (get_deaths, get_person_time, get_state_person_time,
                                                      get_transition_count,
                                                      get_years_lived_with_disability, get_years_of_life_lost,
                                                      TransitionString)

from vivarium_csu_sanofi_multiple_myeloma.constants import models, results
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import (RACE_AND_CYTO_EXPOSURES,
                                                                        RENAL_RISK_EXPOSURE, RISKS,
                                                                        RISK_EXPOSURE_LEVELS, RISK_LEVEL_MAP,
                                                                        RISK_OS_HR, RISK_PFS_HR, UNDIAGNOSED)
from vivarium_csu_sanofi_multiple_myeloma.utilities import LogNormalHazardRate

if TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event
    from vivarium.framework.population import SimulantData


class MultipleMyelomaRiskEffects:
    def __init__(self):
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: 'Builder'):
        self.randomness = builder.randomness.get_stream(self.name)
        required_columns = [
            "age",
            "sex",
            models.MULTIPLE_MYELOMA_MODEL_NAME,
            f'{models.MULTIPLE_MYELOMA_1_STATE_NAME}_event_time'
        ]
        self.required_columns = required_columns
        created_columns = list(RISKS)

        # for age and sex, load and sample PFS ratio and OS ratio, get single scalar per risk
        # wrap scalars in a lookup table (dictionary of lookups e.g. pfs_hazard_ratios)
        draw = builder.configuration.input_data.input_draw_number
        progression_hazard_ratio, mortality_hazard_ratio = make_hazard_ratios(draw)
        self.progression_hazard_ratio = builder.lookup.build_table(
            progression_hazard_ratio,
            key_columns=created_columns,
            parameter_columns=['year']
        )
        self.mortality_hazard_ratio = builder.lookup.build_table(
            mortality_hazard_ratio,
            key_columns=created_columns,
            parameter_columns=['year']
        )

        self.population_view = builder.population.get_view(required_columns + created_columns)

        builder.population.initializes_simulants(self.on_initialize_simulants, creates_columns=created_columns,
                                                 requires_columns=required_columns, requires_streams=[self.name])

        states = list(models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES)
        for state, next_state in zip(states, states[1:] + [None]):
            builder.value.register_value_modifier(
                f'{state}.excess_mortality_rate',
                modifier=self.scale_mortality_hazard,
                requires_columns=required_columns,
            )

            if next_state:
                builder.value.register_value_modifier(
                    f'{state}_to_{next_state}.pfs_rate',
                    modifier=self.scale_progression_hazard,
                    requires_columns=required_columns,
                )
        builder.event.register_listener('time_step__cleanup', self.on_time_step_cleanup)

    def scale_mortality_hazard(self, index: pd.Index, mortality_hazard: pd.Series):
        return mortality_hazard * self.mortality_hazard_ratio(index)

    def scale_progression_hazard(self, index: pd.Index, progression_hazard: pd.Series):
        return progression_hazard * self.progression_hazard_ratio(index)

    def on_initialize_simulants(self, pop_data: 'SimulantData'):
        pop = self.population_view.subview(self.required_columns).get(pop_data.index)
        pop_update = self.set_values_on_diagnosis(pop)
        self.population_view.update(pop_update)

    def on_time_step_cleanup(self, event: 'Event'):
        pop = self.population_view.get(event.index)
        newly_with_condition = pop[f'{models.MULTIPLE_MYELOMA_1_STATE_NAME}_event_time'] == event.time
        pop = pop.loc[newly_with_condition]
        if not pop.empty:
            pop_update = self.set_values_on_diagnosis(pop)
            self.population_view.update(pop_update)

    def set_values_on_diagnosis(self, pop):
        """Assuming mask of simulants newly diagnosed, set values for them for risks at diagnosis."""
        risk_exposure = pd.DataFrame(UNDIAGNOSED, columns=list(RISKS), index=pop.index)
        with_condition = pop[models.MULTIPLE_MYELOMA_MODEL_NAME] != models.SUSCEPTIBLE_STATE_NAME
        over_65 = pop.age > 65
        risk_exposure.loc[with_condition & over_65, RISKS.age_at_diagnosis] = RISK_EXPOSURE_LEVELS.over_65
        risk_exposure.loc[with_condition & ~over_65, RISKS.age_at_diagnosis] = RISK_EXPOSURE_LEVELS.under_65

        risk_exposure.loc[with_condition, RISKS.sex_at_diagnosis] = pop.loc[with_condition, 'sex']

        risk_exposure.loc[with_condition, RISKS.race_and_cytogenetic_risk_at_diagnosis] = self.randomness.choice(
            risk_exposure.loc[with_condition].index, choices=list(RACE_AND_CYTO_EXPOSURES.keys()), p=list(RACE_AND_CYTO_EXPOSURES.values()))

        risk_exposure.loc[with_condition, RISKS.renal_function_at_diagnosis] = self.randomness.choice(
            risk_exposure.loc[with_condition].index, choices=list(RENAL_RISK_EXPOSURE.keys()), p=list(RENAL_RISK_EXPOSURE.values()))
        return risk_exposure


def make_hazard_ratios(draw: int):
    idx = pd.MultiIndex.from_product(RISK_LEVEL_MAP.values(), names=RISK_LEVEL_MAP.keys())
    pfs_hazard_ratio = pd.Series(1.0, name='hazard_ratio', index=idx)
    os_hazard_ratio = pd.Series(1.0, name='hazard_ratio', index=idx)

    pfs_hazard_ratio.loc[(UNDIAGNOSED, UNDIAGNOSED, UNDIAGNOSED, UNDIAGNOSED)] = 1.0
    os_hazard_ratio.loc[(UNDIAGNOSED, UNDIAGNOSED, UNDIAGNOSED, UNDIAGNOSED)] = 1.0

    for risk_level_key in itertools.product(*RISK_LEVEL_MAP.values()):
        pfs, os = 1.0, 1.0
        for risk_level in risk_level_key:
            risk_level_pfs, risk_level_os = sample_pfs_and_os(risk_level, draw)
            pfs *= risk_level_pfs
            os *= risk_level_os
        pfs_hazard_ratio.loc[risk_level_key] = pfs
        os_hazard_ratio.loc[risk_level_key] = os

    pfs_hazard_ratio = pfs_hazard_ratio.reset_index()
    os_hazard_ratio = os_hazard_ratio.reset_index()

    # FIXME: Super-duper hack to make lookup table work.  Need at least one continuous parameter.
    pfs_hazard_ratio['year_start'] = 1990
    pfs_hazard_ratio['year_end'] = 2100
    os_hazard_ratio['year_start'] = 1990
    os_hazard_ratio['year_end'] = 2100

    return pfs_hazard_ratio, os_hazard_ratio


def sample_pfs_and_os(risk_level: str, draw: int):
    random_seed = f'{risk_level}_{draw}'
    rs = np.random.RandomState(get_hash(random_seed))
    survival_percentile = rs.random()
    pfs_hr = RISK_PFS_HR[risk_level]
    if isinstance(pfs_hr, tuple):
        pfs_hr = LogNormalHazardRate(*pfs_hr).get_random_variable(survival_percentile)
    os_hr = RISK_OS_HR[risk_level]
    if isinstance(os_hr, tuple):
        os_hr = LogNormalHazardRate(*os_hr).get_random_variable(survival_percentile)
    return pfs_hr, os_hr
