from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from vivarium.framework.randomness import get_hash

from vivarium_csu_sanofi_multiple_myeloma.constants import models
from vivarium_csu_sanofi_multiple_myeloma.constants.metadata import SCENARIOS
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import (OS_HR, PFS_HR, PROBABILITY_RETREAT,
                                                                        REGISTRY_ENROLL_PROBABILITY)
from vivarium_csu_sanofi_multiple_myeloma.utilities import LogNormalHazardRate

if TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event
    from vivarium.framework.population import SimulantData

TREATMENT_LINES = pd.Index(
    list(models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES),
    name=models.MULTIPLE_MYELOMA_MODEL_NAME
)


def make_treatment_coverage(year, scenario):
    scalar_2019 = (2019-2016)/(2021-2016)
    scalar_2020 = (2020-2016)/(2021-2016)
    coverages = {
        (2016, SCENARIOS.baseline): (
            [0.000, 0.000, 0.000, 0.000, 0.000],
            [0.000, 0.000, 0.000, 0.000, 0.000],
        ),
        (2019, SCENARIOS.baseline): (
            [0.000, 0.000, 0.000, 0.000, 0.000],
            [0.000, 0.198 * scalar_2019, 0.323 * scalar_2019, 0.365 * scalar_2019, 0.3011 * scalar_2019],
        ),
        (2020, SCENARIOS.baseline): (
            [0.000, 0.000, 0.000, 0.000, 0.000],
            [0.000, 0.198 * scalar_2020, 0.323 * scalar_2020, 0.365 * scalar_2020, 0.3011 * scalar_2020]
        ),
        (2021, SCENARIOS.baseline): (
            [0.000, 0.005, 0.010, 0.033, 0.033],
            [0.029, 0.198, 0.323, 0.365, 0.3011],
        ),
        (2025, SCENARIOS.baseline): (
            [0.000, 0.100, 0.090, 0.070, 0.070],
            [0.34, 0.34, 0.34, 0.34, 0.34],
        ),
        (2025, SCENARIOS.alternative): (
            [0.100, 0.100, 0.090, 0.070, 0.070],
            [0.34, 0.34, 0.34, 0.34, 0.34],
        )
    }
    for target_year in (2016, 2019, 2020, 2021):
        coverages[(target_year, SCENARIOS.alternative)] = coverages[(target_year, SCENARIOS.baseline)]

    coverage_data = coverages[(year, scenario)]
    coverage = pd.DataFrame({
        models.TREATMENTS.isatuximab: coverage_data[0],
        models.TREATMENTS.daratumumab: coverage_data[1],
    }, index=TREATMENT_LINES)
    coverage[models.TREATMENTS.residual] = 1 - coverage.sum(axis=1)
    return coverage


def make_hazard_ratios(draw: int):
    index_cols = [models.MULTIPLE_MYELOMA_MODEL_NAME, 'multiple_myeloma_treatment', 'retreated']
    pfs_hazard_ratio = pd.DataFrame(columns=index_cols + ['hazard_ratio']).set_index(index_cols)
    os_hazard_ratio = pd.DataFrame(columns=index_cols + ['hazard_ratio']).set_index(index_cols)

    pfs_hazard_ratio.loc[(models.SUSCEPTIBLE_STATE_NAME, models.TREATMENTS.not_treated, False)] = 1.0
    os_hazard_ratio.loc[(models.SUSCEPTIBLE_STATE_NAME, models.TREATMENTS.not_treated, False)] = 1.0

    for key in PFS_HR:
        random_seed = '_'.join([str(k) for k in key] + [str(draw)])
        rs = np.random.RandomState(get_hash(random_seed))
        survival_percentile = rs.random()
        pfs_hazard_ratio.loc[key] = LogNormalHazardRate(*PFS_HR[key]).get_random_variable(survival_percentile)
        os_hazard_ratio.loc[key] = LogNormalHazardRate(*OS_HR[key]).get_random_variable(survival_percentile)

    for key in set(OS_HR).difference(PFS_HR):
        random_seed = '_'.join([str(k) for k in key] + [str(draw)])
        rs = np.random.RandomState(get_hash(random_seed))
        survival_percentile = rs.random()
        os_hazard_ratio.loc[key] = LogNormalHazardRate(*OS_HR[key]).get_random_variable(survival_percentile)

    pfs_hazard_ratio = pfs_hazard_ratio.reset_index()
    os_hazard_ratio = os_hazard_ratio.reset_index()

    # FIXME: Super-duper hack to make lookup table work.  Need at least one continuous parameter.
    pfs_hazard_ratio['year_start'] = 1990
    pfs_hazard_ratio['year_end'] = 2100
    os_hazard_ratio['year_start'] = 1990
    os_hazard_ratio['year_end'] = 2100

    return pfs_hazard_ratio, os_hazard_ratio


class MultipleMyelomaTreatmentCoverage:

    configuration_defaults = {
        'mm_treatment_scenario': SCENARIOS.baseline,
    }

    @property
    def name(self):
        return self.__class__.__name__

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: 'Builder') -> None:
        self.clock = builder.time.clock()
        self.randomness = builder.randomness.get_stream(self.name)

        scenario = builder.configuration.mm_treatment_scenario
        assert scenario in SCENARIOS
        self.coverage = {}
        for year in (2016, 2019, 2020, 2021, 2025):
            self.coverage[year] = make_treatment_coverage(year, scenario)

        # What treatment are they currently on.
        self.treatment_column = 'multiple_myeloma_treatment'
        # Did they previously receive isatuximab or daratumumab
        self.ever_isa_or_dara_column = 'ever_isa_or_dara'
        self.retreated_column = 'retreated'
        self.registry_evaluation_status = 'registry_evaluation_status'  # 4 potential values: unevaluated, eligible, enrolled
        self.registry_evaluation_date = 'registry_evaluation_date'
        self.ever_isa = 'ever_isa'
        self.registry_start_date = pd.Timestamp('2021-01-01')
        columns_created = [
            self.treatment_column,
            self.ever_isa_or_dara_column,
            self.retreated_column,
            self.registry_evaluation_status,
            self.registry_evaluation_date,
            self.ever_isa
        ]
        columns_required = (
            [models.MULTIPLE_MYELOMA_MODEL_NAME]
            + [f'{s}_event_time' for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES]
        )
        self.population_view = builder.population.get_view(columns_required + columns_created)
        builder.population.initializes_simulants(
            self.on_initialize_simulants,
            creates_columns=columns_created,
            requires_columns=columns_required,
            requires_streams=[self.randomness.key],
        )

        builder.event.register_listener('time_step__cleanup', self.on_time_step_cleanup)

    def on_initialize_simulants(self, pop_data: 'SimulantData') -> None:
        current_coverage = self.get_current_coverage(pop_data.creation_time)
        initial_mm_state = self.population_view.subview([models.MULTIPLE_MYELOMA_MODEL_NAME]).get(pop_data.index)
        pop_update = pd.DataFrame({
            self.treatment_column: models.TREATMENTS.not_treated,
            self.ever_isa_or_dara_column: False,
            self.retreated_column: False,
            self.registry_evaluation_status: 'unevaluated',
            self.registry_evaluation_date: pd.NaT,
            self.ever_isa: False
        }, index=pop_data.index)
        with_mm = initial_mm_state.loc[
            initial_mm_state[models.MULTIPLE_MYELOMA_MODEL_NAME] != models.SUSCEPTIBLE_STATE_NAME,
            models.MULTIPLE_MYELOMA_MODEL_NAME
        ]
        pop_update.loc[with_mm.index, self.treatment_column] = models.TREATMENTS.residual
        self.population_view.update(pop_update)

    def on_time_step_cleanup(self, event: 'Event'):
        pop = self.population_view.get(event.index)

        retreat_mask = self.randomness.get_draw(pop.index, additional_key='retreat') < PROBABILITY_RETREAT
        ever_isa_or_dara = pop[self.ever_isa_or_dara_column].copy()
        ever_isa = pop[self.ever_isa].copy()
        registry_eligible = pd.Series(False, index=pop.index)
        registry_mask = self.randomness.get_draw(pop.index, additional_key='registry') < REGISTRY_ENROLL_PROBABILITY

        proportion_ever_isa_or_dara = 0
        coverage = self.get_current_coverage(event.time)
        lines = TREATMENT_LINES.tolist()

        for current_line, previous_line in zip(lines, [None] + lines[:-1]):
            # First, unpack probabilities for the current and previous line.
            p_isa = coverage.at[current_line, models.TREATMENTS.isatuximab]
            p_dara = coverage.at[current_line, models.TREATMENTS.daratumumab]
            p_resid = coverage.at[current_line, models.TREATMENTS.residual]

            # Our base filter, which we'll partition.
            new_treatment_line = pop[f'{current_line}_event_time'] == event.time

            # First group, never had isa/dara
            naive = new_treatment_line & ~ever_isa_or_dara
            naive_choices = [models.TREATMENTS.isatuximab, models.TREATMENTS.daratumumab, models.TREATMENTS.residual]
            rescale_probabilities = lambda p1, p2, e: (p1 - e * PROBABILITY_RETREAT * p1/(p1 + p2))/ (1 - e)
            p_isa_naive = rescale_probabilities(p_isa, p_dara, proportion_ever_isa_or_dara)
            p_dara_naive = rescale_probabilities(p_dara, p_isa, proportion_ever_isa_or_dara)
            if p_isa_naive + p_dara_naive > 1:
                p_isa_naive = p_isa_naive / (p_isa_naive + p_dara_naive)
                p_dara_naive = p_dara_naive / (p_isa_naive + p_dara_naive)
            p_resid_naive = 1 - p_isa_naive - p_dara_naive
            naive_probs = [p_isa_naive, p_dara_naive, p_resid_naive]
            pop.loc[naive, self.treatment_column] = self.randomness.choice(
                pop.loc[naive].index,
                choices=naive_choices,
                p=naive_probs,
            )
            isa_or_dara = pop[self.treatment_column].isin([
                models.TREATMENTS.isatuximab, models.TREATMENTS.daratumumab
            ])
            isa = pop[self.treatment_column] == models.TREATMENTS.isatuximab
            pop.loc[naive & isa, self.ever_isa] = True
            pop.loc[naive & isa_or_dara, self.ever_isa_or_dara_column] = True

            # These are no-ops.  Here for clarity.
            pop.loc[naive & ~isa_or_dara, self.ever_isa_or_dara_column] = False
            pop.loc[naive, self.retreated_column] = False
            # ever_x = (1 - PROBABILITY_RETREAT) * ever_x-1 + coverage_x
            proportion_ever_isa_or_dara = (1 - PROBABILITY_RETREAT) * proportion_ever_isa_or_dara + p_isa + p_dara

            # Second group, simulants w/prior exposure to isa/dara, and will be retreated this line
            retreat = new_treatment_line & ever_isa_or_dara & retreat_mask
            retreat_choices = [models.TREATMENTS.isatuximab, models.TREATMENTS.daratumumab]
            retreat_probs = [p_isa / (p_isa + p_dara), p_dara / (p_isa + p_dara)]

            pop.loc[retreat, self.treatment_column] = self.randomness.choice(
                pop.loc[retreat].index,
                choices=retreat_choices,
                p=retreat_probs,
            )
            pop.loc[retreat, self.ever_isa_or_dara_column] = True  # This is a no-op.  Here for clarity.
            pop.loc[retreat, self.retreated_column] = True

            # Third group, got 1 dose of isa/dara, but won't receive one this line, may receive again
            no_retreat = new_treatment_line & ever_isa_or_dara & ~retreat_mask

            pop.loc[no_retreat, self.treatment_column] = models.TREATMENTS.residual
            # pop.loc[no_retreat, ever_isa_or_dara] does not change
            # pop.loc[no_retreat, retreated] does not change

            # Build registry mask
            if self.registry_start_date <= event.time:
                registry_eligible = registry_eligible | (~ever_isa & isa)

        if self.registry_start_date <= event.time:
            pop.loc[registry_eligible & registry_mask, self.registry_evaluation_status] = 'enrolled'
            pop.loc[registry_eligible, self.registry_evaluation_date] = event.time
            pop.loc[registry_eligible & ~registry_mask, self.registry_evaluation_status] = 'eligible'
        self.population_view.update(pop)

    def get_current_coverage(self, time: pd.Timestamp) -> pd.DataFrame:
        """Get a df with columns: [TREATMENTS.isatuximab, TREATMENTS.daratumumab, TREATMENTS.residual]
        indexed by multiple myeloma state."""
        if time.year < 2016:
            return self.coverage[2016]
        elif time.year < 2019:
            upper_year = 2019
            lower_year = 2016
        elif time.year < 2020:
            upper_year = 2020
            lower_year = 2019
        elif time.year < 2021:
            upper_year = 2021
            lower_year = 2020
        elif time.year < 2025:
            upper_year = 2025
            lower_year = 2021
        elif time.year >= 2025:
            return self.coverage[2025]
        t = (time - pd.Timestamp(f'{lower_year}-01-01')) / (pd.Timestamp(f'{upper_year}-01-01') - pd.Timestamp(f'{lower_year}-01-01'))

        treatments = [models.TREATMENTS.isatuximab, models.TREATMENTS.daratumumab]
        slope = self.coverage[upper_year][treatments] - self.coverage[lower_year][treatments]
        coverage = self.coverage[lower_year][treatments] + slope * t
        coverage[models.TREATMENTS.residual] = 1 - coverage.sum(axis=1)
        return coverage


class MultipleMyelomaTreatmentEffect:

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def setup(self, builder: 'Builder') -> None:
        draw = builder.configuration.input_data.input_draw_number
        required_columns = [models.MULTIPLE_MYELOMA_MODEL_NAME, 'multiple_myeloma_treatment', 'retreated']
        progression_hazard_ratio, mortality_hazard_ratio = make_hazard_ratios(draw)
        self.progression_hazard_ratio = builder.lookup.build_table(
            progression_hazard_ratio,
            key_columns=required_columns,
            parameter_columns=['year']
        )
        self.mortality_hazard_ratio = builder.lookup.build_table(
            mortality_hazard_ratio,
            key_columns=required_columns,
            parameter_columns=['year']
        )
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

    def scale_mortality_hazard(self, index: pd.Index, mortality_hazard: pd.Series):
        return mortality_hazard * self.mortality_hazard_ratio(index)

    def scale_progression_hazard(self, index: pd.Index, progression_hazard: pd.Series):
        return progression_hazard * self.progression_hazard_ratio(index)
