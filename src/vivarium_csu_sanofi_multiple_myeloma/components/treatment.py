from typing import TYPE_CHECKING

import pandas as pd

from vivarium_csu_sanofi_multiple_myeloma.constants import models
from vivarium_csu_sanofi_multiple_myeloma.constants.metadata import SCENARIOS
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import (PROGRESSION_HAZARD_RATE,
                                                                        MORTALITY_HAZARD_RATE,
                                                                        PROBABILITY_RETREAT)

if TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event
    from vivarium.framework.population import SimulantData

TREATMENT_LINES = pd.Index(
    list(models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES),
    name=models.MULTIPLE_MYELOMA_MODEL_NAME
)


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
        models.TREATMENTS.isatuxamib: coverage_data[0],
        models.TREATMENTS.daratumamab: coverage_data[1],
    }, index=TREATMENT_LINES)
    coverage[models.TREATMENTS.residual] = 1 - coverage.sum(axis=1)
    return coverage


def make_progression_hazard_ratio(draw: int):
    index_cols = [models.MULTIPLE_MYELOMA_MODEL_NAME, 'multiple_myeloma_treatment', 'retreated']
    hazard_ratio = pd.DataFrame(columns=index_cols + ['hazard_ratio']).set_index(index_cols)

    hazard_ratio.loc[(models.SUSCEPTIBLE_STATE_NAME, models.TREATMENTS.not_treated, False)] = 1.0

    line = models.MULTIPLE_MYELOMA_1_STATE_NAME
    hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, False)] = PROGRESSION_HAZARD_RATE.FIRST_LINE_ISA.get_random_variable(draw)
    hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, False)] = PROGRESSION_HAZARD_RATE.FIRST_LINE_DARA.get_random_variable(draw)
    hazard_ratio.loc[(line, models.TREATMENTS.residual, False)] = PROGRESSION_HAZARD_RATE.FIRST_LINE_RESIDUAL.get_random_variable(draw)

    for line in TREATMENT_LINES.tolist()[1:]:
        hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, True)] = PROGRESSION_HAZARD_RATE.LATER_LINE_ISA.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, False)] = PROGRESSION_HAZARD_RATE.LATER_LINE_ISA_RETREAT.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, True)] = PROGRESSION_HAZARD_RATE.LATER_LINE_DARA.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, False)] = PROGRESSION_HAZARD_RATE.LATER_LINE_DARA_RETREAT.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.residual, True)] = PROGRESSION_HAZARD_RATE.LATER_LINE_RESIDUAL.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.residual, False)] = PROGRESSION_HAZARD_RATE.LATER_LINE_RESIDUAL.get_random_variable(draw)

    hazard_ratio = hazard_ratio.reset_index()
    # FIXME: Super-duper hack to make lookup table work.  Need at least one continuous parameter.
    hazard_ratio['year_start'] = 1990
    hazard_ratio['year_end'] = 2100
    return hazard_ratio


def make_mortality_hazard_ratio(draw: int):
    index_cols = [models.MULTIPLE_MYELOMA_MODEL_NAME, 'multiple_myeloma_treatment', 'retreated']
    hazard_ratio = pd.DataFrame(columns=index_cols + ['hazard_ratio']).set_index(index_cols)

    hazard_ratio.loc[(models.SUSCEPTIBLE_STATE_NAME, models.TREATMENTS.not_treated, False)] = 1.0

    line = models.MULTIPLE_MYELOMA_1_STATE_NAME
    hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, False)] = MORTALITY_HAZARD_RATE.FIRST_LINE_ISA.get_random_variable(draw)
    hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, False)] = MORTALITY_HAZARD_RATE.FIRST_LINE_DARA.get_random_variable(draw)
    hazard_ratio.loc[(line, models.TREATMENTS.residual, False)] = MORTALITY_HAZARD_RATE.FIRST_LINE_RESIDUAL.get_random_variable(draw)

    for line in TREATMENT_LINES.tolist()[1:]:
        hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, False)] = MORTALITY_HAZARD_RATE.LATER_LINE_ISA.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.isatuxamib, True)] = MORTALITY_HAZARD_RATE.LATER_LINE_ISA_RETREAT.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, False)] = MORTALITY_HAZARD_RATE.LATER_LINE_DARA.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.daratumamab, True)] = MORTALITY_HAZARD_RATE.LATER_LINE_DARA_RETREAT.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.residual, True)] = MORTALITY_HAZARD_RATE.LATER_LINE_RESIDUAL.get_random_variable(draw)
        hazard_ratio.loc[(line, models.TREATMENTS.residual, False)] = MORTALITY_HAZARD_RATE.LATER_LINE_RESIDUAL.get_random_variable(draw)

    hazard_ratio = hazard_ratio.reset_index()
    # FIXME: Super-duper hack to make lookup table work.  Need at least one continuous parameter.
    hazard_ratio['year_start'] = 1990
    hazard_ratio['year_end'] = 2100
    return hazard_ratio


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
        self.coverage_2021 = make_treatment_coverage(2021, scenario)
        self.coverage_2025 = make_treatment_coverage(2025, scenario)

        # What treatment are they currently on.
        self.treatment_column = 'multiple_myeloma_treatment'
        # Did they previously recieve isatuxamib or daratumumab
        self.retreatment_eligible_column = 'retreatment_eligible'
        self.retreated_column = 'retreated'
        columns_created = [
            self.treatment_column,
            self.retreatment_eligible_column,
            self.retreated_column
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
            self.retreatment_eligible_column: 'unknown',
            self.retreated_column: False
        }, index=pop_data.index)
        with_mm = initial_mm_state.loc[
            initial_mm_state[models.MULTIPLE_MYELOMA_MODEL_NAME] != models.SUSCEPTIBLE_STATE_NAME,
            models.MULTIPLE_MYELOMA_MODEL_NAME
        ]
        current_coverage = pd.DataFrame(
            current_coverage.loc[with_mm.values].values,
            columns=current_coverage.columns,
            index=with_mm.index,
        )
        treatment = self.randomness.choice(
            with_mm.index,
            choices=current_coverage.columns.tolist(),
            p=current_coverage.values,
            additional_key='initial_treatment_status',
        )
        pop_update.loc[with_mm.index, self.treatment_column] = treatment
        dara_or_isa = pop_update[self.treatment_column].isin(
            [models.TREATMENTS.daratumamab, models.TREATMENTS.isatuxamib]
        )
        pop_update.loc[dara_or_isa, self.retreatment_eligible_column] = 'true'
        self.population_view.update(pop_update)

    def on_time_step_cleanup(self, event: 'Event'):
        pop = self.population_view.get(event.index)

        retreat_mask = self.randomness.get_draw(pop.index, additional_key='retreat') < PROBABILITY_RETREAT
        retreatment_eligible = pop[self.retreatment_eligible_column] == 'true'
        retreatment_unknown = pop[self.retreatment_eligible_column] == 'unknown'
        retreatment_ineligible = pop[self.retreatment_eligible_column] == 'false'

        coverage = self.get_current_coverage(event.time)
        lines = TREATMENT_LINES.tolist()
        for current_line, previous_line in zip(lines, [None] + lines[:-1]):
            # First, unpack probabilities for the current and previous line.
            p_isa_c = coverage.at[current_line, models.TREATMENTS.isatuxamib]
            p_dara_c = coverage.at[current_line, models.TREATMENTS.daratumamab]
            p_resid_c = coverage.at[current_line, models.TREATMENTS.residual]
            if previous_line:
                p_isa_p = coverage.at[previous_line, models.TREATMENTS.isatuxamib]
                p_dara_p = coverage.at[previous_line, models.TREATMENTS.daratumamab]
                p_resid_p = coverage.at[previous_line, models.TREATMENTS.residual]
            else:
                p_isa_p, p_dara_p, p_resid_p = 0., 0., 1.

            # Our base filter, which we'll partition.
            new_treatment_line = pop[f'{current_line}_event_time'] == event.time

            # First group, getting their 2nd+ round of isa/dara
            retreat = new_treatment_line & retreatment_eligible & retreat_mask
            retreat_choices = [models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab]
            retreat_probs = [p_isa_c / (p_isa_c + p_dara_c), p_dara_c / (p_isa_c + p_dara_c)]

            pop.loc[retreat, self.treatment_column] = self.randomness.choice(
                pop.loc[retreat].index,
                choices=retreat_choices,
                p=retreat_probs,
            )
            pop.loc[retreat, self.retreatment_eligible_column] = 'true'  # This is a no-op.  Here for clarity.
            pop.loc[retreat, self.retreated_column] = True

            # Second group, got 1 dose of isa/dara, but can't receive again.
            no_retreat = (
                (new_treatment_line & retreatment_eligible & ~retreat_mask)
                | (new_treatment_line & retreatment_ineligible)
            )

            pop.loc[no_retreat, self.treatment_column] = models.TREATMENTS.residual
            pop.loc[no_retreat, self.retreatment_eligible_column] = 'false'
            pop.loc[no_retreat, self.retreated_column] = False  # This is a no-op.  Here for clarity.

            # Third group, getting their first dose of isa/dara and determining if they can be retreated.
            unknown_retreat = new_treatment_line & retreatment_unknown
            unknown_choices = [models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab, models.TREATMENTS.residual]
            # TODO: add a link to the docs when they're live for this algorithm.
            old_to_new_scale = (p_isa_p + p_dara_p) / (p_isa_c + p_dara_c)
            final_scale = (1 - PROBABILITY_RETREAT * old_to_new_scale) / p_resid_p
            p_isa = p_isa_c * final_scale
            p_dara = p_dara_c * final_scale
            if p_isa + p_dara > 1:
                p_isa, p_dara = p_isa / (p_isa + p_dara), p_dara / (p_isa + p_dara)
            p_resid = 1 - p_isa - p_dara
            unknown_treatment_probs = [p_isa, p_dara, p_resid]

            pop.loc[unknown_retreat, self.treatment_column] = self.randomness.choice(
                pop.loc[unknown_retreat].index,
                choices=unknown_choices,
                p=unknown_treatment_probs,
            )
            isa_or_dara = pop[self.treatment_column].isin([
                models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab
            ])
            pop.loc[unknown_retreat & isa_or_dara, self.retreatment_eligible_column] = 'true'
            # These are no-ops.  Here for clarity.
            pop.loc[unknown_retreat & ~isa_or_dara, self.retreatment_eligible_column] = 'unknown'
            pop.loc[unknown_retreat, self.retreated_column] = False

        self.population_view.update(pop)

    def get_current_coverage(self, time: pd.Timestamp) -> pd.DataFrame:
        """Get a df with columns: [TREATMENTS.isatuxamib, TREATMENTS.daratumamab, TREATMENTS.residual]
        indexed by multiple myeloma state."""
        if time.year < 2021:
            return self.coverage_2021
        elif time.year > 2025:
            return self.coverage_2025
        t = (time - pd.Timestamp('2021-01-01')) / (pd.Timestamp('2026-01-01') - pd.Timestamp('2021-01-01'))

        treatments = [models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab]
        slope = self.coverage_2025[treatments] - self.coverage_2021[treatments]
        coverage = self.coverage_2021[treatments] + slope * t
        coverage[models.TREATMENTS.residual] = 1 - coverage.sum(axis=1)
        return coverage


class MultipleMyelomaTreatmentEffect:

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def setup(self, builder: 'Builder') -> None:
        draw = builder.configuration.input_data.input_draw_number
        required_columns = [models.MULTIPLE_MYELOMA_MODEL_NAME, 'multiple_myeloma_treatment', 'retreated']
        progression_hazard_ratio = make_progression_hazard_ratio(draw)
        self.progression_hazard_ratio = builder.lookup.build_table(
            progression_hazard_ratio,
            key_columns=required_columns,
            parameter_columns=['year']
        )
        mortality_hazard_ratio = make_mortality_hazard_ratio(draw)
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
                    f'{state}_to_{next_state}.transition_rate',
                    modifier=self.scale_progression_hazard,
                    requires_columns=required_columns,
                )

    def scale_mortality_hazard(self, index: pd.Index, mortality_hazard: pd.Series):
        return mortality_hazard * self.mortality_hazard_ratio(index)

    def scale_progression_hazard(self, index: pd.Index, progression_hazard: pd.Series):
        return progression_hazard * self.progression_hazard_ratio(index)

