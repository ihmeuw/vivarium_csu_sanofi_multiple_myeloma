from typing import TYPE_CHECKING

import pandas as pd

from vivarium_csu_sanofi_multiple_myeloma.constants import models
from vivarium_csu_sanofi_multiple_myeloma.constants.metadata import SCENARIOS

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

PROBABILITY_RETREAT = 0.15


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
        self.coverage_2025 = make_treatment_coverage(2021, scenario)

        # What treatment are they currently on.
        self.treatment_column = 'multiple_myeloma_treatment'
        # Did they previously recieve isatuxamib or daratumumab
        self.previous_isa_or_dara_column = 'previous_isa_or_dara'
        columns_created = [
            self.treatment_column,
            self.previous_isa_or_dara_column,
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
            self.treatment_column: models.TREATMENTS.no_treatment,
            self.previous_isa_or_dara_column: False,
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
        pop_update.loc[dara_or_isa, self.previous_isa_or_dara_column] = True
        self.population_view.update(pop_update)

    def on_time_step_cleanup(self, event: 'Event'):
        pop = self.population_view.get(event.index)
        retreat = self.randomness.get_draw(pop.index, additional_key='retreat') < PROBABILITY_RETREAT
        previous_isa_or_dara = pop[self.previous_isa_or_dara_column]
        coverage = self.get_current_coverage(event.time)
        for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES:
            line_coverage = coverage.loc[s]
            new_treatment_line = pop[f'{s}_event_time'] == event.time
            pop.loc[new_treatment_line & previous_isa_or_dara & retreat, self.treatment_column] = (
                self.randomness.choice(
                    pop.loc[new_treatment_line & previous_isa_or_dara & retreat].index,
                    choices=[models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab],
                    p=line_coverage.loc[[models.TREATMENTS.isatuxamib, models.TREATMENTS.daratumamab]],
                )
            )
            pop.loc[new_treatment_line & previous_isa_or_dara & ~retreat, self.treatment_column] = (
                models.TREATMENTS.residual
            )
            pop.loc[new_treatment_line & ~previous_isa_or_dara, self.treatment_column] = (
                self.randomness.choice(
                    pop.loc[new_treatment_line & ~previous_isa_or_dara].index,
                    choices=coverage.columns.tolist(),
                    p=line_coverage,
                )
            )
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
