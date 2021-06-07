from datetime import timedelta
import pandas as pd
#from typing import Union
import typing

from vivarium_public_health.disease import (DiseaseState, DiseaseModel, SusceptibleState,
                                            RateTransition as RateTransition_, RecoveredState)

from vivarium_csu_sanofi_multiple_myeloma.constants import models, data_keys, data_values
import vivarium_csu_sanofi_multiple_myeloma.paths as paths


if typing.TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event


class DiseaseStateHazard(DiseaseState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.time_since_entrance_col = f'{self.state_id}_time_since_entrance'

    @property
    def columns_created(self):
        return super().columns_created + [self.time_since_entrance_col]

    def setup(self, builder: 'Builder'):
        super().setup(builder)
        self.step_size = builder.time.step_size()
        draw = builder.configuration.input_data.input_draw_number
        builder.event.register_listener('time_step', self.on_time_step)

        # Load in time-variant hazard rate
        hazard_rate_data = pd.read_csv(paths.MORTALITY_FIRST_LINE_PATH)
        hazard_rate_data[f'{self.time_since_entrance_col}_start'] = hazard_rate_data[
            f'{self.time_since_entrance_col}_start'].astype(int).multiply(self.step_size().days)
        hazard_rate_data[f'{self.time_since_entrance_col}_end'] = hazard_rate_data[
            f'{self.time_since_entrance_col}_end'].astype(int).multiply(self.step_size().days)

        # FIXME: Hack for draw-level hazard rate
        hazard_rate_data = hazard_rate_data[
            [f'{self.time_since_entrance_col}_start', f'{self.time_since_entrance_col}_end', f'draw_{draw}']]

        # noinspection PyAttributeOutsideInit
        self.hazard_rate = builder.lookup.build_table(
            hazard_rate_data, parameter_columns=[self.time_since_entrance_col])

    def on_initialize_simulants(self, pop_data):
        super().on_initialize_simulants(pop_data)
        entrance_time = self.population_view.subview([self.event_time_column]).get(pop_data.index)
        time_since_entrance = pd.Series(
            data=(self.clock()-entrance_time.iloc[:, 0]).dt.days,
            index=pop_data.index,
            name=self.time_since_entrance_col)
        self.population_view.update(time_since_entrance)

    def on_time_step(self, event: 'Event'):
        entrance_time = self.population_view.subview([self.event_time_column]).get(event.index)
        time_since_entrance = pd.Series(
            data=(event.time-entrance_time.iloc[:, 0]).dt.days,
            index=event.index,
            name=self.time_since_entrance_col)
        self.population_view.update(time_since_entrance)

    def compute_excess_mortality_rate(self, index):
        excess_mortality_rate = pd.Series(0, index=index)
        with_condition = self.with_condition(index)
        excess_mortality_rate.loc[with_condition] = self.hazard_rate(with_condition)
        return excess_mortality_rate


def MultipleMyeloma():
    susceptible = SusceptibleState(models.MULTIPLE_MYELOMA_MODEL_NAME)
    mm = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_STATE_NAME
    )

    # Add transitions for Susceptible state
    susceptible.allow_self_transitions()
    susceptible.add_transition(
        mm,
        source_data_type='rate',
        get_data_functions={
            'incidence_rate': lambda _, builder: builder.data.load(data_keys.MULTIPLE_MYELOMA.INCIDENCE_RATE)
        }
    )

    return DiseaseModel('multiple_myeloma', states=[susceptible, mm])
