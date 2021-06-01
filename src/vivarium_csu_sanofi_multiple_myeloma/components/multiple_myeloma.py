import pandas as pd
from vivarium_public_health.disease import (DiseaseState, DiseaseModel, SusceptibleState,
                                            RateTransition as RateTransition_, RecoveredState)

from vivarium_csu_sanofi_multiple_myeloma.constants import models, data_keys, data_values
import vivarium_csu_sanofi_multiple_myeloma.paths as paths

import typing

if typing.TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event


class DiseaseStateHazard(DiseaseState):
    def __init__(self):
        self.time_since_entrance_col = f'{self.state_id}_time_since_entrance'

    @property
    def columns_created(self):
        return super().columns_created + [self.time_since_entrance_col]

    def setup(self, builder: 'Builder'):
        super().setup(builder)
        # hazard_rate_data = pd.DataFrame(
        #     columns=[f'{self.time_since_entrance_col}_start',
        #              f'{self.time_since_entrance_col}_end', 'rate'])
        # load in YX dataset - time start, time end, rate ( or synthetic)
        hazard_rate_data = pd.read_csv(paths.MORTALITY_FIRST_LINE_PATH)

        # noinspection PyAttributeOutsideInit
        self.hazard_rate = builder.lookup.build_table(
            hazard_rate_data, parameter_columns=[self.time_since_entrance_col])

    def on_initialize_simulants(self, pop_data):
        super().on_initialize_simulants(pop_data)
        entrance_time = self.population_view.subview([self.event_time_column]).get(pop_data.index)
        time_since_entrance = pd.Series(
            data=self.clock()-entrance_time, index=pop_data.index, name=self.time_since_entrance_col)
        self.population_view.update(time_since_entrance)

    # TODO: add on time step to update the time since entrance (__prepare?)
    def on_time_step(self, event: 'Event'):
        super().on_time_step(event)
        entrance_time = self.population_view.subview([self.event_time_column]).get(event.index)
        time_since_entrance = pd.Series(
            data=self.clock()-entrance_time, index=event.index, name=self.time_since_entrance_col)
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
            'incidence_rate': lambda _, builder: builder.data.load(data_keys.CERVICAL_CANCER.HRHPV_INCIDENCE_RATE)
        }
    )

    return DiseaseModel('multiple_myeloma', states=[susceptible, mm])
