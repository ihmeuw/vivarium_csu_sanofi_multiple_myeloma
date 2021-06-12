from datetime import timedelta
import pandas as pd
from typing import Union
import typing

from vivarium_public_health.disease import (DiseaseState, DiseaseModel, SusceptibleState,
                                            RateTransition as RateTransition_, RecoveredState)

from vivarium_csu_sanofi_multiple_myeloma.constants import models, data_keys, data_values
import vivarium_csu_sanofi_multiple_myeloma.paths as paths


if typing.TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event


class RateTransition(RateTransition_):
    def setup(self, builder):
        super().setup(builder)
        rate_data, _ = self.load_transition_rate_data(builder)
        self.base_rate = builder.lookup.build_table(
            rate_data, parameter_columns=[self.input_state.time_since_entrance_col])

    def load_transition_rate_data(self, builder):
        rate_data = load_hazard_rate(builder, self.input_state.state_id, "incidence")
        pipeline_name = f'{self.input_state.state_id}_to_{self.output_state.state_id}.transition_rate'
        return rate_data, pipeline_name


class DiseaseStateHazard(DiseaseState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.time_since_entrance_col = f'{self.state_id}_time_since_entrance'

    @property
    def columns_created(self):
        return super().columns_created + [self.time_since_entrance_col]

    def setup(self, builder: 'Builder'):
        super().setup(builder)
        builder.event.register_listener('time_step', self.on_time_step)

        hazard_rate_data = load_hazard_rate(builder, self.state_id, "mortality")
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

    def add_transition(self, output, source_data_type=None, get_data_functions=None, **kwargs):
        if source_data_type == 'hazard_rate':
            return RateTransition(self, output)
        else:
            return super().add_transition(output, source_data_type, get_data_functions, **kwargs)


def MultipleMyeloma():
    susceptible = SusceptibleState(models.MULTIPLE_MYELOMA_MODEL_NAME)
    mm_1 = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_1_STATE_NAME
    )

    # Add transitions for Susceptible state
    susceptible.allow_self_transitions()
    susceptible.add_transition(
        mm_1,
        source_data_type='rate',
        get_data_functions={
            'incidence_rate': lambda _, builder: builder.data.load(data_keys.MULTIPLE_MYELOMA.INCIDENCE_RATE)
        }
    )

    mm_2 = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_2_STATE_NAME
    )
    mm_1.allow_self_transitions()
    mm_1.add_transition(
        mm_2,
        source_data_type="hazard_rate"
    )

    mm_3 = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_3_STATE_NAME
    )
    mm_2.allow_self_transitions()
    mm_2.add_transition(
        mm_3,
        source_data_type="hazard_rate"
    )

    mm_4 = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_4_STATE_NAME
    )
    mm_3.allow_self_transitions()
    mm_3.add_transition(
        mm_4,
        source_data_type="hazard_rate"
    )

    mm_5 = DiseaseStateHazard(
        models.MULTIPLE_MYELOMA_5_STATE_NAME
    )
    mm_4.allow_self_transitions()
    mm_4.add_transition(
        mm_5,
        source_data_type="hazard_rate"
    )

    mm_5.allow_self_transitions()

    return DiseaseModel('multiple_myeloma', states=[susceptible, mm_1, mm_2, mm_3, mm_4, mm_5])


def load_hazard_rate(builder: 'Builder', state_id, measure):
    # Load in time-variant hazard rate
    step_size = builder.time.step_size()
    draw = builder.configuration.input_data.input_draw_number
    data_map = {
        (models.MULTIPLE_MYELOMA_1_STATE_NAME, "mortality"): paths.MORTALITY_FIRST_LINE_PATH,
        (models.MULTIPLE_MYELOMA_2_STATE_NAME, "mortality"): paths.MORTALITY_SECOND_LINE_PATH,
        (models.MULTIPLE_MYELOMA_3_STATE_NAME, "mortality"): paths.MORTALITY_THIRD_LINE_PATH,
        (models.MULTIPLE_MYELOMA_4_STATE_NAME, "mortality"): paths.MORTALITY_FOURTH_LINE_PATH,
        (models.MULTIPLE_MYELOMA_5_STATE_NAME, "mortality"): paths.MORTALITY_FIFTH_LINE_PATH,
        (models.MULTIPLE_MYELOMA_1_STATE_NAME, "incidence"): paths.INCIDENCE_FIRST_LINE_PATH,
        (models.MULTIPLE_MYELOMA_2_STATE_NAME, "incidence"): paths.INCIDENCE_SECOND_LINE_PATH,
        (models.MULTIPLE_MYELOMA_3_STATE_NAME, "incidence"): paths.INCIDENCE_THIRD_LINE_PATH,
        (models.MULTIPLE_MYELOMA_4_STATE_NAME, "incidence"): paths.INCIDENCE_FOURTH_LINE_PATH,
        (models.MULTIPLE_MYELOMA_5_STATE_NAME, "incidence"): paths.INCIDENCE_FIFTH_LINE_PATH,
    }
    hazard_rate_data = pd.read_csv(data_map.get((state_id, measure), paths.MORTALITY_FIRST_LINE_PATH))
    hazard_rate_data[f'{state_id}_time_since_entrance_start'] = hazard_rate_data[
        f'time_since_entrance_start'].astype(int).multiply(step_size().days)
    hazard_rate_data[f'{state_id}_time_since_entrance_end'] = hazard_rate_data[
        f'time_since_entrance_end'].astype(int).multiply(step_size().days)
    hazard_rate_data['rate'] = hazard_rate_data[f'draw_{draw}']
    # FIXME: Hack for draw-level hazard rate
    hazard_rate_data = hazard_rate_data[
        [f'{state_id}_time_since_entrance_start', f'{state_id}_time_since_entrance_end', 'rate']]
    return hazard_rate_data
