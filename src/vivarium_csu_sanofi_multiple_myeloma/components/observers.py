import itertools
from collections import Counter

from typing import Callable, Dict, Iterable, List, Tuple

import pandas as pd
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event
from vivarium.framework.population import SimulantData

from vivarium_public_health.metrics import MortalityObserver as MortalityObserver_, DiseaseObserver as DiseaseObserver_
from vivarium_public_health.metrics.utilities import (get_age_bins, get_deaths, get_person_time, get_state_person_time,
                                                      get_transition_count, get_group_counts, get_output_template,
                                                      get_years_of_life_lost, TransitionString)

from vivarium_csu_sanofi_multiple_myeloma.constants import data_values, models, results
from vivarium_csu_sanofi_multiple_myeloma.constants.data_values import RISKS, RISK_LEVEL_MAP, UNDIAGNOSED


class ResultsStratifier:
    """Centralized component for handling results stratification.

    This should be used as a sub-component for observers.  The observers
    can then ask this component for population subgroups and labels during
    results production and have this component manage adjustments to the
    final column labels for the subgroups.

    """

    def __init__(self, observer_name: str, stratify_by_treatment: bool = True, stratify_by_risks: bool = True):
        self.name = f'{observer_name}_results_stratifier'
        self.stratify_by_treatment = stratify_by_treatment
        self.stratify_by_risks = stratify_by_risks

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: Builder):
        """Perform this component's setup."""
        # The only thing you should request here are resources necessary for results stratification.
        self.pipelines = {}
        columns_required = [
            'age',
            'multiple_myeloma_treatment',
            'retreated',
            'race_and_cytogenetic_risk_at_diagnosis',
            'renal_function_at_diagnosis'
        ]

        def get_treatment_state_function(state):
            return lambda: self.population_values['multiple_myeloma_treatment'] == state

        def get_retreatment_state_function(state):
            return lambda: self.population_values['retreated'] == state

        def get_race_and_cytogenetic_risk_at_diagnosis_function(state):
            return lambda: self.population_values['race_and_cytogenetic_risk_at_diagnosis'] == state

        def get_renal_function_at_diagnosis(state):
            return lambda: self.population_values['renal_function_at_diagnosis'] == state

        self.stratification_levels = {}

        if self.stratify_by_treatment:
            self.stratification_levels['treatment_state'] = {
                treatment_state: get_treatment_state_function(treatment_state)
                for treatment_state in models.TREATMENTS
            }
            self.stratification_levels['retreated'] = {
                retreatment_state: get_retreatment_state_function(retreatment_state)
                for retreatment_state in (True, False)
            }

        if self.stratify_by_risks:
            self.stratification_levels['race_and_cytogenetic_risk_at_diagnosis'] = {
                race_and_cytogenetic_risk_at_diagnosis: get_race_and_cytogenetic_risk_at_diagnosis_function(
                    race_and_cytogenetic_risk_at_diagnosis)
                for race_and_cytogenetic_risk_at_diagnosis
                in RISK_LEVEL_MAP[RISKS.race_and_cytogenetic_risk_at_diagnosis] + [UNDIAGNOSED]
            }
            self.stratification_levels['renal_function_at_diagnosis'] = {
                renal_function_at_diagnosis: get_renal_function_at_diagnosis(renal_function_at_diagnosis)
                for renal_function_at_diagnosis in RISK_LEVEL_MAP[RISKS.renal_function_at_diagnosis] + [UNDIAGNOSED]
            }

        self.population_view = builder.population.get_view(columns_required)
        self.pipeline_values = {pipeline: None for pipeline in self.pipelines}
        self.population_values = None
        self.stratification_groups = None

        builder.population.initializes_simulants(self.on_initialize_simulants,
                                                 requires_columns=columns_required,
                                                 requires_values=list(self.pipelines.keys()))

        builder.event.register_listener('time_step__prepare', self.on_timestep_prepare)

    # noinspection PyAttributeOutsideInit
    def on_initialize_simulants(self, pop_data: SimulantData):
        self.pipeline_values = {name: pipeline(pop_data.index) for name, pipeline in self.pipelines.items()}
        self.population_values = self.population_view.get(pop_data.index)
        self.stratification_groups = self.get_stratification_groups(pop_data.index)

    def get_stratification_groups(self, index):
        stratification_groups = pd.Series('', index=index)
        all_stratifications = self.get_all_stratifications()
        for stratification in all_stratifications:
            stratification_group_name = '_'.join([f'{metric["metric"]}_{metric["category"]}'
                                                  for metric in stratification])
            mask = pd.Series(True, index=index)
            for metric in stratification:
                mask &= self.stratification_levels[metric['metric']][metric['category']]()
            stratification_groups.loc[mask] = stratification_group_name
        return stratification_groups

    def append_new_entrants(self, existing_data: pd.Series, new_index: pd.Index, getter: Callable):
        intersection = existing_data.loc[new_index.intersection(existing_data.index)]

        new_entrants_index = new_index.difference(self.stratification_groups.index)
        new_entrants_stratifications = getter(new_entrants_index)
        return intersection.append(new_entrants_stratifications)

    # noinspection PyAttributeOutsideInit
    def on_timestep_prepare(self, event: Event):
        self.pipeline_values = {name: pipeline(event.index) for name, pipeline in self.pipelines.items()}
        self.population_values = self.population_view.get(event.index)
        self.stratification_groups = self.get_stratification_groups(event.index)

    def get_all_stratifications(self) -> List[Tuple[Dict[str, str], ...]]:
        """
        Gets all stratification combinations. Returns a List of Stratifications. Each Stratification is represented as a
        Tuple of Stratification Levels. Each Stratification Level is represented as a Dictionary with keys 'metric' and
        'category'. 'metric' refers to the stratification level's name, and 'category' refers to the stratification
        category.

        If no stratification levels are defined, returns a List with a single empty Tuple
        """
        # Get list of lists of metric and category pairs for each metric
        groups = [[{'metric': metric, 'category': category} for category, category_mask in category_maps.items()]
                  for metric, category_maps in self.stratification_levels.items()]
        # Get product of all stratification combinations
        return list(itertools.product(*groups))

    @staticmethod
    def get_stratification_key(stratification: Iterable[Dict[str, str]]) -> str:
        return ('' if not stratification
                else '_'.join([f'{metric["metric"]}_{metric["category"]}' for metric in stratification]))

    def group(self, pop: pd.DataFrame) -> Iterable[Tuple[Tuple[str, ...], pd.DataFrame]]:
        """Takes the full population and yields stratified subgroups.

        Parameters
        ----------
        pop
            The population to stratify.

        Yields
        ------
            A tuple of stratification labels and the population subgroup
            corresponding to those labels.

        """
        stratification_groups = self.append_new_entrants(self.stratification_groups, pop.index,
                                                         self.get_stratification_groups)
        stratifications = self.get_all_stratifications()
        for stratification in stratifications:
            stratification_key = self.get_stratification_key(stratification)
            if pop.empty:
                pop_in_group = pop
            else:
                pop_in_group = pop.loc[(stratification_groups == stratification_key)]
            yield (stratification_key,), pop_in_group

    @staticmethod
    def update_labels(measure_data: Dict[str, float], labels: Tuple[str, ...]) -> Dict[str, float]:
        """Updates a dict of measure data with stratification labels.

        Parameters
        ----------
        measure_data
            The measure data with unstratified column names.
        labels
            The stratification labels. Yielded along with the population
            subgroup the measure data was produced from by a call to
            :obj:`ResultsStratifier.group`.

        Returns
        -------
            The measure data with column names updated with the stratification
            labels.

        """
        stratification_label = f'_{labels[0]}' if labels[0] else ''
        measure_data = {f'{k}{stratification_label}': v for k, v in measure_data.items()}
        return measure_data


class MortalityObserver(MortalityObserver_):

    def __init__(self, stratify_by_treatment: str = 'True', stratify_by_risks: str = 'True'):
        super().__init__()
        self.stratifier = ResultsStratifier(self.name, stratify_by_treatment == 'True', stratify_by_risks == 'True')

    @property
    def sub_components(self) -> List[ResultsStratifier]:
        return [self.stratifier]

    def metrics(self, index: pd.Index, metrics: Dict[str, float]) -> Dict[str, float]:
        pop = self.population_view.get(index)
        pop.loc[pop.exit_time.isnull(), 'exit_time'] = self.clock()

        measure_getters = (
            (get_deaths, (self.causes + [models.SUSCEPTIBLE_STATE_NAME],)),
            (get_person_time, ()),
            (get_years_of_life_lost, (self.life_expectancy, self.causes + [models.SUSCEPTIBLE_STATE_NAME])),
        )

        for labels, pop_in_group in self.stratifier.group(pop):
            base_args = (pop_in_group, self.config.to_dict(), self.start_time, self.clock(), self.age_bins)

            for measure_getter, extra_args in measure_getters:
                measure_data = measure_getter(*base_args, *extra_args)
                measure_data = self.stratifier.update_labels(measure_data, labels)
                metrics.update(measure_data)

        the_living = pop[(pop.alive == 'alive') & pop.tracked]
        the_dead = pop[pop.alive == 'dead']
        metrics[results.TOTAL_YLLS_COLUMN] = self.life_expectancy(the_dead.index).sum()
        metrics['total_population_living'] = len(the_living)
        metrics['total_population_dead'] = len(the_dead)

        return metrics


class MultipleMyelomaTreatmentObserver:

    @property
    def name(self) -> str:
        return self.__class__.__name__

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: 'Builder') -> None:
        self.counts = Counter()
        self.population_view = builder.population.get_view(
            [f'{s}_event_time' for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES]
            + ['multiple_myeloma_treatment']
        )
        builder.value.register_value_modifier('metrics', self.metrics)
        builder.event.register_listener('collect_metrics', self.on_collect_metrics)

    def on_collect_metrics(self, event: 'Event') -> None:
        treatment_template = 'line_{treatment_line}_treatment_{treatment}_year_{year}'
        pop = self.population_view.get(event.index)
        counts = {}
        for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES:
            treatment_line = s.split('_')[-1]
            had_mm_event = pop[f'{s}_event_time'] == event.time
            for t in models.TREATMENTS:
                key = treatment_template.format(treatment_line=treatment_line, treatment=t, year=event.time.year)
                got_treatment = pop['multiple_myeloma_treatment'] == t
                counts[key] = (had_mm_event & got_treatment).sum()
        self.counts.update(counts)

    def metrics(self, index: pd.Index, metrics: Dict[str, float]) -> Dict[str, float]:
        metrics.update(self.counts)
        return metrics

    def __repr__(self):
        return f'{self.__class__.__name__}()'


class SurvivalObserver:

    configuration_defaults = {
        'metrics': {
            'observation_start': {
                'year': 2021,
                'month': 1,
                'day': 1,
            }
        }
    }

    @property
    def name(self) -> str:
        return self.__class__.__name__

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: 'Builder') -> None:
        config = builder.configuration.metrics
        self.observation_start = pd.Timestamp(**config.observation_start)
        # Pull directly from config rather than the more flexible
        # builder.time.step_size(). This component depends on a
        # constant step size.
        self.step_size = builder.configuration.time.step_size
        self.bins = (pd.interval_range(0, 60*28, freq=self.step_size)
                     .append(pd.Index([pd.Interval(60*28, 1000*28)])))

        self.counts = Counter()
        risk_template = '_'.join([f'{s}_{{{s}}}' for s in data_values.RISKS])
        count_template = 'alive_at_period_{period_start}_line_{treatment_line}' + '_' + risk_template
        progression_template = 'progressed_by_period_{period_end}_line_{treatment_line}' + '_' + risk_template
        death_template = 'died_by_period_{period_end}_line_{treatment_line}' + '_' + risk_template
        self.sim_end_template = 'sim_end_on_period_{period_end}_line_{treatment_line}' + '_' + risk_template
        self.templates = [
            ('alive', count_template),
            ('progressed', progression_template),
            ('died', death_template),
        ]

        self.population_view = builder.population.get_view(
            [f'{s}_event_time' for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES]
            + [f'{s}_time_since_entrance' for s in models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES]
            + [models.MULTIPLE_MYELOMA_MODEL_NAME, 'alive', 'exit_time']
            + list(data_values.RISKS)
        )

        builder.event.register_listener('collect_metrics', self.on_collect_metrics)
        builder.event.register_listener('simulation_end', self.on_simulation_end)

        builder.value.register_value_modifier('metrics', self.metrics)

    def on_collect_metrics(self, event: 'Event') -> None:
        pop = self.get_denominator_pop(event)
        states = list(models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES)

        for current_state, next_state in zip(states, states[1:] + [states[-1]]):
            state_denominator = self.subset_state_denominator(pop, current_state, next_state, event)
            if not state_denominator.empty:
                print(f'current state: {current_state}, state denominator size: {len(state_denominator)}')
            for risk_status in itertools.product(*data_values.RISK_LEVEL_MAP.values()):
                denominator = self.subset_risk_denominator(state_denominator, risk_status)

                alive_at_start = (denominator
                                  .groupby('group')
                                  .multiple_myeloma
                                  .count()
                                  .rename('alive'))
                died_by_end = (denominator[denominator['exit_time'] == event.time]
                               .groupby('group')
                               .multiple_myeloma
                               .count()
                               .rename('died'))
                progressed_by_end = (denominator[denominator[f'{next_state}_event_time'] == event.time]
                                     .groupby('group')
                                     .multiple_myeloma
                                     .count()
                                     .rename('progressed'))
                survival_results = pd.concat([alive_at_start, died_by_end, progressed_by_end], axis=1)
                survival_results.index = pd.IntervalIndex(survival_results.index)
                if not denominator.empty:
                    print(f'risk_status: {risk_status}, denominator size: {len(denominator)}, summary: {survival_results.sum().to_dict()}')
                treatment_line = current_state.split('_')[-1]
                for interval, interval_data in survival_results.iterrows():
                    for measure, template in self.templates:
                        key = template.format(
                            treatment_line=treatment_line,
                            period_start=interval.left,
                            period_end=interval.right,
                            **dict(zip(data_values.RISKS, risk_status))
                        )
                        self.counts[key] += interval_data.loc[measure]

    def on_simulation_end(self, event: 'Event'):
        pop = self.get_denominator_pop(event)
        states = list(models.MULTIPLE_MYELOMA_WITH_CONDITION_STATES)

        for current_state, next_state in zip(states, states[1:] + [states[-1]]):
            for risk_status in itertools.product(*data_values.RISK_LEVEL_MAP.values()):
                denominator = self.subset_risk_denominator(pop, risk_status)
                denominator = self.subset_state_denominator(denominator, current_state, next_state, event)
                right_censored_mask = ~(
                    (denominator['exit_time'] != event.time)
                    | (denominator[f'{next_state}_event_time'] == event.time)
                )
                right_censored = (denominator[right_censored_mask]
                                  .groupby('group')
                                  .multiple_myeloma
                                  .count()
                                  .rename('right_censored'))
                treatment_line = current_state.split('_')[-1]
                for interval, count in right_censored.iteritems():
                    key = self.sim_end_template.format(
                        treatment_line=treatment_line,
                        period_end=interval.right,
                        **dict(zip(data_values.RISKS, risk_status))
                    )
                    self.counts[key] += count

    def metrics(self, index: pd.Index, metrics: Dict[str, float]) -> Dict[str, float]:
        metrics.update(self.counts)
        return metrics

    def get_denominator_pop(self, event: 'Event'):
        pop = self.population_view.get(event.index)
        living = pop['alive'] == 'alive'
        died_this_step = (pop['alive'] == 'dead') & (pop['exit_time'] == event.time)
        in_denominator = living | died_this_step
        return pop.loc[in_denominator]

    def subset_risk_denominator(self, pop: pd.DataFrame, risk_status):
        risk_mask = pd.Series(True, index=pop.index)
        for risk, risk_level in zip(data_values.RISKS, risk_status):
            risk_mask &= pop[risk] == risk_level
        return pop.loc[risk_mask]

    def subset_state_denominator(self, pop: pd.DataFrame, current_state: str, next_state: str, event: 'Event'):
        left_censored = (pop[f'{current_state}_event_time'].notnull()
                         & (pop[f'{current_state}_event_time'] < self.observation_start))
        in_current_state_denominator = ~left_censored & (
            # In the current state and didn't get there this time step
            ((pop[models.MULTIPLE_MYELOMA_MODEL_NAME] == current_state)
             & (pop[f'{current_state}_event_time'] < event.time))
            # or in the next state, but not til the start of the next step.
            | ((pop[models.MULTIPLE_MYELOMA_MODEL_NAME] == next_state)
               & (pop[f'{next_state}_event_time'] == event.time))
        )

        denominator = pop.loc[in_current_state_denominator].copy()
        denominator['group'] = pd.cut(denominator[f'{current_state}_time_since_entrance'], self.bins)
        return denominator

                               
class DiseaseObserver(DiseaseObserver_):

    def __init__(self, disease: str, stratify_by_treatment: str = 'True', stratify_by_risks: str = 'True'):
        super().__init__(disease)
        self.stratifier = ResultsStratifier(self.name, stratify_by_treatment == 'True', stratify_by_risks == 'True')

    @property
    def sub_components(self) -> List[ResultsStratifier]:
        return [self.stratifier]

    @property
    def name(self):
        return f'{self.disease}_disease_observer'

    def on_time_step_prepare(self, event: Event):
        pop = self.population_view.get(event.index)
        # Ignoring the edge case where the step spans a new year.
        # Accrue all counts and time to the current year.
        for labels, pop_in_group in self.stratifier.group(pop):
            for state in self.states:
                # noinspection PyTypeChecker
                state_person_time_this_step = get_state_person_time(pop_in_group, self.config, self.disease, state,
                                                                    self.clock().year, event.step_size, self.age_bins)
                state_person_time_this_step = self.stratifier.update_labels(state_person_time_this_step, labels)
                self.person_time.update(state_person_time_this_step)

        # This enables tracking of transitions between states
        prior_state_pop = self.population_view.get(event.index)
        prior_state_pop[self.previous_state_column] = prior_state_pop[self.disease]
        self.population_view.update(prior_state_pop)

    def on_collect_metrics(self, event: Event):
        pop = self.population_view.get(event.index)
        for labels, pop_in_group in self.stratifier.group(pop):
            for transition in self.transitions:
                transition = TransitionString(transition)
                # noinspection PyTypeChecker
                transition_counts_this_step = get_transition_count(pop_in_group, self.config, self.disease, transition,
                                                                   event.time, self.age_bins)
                transition_counts_this_step = self.stratifier.update_labels(transition_counts_this_step, labels)
                self.counts.update(transition_counts_this_step)

    def __repr__(self) -> str:
        return f"DiseaseObserver({self.disease})"


class RegistryObserver:
    @property
    def name(self):
        return self.__class__.__name__

    # noinspection PyAttributeOutsideInit
    def setup(self, builder):
        self.clock = builder.time.clock()
        self.age_bins = get_age_bins(builder)

        required_columns = [
            'registry_evaluation_status',
            'registry_evaluation_date',
            'age',
            'sex',
            'alive',
            *data_values.RISKS
        ]

        self.population_view = builder.population.get_view(required_columns)

        builder.value.register_value_modifier('metrics', self.metrics)
        builder.event.register_listener('collect_metrics', self.on_collect_metrics)

        self.registry_entrances = Counter()
        self.registry_prevalence = Counter()

    def on_collect_metrics(self, event: Event):
        # count by stratification groups
        pop = self.population_view.get(event.index)
        self.registry_entrances.update(get_entrances(pop, self.age_bins, event.time, 'eligible'))
        self.registry_entrances.update(get_entrances(pop, self.age_bins, event.time, 'enrolled'))
        self.registry_prevalence.update(get_registry_population(pop, self.age_bins, event.time, event.step_size.days/365))

    def metrics(self, index, metrics):
        metrics.update(self.registry_entrances)
        metrics.update(self.registry_prevalence)
        return metrics

def get_entrances(pop: pd.DataFrame, age_bins: pd.DataFrame, event_time: pd.Timestamp,
                  evaluation_status: str) -> Dict[str, int]:

    mask = ((pop.alive == 'alive') & (pop.registry_evaluation_date == event_time)
            & (pop.registry_evaluation_status == evaluation_status))
    config = {'by_age': True, 'by_sex': True, 'by_year': True}
    result = {}
    risk_template = '_'.join([f'{s}_{{}}' for s in [
        data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis, data_values.RISKS.renal_function_at_diagnosis]])

    for rcr, rf in itertools.product(data_values.RISK_LEVEL_MAP[data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis],
                                            data_values.RISK_LEVEL_MAP[data_values.RISKS.renal_function_at_diagnosis]):
        risk_level_mask = ((pop[data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis] == rcr)
                           & (pop[data_values.RISKS.renal_function_at_diagnosis] == rf) & mask)
        risk_level_pop = pop.loc[risk_level_mask]
        base_key = get_output_template(**config).substitute(measure=f'registry_status_newly_{evaluation_status}', year=event_time.year)
        group_counts = get_group_counts(risk_level_pop, "", base_key, config, age_bins)
        group_counts = {key + '_' + risk_template.format(rcr, rf): value for key, value in group_counts.items()}
        result.update(group_counts)

    return result


def get_registry_population(pop: pd.DataFrame, age_bins: pd.DataFrame, event_time: pd.Timestamp, scale: float) -> Dict[str, int]:
    mask = (pop.alive == 'alive') & (pop.registry_evaluation_status == 'enrolled')
    config = {'by_age': True, 'by_sex': True, 'by_year': True}
    result = {}
    risk_template = '_'.join([f'{s}_{{}}' for s in [
        data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis, data_values.RISKS.renal_function_at_diagnosis]])
    for rcr, rf in itertools.product(
            data_values.RISK_LEVEL_MAP[data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis],
            data_values.RISK_LEVEL_MAP[data_values.RISKS.renal_function_at_diagnosis]):
        risk_level_mask = ((pop[data_values.RISKS.race_and_cytogenetic_risk_at_diagnosis] == rcr)
                           & (pop[data_values.RISKS.renal_function_at_diagnosis] == rf) & mask)
        risk_level_pop = pop.loc[risk_level_mask]
        measure = 'registry_status_enrolled'
        base_key = get_output_template(**config).substitute(measure=measure, year=event_time.year)
        group_counts = get_group_counts(risk_level_pop, "", base_key, config, age_bins, lambda x: len(x)*scale)
        group_counts = {key + '_' + risk_template.format(rcr, rf): value for key, value in group_counts.items()}
        result.update(group_counts)

    return result
