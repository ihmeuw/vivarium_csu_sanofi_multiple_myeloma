import itertools
from collections import Counter
from typing import Callable, Dict, Iterable, List, Tuple, TYPE_CHECKING

import pandas as pd
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event
from vivarium.framework.population import SimulantData
from vivarium_public_health.metrics import MortalityObserver as MortalityObserver_
from vivarium_public_health.metrics.utilities import (get_deaths, get_person_time, get_transition_count,
                                                      get_years_lived_with_disability, get_years_of_life_lost,
                                                      TransitionString)

from vivarium_csu_sanofi_multiple_myeloma.constants import models, results


class ResultsStratifier:
    """Centralized component for handling results stratification.

    This should be used as a sub-component for observers.  The observers
    can then ask this component for population subgroups and labels during
    results production and have this component manage adjustments to the
    final column labels for the subgroups.

    """

    def __init__(self, observer_name: str, stratify_by_treatment: bool = True):
        self.name = f'{observer_name}_results_stratifier'
        self.stratify_by_treatment = stratify_by_treatment

    # noinspection PyAttributeOutsideInit
    def setup(self, builder: Builder):
        """Perform this component's setup."""
        # The only thing you should request here are resources necessary for results stratification.
        self.pipelines = {}
        columns_required = [
            'age',
            'multiple_myeloma_treatment',
            'retreated'
        ]

        def get_treatment_state_function(state):
            return lambda: self.population_values['multiple_myeloma_treatment'] == state

        def get_retreatment_state_function(state):
            return lambda: self.population_values['retreated'] == state

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

    def __init__(self, stratify_by_treatment: str = 'True'):
        super().__init__()
        self.stratifier = ResultsStratifier(self.name, stratify_by_treatment == 'True')

    @property
    def sub_components(self) -> List[ResultsStratifier]:
        return [self.stratifier]

    def metrics(self, index: pd.Index, metrics: Dict[str, float]) -> Dict[str, float]:
        pop = self.population_view.get(index)
        pop.loc[pop.exit_time.isnull(), 'exit_time'] = self.clock()

        measure_getters = (
            (get_deaths, (self.causes,)),
            (get_person_time, ()),
            (get_years_of_life_lost, (self.life_expectancy, self.causes)),
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
