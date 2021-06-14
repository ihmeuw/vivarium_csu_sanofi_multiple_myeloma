from collections import Counter
from typing import Dict, TYPE_CHECKING

import pandas as pd

from vivarium_csu_sanofi_multiple_myeloma.components.treatment import TREATMENTS
from vivarium_csu_sanofi_multiple_myeloma.constants import models

if TYPE_CHECKING:
    from vivarium.framework.engine import Builder
    from vivarium.framework.event import Event


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
            for t in TREATMENTS:
                key = treatment_template.format(treatment_line=treatment_line, treatment=t, year=event.time.year)
                got_treatment = pop['multiple_myeloma_treatment'] == t
                counts[key] = (had_mm_event & got_treatment).sum()
        self.counts.update(counts)

    def metrics(self, index: pd.Index, metrics: Dict[str, float]) -> Dict[str, float]:
        metrics.update(self.counts)
        return metrics

    def __repr__(self):
        return f'{self.__class__.__name__}()'
