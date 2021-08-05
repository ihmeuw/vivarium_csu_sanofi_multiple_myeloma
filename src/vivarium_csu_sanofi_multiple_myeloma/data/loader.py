"""Loads, standardizes and validates input data for the simulation.

Abstract the extract and transform pieces of the artifact ETL.
The intent here is to provide a uniform interface around this portion
of artifact creation. The value of this interface shows up when more
complicated data needs are part of the project. See the BEP project
for an example.

`BEP <https://github.com/ihmeuw/vivarium_gates_bep/blob/master/src/vivarium_gates_bep/data/loader.py>`_

.. admonition::

   No logging is done here. Logging is done in vivarium inputs itself and forwarded.
"""
import pandas as pd

from gbd_mapping import causes, covariates, risk_factors
from vivarium.framework.artifact import EntityKey
from vivarium_gbd_access import gbd
from vivarium_inputs import globals as vi_globals, interface, utilities as vi_utils, utility_data
from vivarium_inputs.mapping_extension import alternative_risk_factors

from vivarium_csu_sanofi_multiple_myeloma.constants import data_keys, models


def get_data(lookup_key: str, location: str) -> pd.DataFrame:
    """Retrieves data from an appropriate source.

    Parameters
    ----------
    lookup_key
        The key that will eventually get put in the artifact with
        the requested data.
    location
        The location to get data for.

    Returns
    -------
        The requested data.

    """
    mapping = {
        data_keys.POPULATION.LOCATION: load_population_location,
        data_keys.POPULATION.STRUCTURE: load_population_structure,
        data_keys.POPULATION.AGE_BINS: load_age_bins,
        data_keys.POPULATION.DEMOGRAPHY: load_demographic_dimensions,
        data_keys.POPULATION.TMRLE: load_theoretical_minimum_risk_life_expectancy,
        data_keys.POPULATION.ACMR: load_standard_data,

        data_keys.MULTIPLE_MYELOMA.PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_1_STATE_NAME),
        data_keys.MULTIPLE_MYELOMA.INCIDENCE_RATE: load_standard_data,
        data_keys.MULTIPLE_MYELOMA.DISABILITY_WEIGHT: load_rrmm_disability_weight,
        data_keys.MULTIPLE_MYELOMA.CSMR: load_acmr,
        data_keys.MULTIPLE_MYELOMA.RESTRICTIONS: load_metadata,

        data_keys.MULTIPLE_MYELOMA.LINE_1_PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_1_STATE_NAME),
        data_keys.MULTIPLE_MYELOMA.LINE_2_PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_2_STATE_NAME),
        data_keys.MULTIPLE_MYELOMA.LINE_3_PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_3_STATE_NAME),
        data_keys.MULTIPLE_MYELOMA.LINE_4_PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_4_STATE_NAME),
        data_keys.MULTIPLE_MYELOMA.LINE_5_PREVALENCE: load_rrmm_prevalence(models.MULTIPLE_MYELOMA_5_STATE_NAME),

        data_keys.MULTIPLE_MYELOMA.LINE_1_DISABILITY_WEIGHT: load_rrmm_disability_weight,
        data_keys.MULTIPLE_MYELOMA.LINE_2_DISABILITY_WEIGHT: load_rrmm_disability_weight,
        data_keys.MULTIPLE_MYELOMA.LINE_3_DISABILITY_WEIGHT: load_rrmm_disability_weight,
        data_keys.MULTIPLE_MYELOMA.LINE_4_DISABILITY_WEIGHT: load_rrmm_disability_weight,
        data_keys.MULTIPLE_MYELOMA.LINE_5_DISABILITY_WEIGHT: load_rrmm_disability_weight,

        data_keys.MULTIPLE_MYELOMA.LINE_1_EMR: load_rrmm_emr,
        data_keys.MULTIPLE_MYELOMA.LINE_2_EMR: load_rrmm_emr,
        data_keys.MULTIPLE_MYELOMA.LINE_3_EMR: load_rrmm_emr,
        data_keys.MULTIPLE_MYELOMA.LINE_4_EMR: load_rrmm_emr,
        data_keys.MULTIPLE_MYELOMA.LINE_5_EMR: load_rrmm_emr,

        data_keys.MULTIPLE_MYELOMA.GBD_CSMR: load_gbd_csmr,
        data_keys.MULTIPLE_MYELOMA.SUSCEPTIBLE_EMR: load_susceptible_emr,
    }
    return mapping[lookup_key](lookup_key, location)


def load_population_location(key: str, location: str) -> str:
    if key != data_keys.POPULATION.LOCATION:
        raise ValueError(f'Unrecognized key {key}')
    return location


def load_population_structure(key: str, location: str) -> pd.DataFrame:
    return interface.get_population_structure(location)


def load_age_bins(key: str, location: str) -> pd.DataFrame:
    return interface.get_age_bins()


def load_demographic_dimensions(key: str, location: str) -> pd.DataFrame:
    return interface.get_demographic_dimensions(location)


def load_theoretical_minimum_risk_life_expectancy(key: str, location: str) -> pd.DataFrame:
    return interface.get_theoretical_minimum_risk_life_expectancy()


def load_standard_data(key: str, location: str) -> pd.DataFrame:
    key = EntityKey(key)
    entity = get_entity(key)
    return interface.get_measure(entity, key.measure, location).droplevel('location')


def load_metadata(key: str, location: str):
    key = EntityKey(key)
    entity = get_entity(key)
    entity_metadata = entity[key.measure]
    if hasattr(entity_metadata, 'to_dict'):
        entity_metadata = entity_metadata.to_dict()
    return entity_metadata


def _load_em_from_meid(location, meid, measure):
    location_id = utility_data.get_location_id(location)
    data = gbd.get_modelable_entity_draws(meid, location_id)
    data = data[data.measure_id == vi_globals.MEASURES[measure]]
    data = vi_utils.normalize(data, fill_value=0)
    data = data.filter(vi_globals.DEMOGRAPHIC_COLUMNS + vi_globals.DRAW_COLUMNS)
    data = vi_utils.reshape(data)
    data = vi_utils.scrub_gbd_conventions(data, location)
    data = vi_utils.split_interval(data, interval_column='age', split_column_prefix='age')
    data = vi_utils.split_interval(data, interval_column='year', split_column_prefix='year')
    return vi_utils.sort_hierarchical_data(data)


# Project-specific data functions
def load_acmr(key: str, location) -> pd.DataFrame:
    return load_standard_data('cause.all_causes.cause_specific_mortality_rate', location)


def load_gbd_csmr(key: str, location) -> pd.DataFrame:
    return load_standard_data('cause.multiple_myeloma.cause_specific_mortality_rate', location)


def load_rrmm_prevalence(state: str) -> pd.DataFrame:
    state_multipliers = {
        models.MULTIPLE_MYELOMA_1_STATE_NAME: 1.,
        models.MULTIPLE_MYELOMA_2_STATE_NAME: 0.,
        models.MULTIPLE_MYELOMA_3_STATE_NAME: 0.,
        models.MULTIPLE_MYELOMA_4_STATE_NAME: 0.,
        models.MULTIPLE_MYELOMA_5_STATE_NAME: 0.,
    }

    def _load_rrmm_prevalence(key, location) -> pd.DataFrame:
        return load_standard_data('cause.multiple_myeloma.prevalence', location) * state_multipliers[state]
    return _load_rrmm_prevalence


def load_rrmm_disability_weight(key: str, location):
    return load_standard_data('cause.multiple_myeloma.disability_weight', location)


def load_rrmm_emr(key: str, location):
    return load_standard_data('cause.multiple_myeloma.excess_mortality_rate', location)


def load_susceptible_emr(key: str, location):
    INDEX_COLUMNS = ['sex', 'age_start', 'age_end', 'year_start', 'year_end']
    return (load_acmr(key, location).reset_index().set_index(INDEX_COLUMNS)
            - load_gbd_csmr(key, location).reset_index().set_index(INDEX_COLUMNS))


def get_entity(key: str):
    # Map of entity types to their gbd mappings.
    type_map = {
        'cause': causes,
        'covariate': covariates,
        'risk_factor': risk_factors,
        'alternative_risk_factor': alternative_risk_factors
    }
    key = EntityKey(key)
    return type_map[key.type][key.name]
