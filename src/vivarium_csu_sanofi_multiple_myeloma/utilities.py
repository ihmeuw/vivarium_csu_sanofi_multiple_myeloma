import click
import numpy as np
import pandas as pd

from typing import Union, List
from pathlib import Path
from loguru import logger
from scipy.stats import norm, lognorm

from vivarium_public_health.risks.data_transformations import pivot_categorical

from vivarium_csu_sanofi_multiple_myeloma.constants import metadata


def len_longest_location() -> int:
    """Returns the length of the longest location in the project.

    Returns
    -------
       Length of the longest location in the project.
    """
    return len(max(metadata.LOCATIONS, key=len))


def sanitize_location(location: str):
    """Cleans up location formatting for writing and reading from file names.

    Parameters
    ----------
    location
        The unsanitized location name.

    Returns
    -------
        The sanitized location name (lower-case with white-space and
        special characters removed.

    """
    # FIXME: Should make this a reversible transformation.
    return location.replace(" ", "_").replace("'", "_").lower()


def delete_if_exists(*paths: Union[Path, List[Path]], confirm=False):
    paths = paths[0] if isinstance(paths[0], list) else paths
    existing_paths = [p for p in paths if p.exists()]
    if existing_paths:
        if confirm:
            # Assumes all paths have the same root dir
            root = existing_paths[0].parent
            names = [p.name for p in existing_paths]
            click.confirm(f"Existing files {names} found in directory {root}. Do you want to delete and replace?",
                          abort=True)
        for p in existing_paths:
            logger.info(f'Deleting artifact at {str(p)}.')
            p.unlink()


def read_data_by_draw(artifact_path: str, key : str, draw: int) -> pd.DataFrame:
    """Reads data from the artifact on a per-draw basis. This
    is necessary for Low Birthweight Short Gestation (LBWSG) data.

    Parameters
    ----------
    artifact_path
        The artifact to read from.
    key
        The entity key associated with the data to read.
    draw
        The data to retrieve.

    """
    key = key.replace(".", "/")
    with pd.HDFStore(artifact_path, mode='r') as store:
        index = store.get(f'{key}/index')
        draw = store.get(f'{key}/draw_{draw}')
    draw = draw.rename("value")
    data = pd.concat([index, draw], axis=1)
    data = data.drop(columns='location')
    data = pivot_categorical(data)
    return data


class LogNormalHazardRate:
    """Defines an instance of a lognormal hazard rate normal distribution.
    Parameters
    ----------
    hr
        mean of distribution
    hr_lower
        lower bound the lognormal distribution
    hr_upper
        upper bound of lognormal distribution
    Returns
    -------
        An object with parameters for scipy.stats.lognorm
    """

    def __init__(self, hr: float, hr_lower: float, hr_upper: float):
        self.hr = hr
        self.hr_upper = hr_upper
        q_975_stdnorm = norm().ppf(0.975)
        mu = np.log(self.hr)
        sigma = (np.log(self.hr_upper) - mu) / q_975_stdnorm
        self.hr_distribution = lognorm(s=sigma, scale=self.hr)

    def get_random_variable(self, percentile: float) -> float:
        """Gets a single random draw from a log normal hazard rate distribution.
        Parameters
        ----------
        percentile
            Percentile for sample.
        Returns
        -------
            The random variate from the log normal hazard rate distribution.
        """
        return self.hr_distribution.ppf(percentile)


