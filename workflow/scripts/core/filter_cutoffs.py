import numpy as np

from functools import wraps
from typing import Callable, Tuple
from numpy.typing import ArrayLike


def _validate_array(func: Callable) -> Callable:
    '''Decorator for sprase normalization methods'''

    @wraps(func)
    def wrapper(arr, **kwargs):
        if not np.isscalar(arr):
            return func(arr, **kwargs)
        else:
            raise ValueError("Argument must be an array!")

    return wrapper


@_validate_array
def get_iqr_cutoffs(arr: ArrayLike, k=1.5) -> Tuple[float]:
    q1, q3 = np.percentile(arr, [25, 75])
    iqr = q3 - q1    
    lower_cutoff = q1 - k*iqr
    upper_cutoff = q3 + k*iqr

    return (lower_cutoff, upper_cutoff)
    

@_validate_array
def get_zscore_cutoffs(arr: ArrayLike, threshold: float = 1.96) -> Tuple[float]:
    std = arr.std()
    mean = arr.mean()
    lower_cutoff = mean - threshold*std
    upper_cutoff = mean + threshold*std

    return (lower_cutoff, upper_cutoff)
    

@_validate_array
def get_mad_cutoffs(arr: ArrayLike, threshold: float = 2.5) -> Tuple[float]:        
    pass