import typing
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import product

from string import ascii_uppercase

COLS_96 = range(1, 13)
ROWS_96 = ascii_uppercase[:8]
COLS_384 = range(1, 25)
ROWS_384 = ascii_uppercase[:16]
INDEX_96_C = pd.Index([f"{r}{c}" for r,c in product(ROWS_96, COLS_96)])
INDEX_96_F = pd.Index(INDEX_96_C.values.reshape(len(ROWS_96), -1).ravel(order="F"))
INDEX_384_C = pd.Index([f"{r}{c}" for r,c in product(ROWS_384, COLS_384)])
INDEX_384_F = pd.Index(INDEX_384_C.values.reshape(len(ROWS_384), -1).ravel(order="F"))

assert len(INDEX_96_C) == len(INDEX_96_F) == 96
assert len(INDEX_384_C) == len(INDEX_384_F) == 384


def well_sort(item):
    return (item[0], int(item[1:]))


def sort_index(series):
    return series.reindex(sorted(series.index, key=well_sort))


def flatten_plate_map(data, colwise=False):
    ravel = "F" if colwise else "C"
    if isinstance(data, pd.Series) or isinstance(data, pd.Index):
        data = data.values
    return data.ravel(order=ravel)


def unflatten_plate_map(data, colwise=False):
    index = INDEX_384_C if len(data) == 384 else INDEX_96_C
    dims = (list(ROWS_384), COLS_384) if len(data) == 384 else (list(ROWS_96), COLS_96)
    plate = pd.DataFrame(index=dims[0], columns=dims[1])
    plate.loc[:,:] = data.reindex(index).to_numpy().reshape(plate.shape)
    return plate


def construct_96(data_96, name, colwise=False):
    index = INDEX_96_F if colwise else INDEX_96_C
    return pd.Series(data_96, name=name, index=index)


def construct_384(data_384, name, colwise=False):
    index = INDEX_384_F if colwise else INDEX_384_C
    return pd.Series(data_384, name=name, index=index)

def index_by_quad(quad, colwise=False):
    assert quad < 4
    ravel = "F" if colwise else "C"
    source = INDEX_384_C.values.reshape(len(ROWS_384), -1)
    i, j = quad//2, quad%2
    return pd.Index(source[i::2, j::2].ravel(order=ravel))


def convert_96_to_384(data_96, name, quad=None, colwise=False):
    if quad is not None and len(data_96) > 96:
        raise ValueError(f"96 well data for quad={quad} has more than 96 values")
    elif quad is None and len(data_96) == 96:
        raise ValueError(f"96 well data with no quad specified has only 96 values")

    if quad is None:
        index = pd.Index(sum([index_by_quad(k, colwise).tolist() for k in range(4)], []))
        quads = np.repeat(range(4), 96)
    else:
        index = index_by_quad(quad, colwise)
        quads = np.ones(96, dtype=np.uint8) * (quad + 1)

    s = pd.DataFrame(data_96, index=index, columns=[name])
    s["Quadrant"] = quads

    return s


def split_row_col(s):
    if isinstance(s, list) or isinstance(s, pd.Index):
        s = pd.Series(s)
    return s.str.extract("([A-Z])(\d+)", expand=True)


def log_assert(assertion, message):
    try:
        assert assertion
    except AssertionError as err:
        logger.exception(message)
        raise err


def assert_path_exists(name: str, path: Path) -> None:
    if path is not None and not path.exists():
        logger.error(f"{name.replace('_',' ').capitalize()}: [{path}] does not exist!")
        exit(2)


def assert_path_does_not_exist(name: str, path: Path) -> None:
    if path is not None and path.exists():
        logger.error(f"{name.replace('_',' ').capitalize()}: [{path}] already exist!")
        exit(2)
