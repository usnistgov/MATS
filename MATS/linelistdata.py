"""
Routines to load example data locally or from web
"""

from functools import partial
from pathlib import Path
import pandas as pd


# def _attribute_getter(self, name):
#     return self[name]


class LoadLineListData(object):
    """
    Helper class to read in supplied LineList DataFrames


    Attributes
    ----------
    names : list
        names of LineList files available

    Methods
    -------
    __getitem__ : get file DataFrame by position/name
        `self[x]` returns a copy of the linelist frame corresponding to
        `self.names[x]` if `x` is an integer, otherwise returns Frame corresponding to `name = x`
    """
    # _names = [
    #     "O2_ABand_Drouin_2017_linelist",
    #     "JQSRT2021_SDNGP_2015",
    #     "CO2_30012",
    # ]

    # _prefix_web = (
    #     "https://raw.githubusercontent.com/usnistgov/MATS/master/MATS/Linelists/"
    # )
    _prefix_local = Path(__file__).parent / "Linelists"

    def __init__(self, paths=None):

        paths_default = list(self._prefix_local.glob('*.csv'))
        if paths is None:
            paths = []
        else:
            paths = [Path(p) for p in paths]
        paths = paths_default + paths

        for p in paths:
            if not p.exists():
                raise ValueError(f'{p} does not exist')

        self._paths_dict = {p.with_suffix('').name : p for p in paths}
        self._cache = {}

    def _ipython_key_completions_(self):
        return self.names

    @property
    def names(self):
        return list(self._paths_dict.keys())

    @property
    def paths(self):
        return self._paths_dict

    def _get_file(self, name):
        if name not in self._paths_dict:
            raise ValueError("file name must be in {}".format(self.names))
        path = self._paths_dict[name]
        if name not in self._cache:
            self._cache[name] = pd.read_csv(path)

        # always return a copy.
        # MATS manipulates in place the linelist
        return self._cache[name].copy()

    def __getitem__(self, index):
        if isinstance(index, int):
            name = self.names[index]
        elif isinstance(index, str):
            name = index
        else:
            raise ValueError("bad index {}".format(index))
        return self._get_file(name)


# Global example data loader
linelistdata = LoadLineListData()
