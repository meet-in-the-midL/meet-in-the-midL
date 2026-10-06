from enum import *
import os
from pathlib import Path as OSPath
from sage.all import *

from maths_utils import *

# ============================================================================ #

load : Callable[[Any], Any]
RELATIVE_DATAPATH = "data"

class DataSets(Enum):
    GRAPH = auto()
    EIGENVALUE = auto()
    PATHS = auto()
    MITM = auto()

def data_path(
        p : Prime,
        l_or_L : Prime | list[Prime],
        datatype : DataSets,
        *args : Any,
    ) -> str:
    sub_path = f"/".join(map(str, args))
    if datatype == DataSets.GRAPH:
        assert isinstance(l_or_L, (int, Integer))
        return os.path.join(os.getcwd(), RELATIVE_DATAPATH, "graphs", sub_path, f"{p}/{l_or_L}.sobj")
    assert isinstance(l_or_L, list)
    l_or_L.sort()
    l_str = L_id(l_or_L)
    match datatype:
        case DataSets.EIGENVALUE:
            return os.path.join(os.getcwd(), RELATIVE_DATAPATH, "eigs", sub_path, f"{p}/{l_str}.sobj")
        case DataSets.PATHS:
            return os.path.join(os.getcwd(), RELATIVE_DATAPATH, "paths", sub_path, f"{p}/{l_str}.sobj")
        case DataSets.MITM:
            return os.path.join(os.getcwd(), RELATIVE_DATAPATH, "mitm", sub_path, f"{p}/{l_str}.sobj")
        
def fetch_data(
        p : Prime,
        l_or_L : Prime | list[Prime],
        datatype : DataSets,
        *args : Any,
    ) -> Any:
    path = data_path(p, l_or_L, datatype, *args)
    out = load(path)
    return out

def store_data(
        p : Prime,
        l_or_L : Prime | list[Prime],
        datatype : DataSets,
        data : Any,
        *args : Any,
    ) -> None:
    path = data_path(p, l_or_L, datatype, *args)
    OSPath(path).parent.mkdir(parents = True, exist_ok = True)
    save(data, path)

