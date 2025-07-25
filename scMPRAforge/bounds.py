from typing import Dict
from dataclasses import dataclass

@dataclass(frozen=True)
class Bounds:
    """
    Describes the boundaries of an experiment.
    """
    metadata:str
    min_mpra_umi:float
    max_mpra_umi:float
    by_cre_theta:float
    by_cell_type_theta:float

    @classmethod
    def from_ortho(cls):
        pass

