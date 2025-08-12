from typing import Dict
from dataclasses import dataclass

@dataclass()
class Bounds:
    """
    Describes the boundaries of an experiment.
    Specifically the UMI portion.
    """

    metadata:dict=None
    min_mpra_umi:float=None
    max_mpra_umi:float=None
    by_cre_theta:float=None
    by_cell_type_theta:float=None
    num_cres:int=None
    num_cell_types:int=None
    #mean_mpra_barcodes_per_cre:float
    #stdev_mpra_barcodes_per_cre:float
    #barcode_deviation:float
    #barcode deviation stores how much each barcode affects umi count
    #ehh... do we really want to model this? probably skip it. 

    @classmethod
    def from_ortho(cls):
        pass

#class wet_bounds
#   read_depth