from typing import Dict
from dataclasses import dataclass
import numpy as np

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
    cells_per_cell_type:dict=None
    transfection_nb_mu:float=None
    transfection_nb_alpha:float=None

    #moi will be multiplicity of infection: number of mprabc per cell
    


    #mean_mpra_barcodes_per_cre:float
    #stdev_mpra_barcodes_per_cre:float
    #barcode_deviation:float
    #barcode deviation stores how much each barcode affects umi count
    #ehh... do we really want to model this? probably skip it. 

    @classmethod
    def from_ortho(cls,inp):
        """
        Takes an ortho object and abstracts out its bounds.

        This is an aggregation function and will hang if ortho is not done fitting yet. 
        
        This function requires that the ortho still have its training data 
        so we can extract things like "number of cells per cell-type" and 
        "number of MPRA barcodes per cell".
        """
        
        ret=cls()

        #get the min & max nb parameters across all models
        mins=[]
        maxes=[]
        
        #for each model class
        for var in ["by_cell_type_parameters","by_cre_parameters"]:
            ## nb min & max ##
            for key in getattr(inp,var).nb:
                current=getattr(inp,var).nb[key].result()
                maxes.append(float(current.max()))
                mins.append(float(current.min()))
            
            thetas=[]
            for key in getattr(inp,var).theta:
                current=getattr(inp,var).theta[key].result()
                thetas.append(current)
            
            ## theta means ##
            #we could munge the strings & use setattr but i think this is more readable
            if var=="by_cell_type_parameters":
                ret.by_cell_type_theta=np.mean(thetas)
            elif var=="by_cre_parameters":
                ret.by_cre_theta=np.mean(thetas)

        
        #min & max nb and add to the return object
        ret.min_mpra_umi=min(mins)
        ret.max_mpra_umi=max(maxes)

        ret.theta=np.mean([ret.by_cell_type_theta,ret.by_cre_theta])

        ## Parameters from training data ##    
        tfection_params=inp.training_data.describe_transfection()
        ret.transfection_nb_mu=tfection_params["nb_mu"]
        ret.transfection_nb_alpha=tfection_params["nb_alpha"]
        
        return ret

#class wet_bounds
#   read_depth