from typing import Dict

import numpy as np
import pandas as pd

import os

from scipy.stats import nbinom

import matplotlib.pyplot as plt
import seaborn as sns




#class wet_bounds
#   read_depth

from pathlib import Path
working_dir = Path(__file__).resolve().parent

#SHENDURE_BOUNDS=Bounds.from_tgz(working_dir/"presets/shendure_bounds.tgz")