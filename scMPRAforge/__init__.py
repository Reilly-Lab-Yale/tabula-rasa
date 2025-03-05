
#just expose everything to the API.
#minimal overhead for faster dev speed
#change later to expose individual functions explicitly

#need to also import core and util dependencies here
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


from .core import *  # Imports everything public from core.py
from .utils import *  # Imports everything public from utils.py

__all__ = [name for name in dir() if not name.startswith("_")]  # Expose all non-private functions


__version__ = "-1"

#Delete symbols used internally
del pd, sns, plt
