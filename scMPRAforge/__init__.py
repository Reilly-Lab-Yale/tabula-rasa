
#just expose everything to the API.
#minimal overhead for faster dev speed
#change later to expose individual functions explicitly

#need to also import core and util dependencies here
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging

# Create a logger
logger = logging.getLogger("scMPRAforge")
logger.setLevel(logging.INFO)  # Default to showing INFO and above

# Create a console handler that prints messages to stdout
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Define a simple formatter
formatter = logging.Formatter("%(name)s: %(levelname)s: %(message)s")
console_handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(console_handler)


from .core import *  # Imports everything public from core.py
from .utils import *  # Imports everything public from utils.py

__all__ = [name for name in dir() if not name.startswith("_")]  # Expose all non-private functions


__version__ = "-1"

#Delete symbols used internally
del pd, sns, plt, logging
