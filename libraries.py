import numpy as np
import scipy
import astropy
import bz2
import pandas as pd
pd.set_option('display.max_rows', 500)
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from astropy.io import fits
from astropy.table import Table
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import units as u
import warnings
warnings.filterwarnings("ignore")
from pylab import figure, cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import kstest
from scipy.stats import binned_statistic
import random
from statsmodels.stats.diagnostic import anderson_statistic as ad_stat
from scipy import stats
from collections import defaultdict
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
