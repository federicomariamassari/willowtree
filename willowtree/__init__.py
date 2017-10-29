'''
__init__.py is automatically run with the Python command 'import willowtree'.

The notation 'from .<module> import <attribute>' reads: search for <module> in
the current directory (.) and import <attribute>.

<attribute> will now be available to call as willowtree.<attribute>, instead
of as willowtree.<module>.<attribute>.

For example, on Terminal (macOS) or Command Prompt (Windows), the following
commands will both run function maketree:

$ python3
>>> import willowtree
>>> willowtree.maketree()

or, using wt as an alias:

>>> import willowtree as wt
>>> wt.maketree()

which is obviously less burdensome than willowtree.maketree.maketree().
'''
import time
import numpy as np
from scipy import stats, optimize
import matplotlib.pyplot as plt
import seaborn as sns

from .__version__ import __version__
from .maketree import maketree
from .sampling import sampling
from .lp import lp
from .graph import graph
