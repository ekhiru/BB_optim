from scipy.spatial import distance
import time
import sys
import os
import numpy as np
import pandas as pd
import mallows_kendall as mk
from pfsp import PFSP
import cego

instance = PFSP.read_instance("pfsp/rec05.txt")
df = cego.runCEGO(instance = instance, m_ini = 10, m = 100, rep=1, best_known_sol=instance.best_sol, worst_known_sol=instance.worst_sol, budgetGA=100)
print(df)
