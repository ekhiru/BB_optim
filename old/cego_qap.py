import time
import sys
import os
import pandas as pd
import mallows_kendall as mk
from qap import QAP
import cego

instance = QAP.read_instance("qap/nug12.dat", "qap/optimal.txt")
df = cego.runCEGO(instance = instance, m_ini = 10, m = 100, rep=1, best_known_sol=instance.best_sol, worst_known_sol=instance.worst_sol, budgetGA=100)
print(df)
