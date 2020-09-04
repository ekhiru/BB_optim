import numpy as np
from problem import Problem

class QAP(Problem):
    # Class attributes
    problem_name = "QAP"
    @classmethod
    def read_instance(cls, filename, opt_filename = None):
        print(f"Reading instance from {filename}")
        with open(filename) as f:
            n = int(f.readline().strip())
            A = np.loadtxt(f, max_rows=n)
            B = np.loadtxt(f, max_rows=n)
        best_sol = None
        if opt_filename is not None:
            with open(opt_filename) as f:
                for line in f:
                    name, best_sol = line.strip().split(" ")
                    if filename.find(name) >= 0: 
                        best_sol = np.fromstring(best_sol, dtype=int, sep=",")
                        best_sol -= 1 # 0-indexed
                        print(f"Reading best solution {best_sol} from {opt_filename}")
                        break
        
        return QAP(A, B, instance_name = filename, best_sol = best_sol)

    def __init__(self, A, B, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        self.A = np.asarray(A)
        self.B = np.asarray(B)
        assert self.A.shape[0] == self.A.shape[1]
        assert self.B.shape[0] == self.B.shape[1]
        assert self.A.shape[0] == self.B.shape[0]
        self.n = self.A.shape[0]
        self.instance_name = instance_name
        super().__init__(best_sol = best_sol, worst_sol = worst_sol)
        
    def fitness_nosave(self, x):
        return np.sum(self.A * self.B[np.ix_(x, x)])
    
