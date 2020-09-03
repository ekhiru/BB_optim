import numpy as np

from rpy2.robjects import FloatVector
import rpy2.rinterface as ri

class QAP:
    # Class attributes
    problem_name = "QAP"
    @classmethod
    def read_instance(cls, filename, opt_filename):
        print(f"Reading instance from {filename}")
        with open(filename) as f:
            n = int(f.readline().strip())
            A = np.loadtxt(f, max_rows=n)
            B = np.loadtxt(f, max_rows=n)
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
        self.evaluations = []
        self.solutions = []
        self.best_sol = best_sol
        self.worst_sol = worst_sol

    
    def fitness_nosave(self, x):
        return np.sum(self.A * self.B[np.ix_(x, x)])
       
    # Minimized
    def fitness(self, x):
        f = self.fitness_nosave(x)
        self.solutions.append(x)
        self.evaluations.append(f)
        return f

    # Returns a closure function that can be called from R.
    # WARNING: this function minimizes for CEGO
    # FIXME: Can we make this a function shared by all problems instead of copy-pasting?
    def make_r_fitness(self):
        @ri.rternalize
        def r_fitness(x):
            xpy = np.asarray(x) - 1 # R vectors are 1-indexed
            y = self.fitness(xpy)
            return FloatVector(np.asarray(y))
        return r_fitness



  
