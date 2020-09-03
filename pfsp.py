import numpy as np

from rpy2.robjects import FloatVector
import rpy2.rinterface as ri

class PFSP:
    # Class attributes
    problem_name = "PFSP"
    @classmethod
    def read_instance(cls, filename, opt_filename = None):
        print(f"Reading instance from {filename}")
        with open(filename) as f:
            header = f.readline()
            print(f"header: {header}")
            # jobs, machines
            n, m = f.readline().strip().split()
            n, m = int(n), int(m)
            P = np.loadtxt(f, max_rows=n)
            # Processing times: each even column is useless.
            P = P[:,1::2]
            assert P.shape[0] == n
            assert P.shape[1] == m

        best_sol = None
        if opt_filename is not None:
            with open(opt_filename) as f:
                for line in f:
                    name, best_sol = line.strip().split()
                    if filename.find(name) >= 0: 
                        best_sol = np.fromstring(best_sol, dtype=int, sep=",")
                        best_sol -= 1 # 0-indexed
                        print(f"Reading best solution {best_sol} from {opt_filename}")
                        break
                    
        return PFSP(P, instance_name = filename, best_sol = best_sol)
        
    def __init__(self, P, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        # Processing times matrix
        self.P = np.asarray(P)
        # jobs, machines
        self.n, self.m = self.P.shape
        
        self.instance_name = instance_name
        self.evaluations = []
        self.solutions = []
        self.best_sol = best_sol
        self.worst_sol = worst_sol

    def completion_times(self, x):
        C = self.P[x, :]
        C[:, 0] = C[:, 0].cumsum()
        C[0, :] = C[0, :].cumsum()
        for i in range(1, self.n):
            for j in range(1, self.m):
                C[i, j] += max(C[i - 1, j], C[i, j - 1])
        return C

    def makespan(self, x):
        C = self.completion_times(x)
        return C[self.n - 1, self.m - 1]
    
    def fitness_nosave(self, x):
        return self.makespan(x)
       
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
