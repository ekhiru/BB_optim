import numpy as np
from problem import Problem

class PFSP(Problem):
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

        return cls(P, instance_name = filename, best_sol = best_sol)

    def __init__(self, P, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        # Processing times matrix
        self.P = np.asarray(P)
        # jobs, machines
        self.n, self.m = self.P.shape
        super().__init__(best_sol = best_sol, worst_sol = worst_sol,
                         instance_name = instance_name)

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

    def sum_completion(self, x):
        """Also known as sum of completion times"""
        return self.completion_times(x).sum()
        
    def fitness_nosave(self, x):
        # In case it is not numpy array.
        x = np.asarray(x, dtype=int)
        assert self.check_permutation(x), f"{x}"
        return self.objective(x)

class PFSP_Cmax(PFSP):
    # Class attributes
    problem_name = "PFSP-Cmax"

    def __init__(self, P, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        super().__init__(P, best_sol = best_sol, worst_sol = worst_sol, instance_name = instance_name)
        self.objective = self.makespan

class PFSP_Csum(PFSP):
    # Class attributes
    problem_name = "PFSP-Csum"

    def __init__(self, P, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        super().__init__(P, best_sol = best_sol, worst_sol = worst_sol, instance_name = instance_name)
        self.objective = self.sum_completion

