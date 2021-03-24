import numpy as np
from scipy.spatial.distance import cdist
from problem import Problem

class TSP(Problem):
    # Class attributes
    problem_name = "TSP"

    @classmethod
    def read_instance(cls, filename, opt_filename = None):
        print(f"Reading instance from {filename}")
        n = -1
        with open(filename) as f:
            instance_name = "unknown"
            while True:
                header = f.readline().split(':')
                field = header[0].strip()
                if field == "DIMENSION":
                    n = int(header[1].strip())
                elif field == "EDGE_WEIGHT_TYPE":
                    assert header[1].strip() == "EUC_2D"
                elif field == "NAME":
                    instance_name = header[1].strip()
                elif field == "NODE_COORD_SECTION":
                    break

            print(f"{instance_name} {n} EUC_2D")
            # First column is index
            coords = np.loadtxt(f, max_rows=n, usecols=(1,2), dtype=int)
        dist = cdist(coords,coords).round(0).astype(int)

        best_sol = None
        if opt_filename is not None:
            with open(opt_filename) as f:
                while True:
                    header = f.readline().split(':')
                    field = header[0].strip()
                    if field == "DIMENSION":
                        assert n == int(header[1].strip())
                    elif field == "TOUR_SECTION":
                        break
                best_sol = np.loadtxt(f, max_rows=n, dtype=int)
                best_sol -= 1 # 0-indexed
                print(f"Reading best solution {best_sol} from {opt_filename}")

        return cls(dist, instance_name = filename, best_sol = best_sol)
        

    def __init__(self, dist, best_sol = None, worst_sol = None, instance_name = "(generated)"):
        # Distance matrix
        self._dist = np.asarray(dist)
        assert self._dist.shape[0] == self._dist.shape[1]
        self.n = self._dist.shape[0]
        super().__init__(best_sol = best_sol, worst_sol = worst_sol,
                         instance_name = instance_name)

    def objective(self,x):
        return self._dist[np.roll(x,1), x].sum()
        
    def fitness_nosave(self, x):
        # In case it is not numpy array.
        x = np.asarray(x, dtype=int)
        assert self.check_permutation(x), f"{x}"
        return self.objective(x)


