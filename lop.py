from problem import Problem
import mallows_kendall as mk
import numpy as np
import re

def synthetic_LOP(n, m, phi):
    instance = np.zeros((n,n))
    s = np.asarray(mk.samplingMM(m, n, phi=phi, k=None))
    central = np.random.permutation(n)
    s = s[:,central] #compose
    for i in range(n):
        for j in range(i+1,n):
            instance[i,j] = (s[:,i] < s[:,j]).sum() / m
            instance[j,i] = 1 - instance[i,j]
    return instance, central

# Linear Ordering Problem
class LOP(Problem):
    # Class attributes
    problem_name = "LOP"

    @classmethod
    def generate_synthetic(cls, n, m, phi, seed = None):
        if seed is None:
            seed = np.random.randint(1, 123456789, 1)[0]
        np.random.seed(seed)
        instance, best_sol = synthetic_LOP(n, m, phi)
        # best_sol = list(range(n))
        worst_sol = best_sol[::-1]
        return cls(n, instance, best_sol = best_sol, worst_sol = worst_sol,
                   instance_name = f"LOP-synthetic,seed={seed},n={n},m={m},phi={phi}")

    @classmethod
    def read_instance(cls, filename, opt_filename = "lop/best_knowns.csv"):
        if "synthetic" in filename:
            seed, n, m, phi = re.search("seed=([0-9]+),n=([0-9]+),m=([0-9]+),phi=([^ ]+)", filename).group(1,2,3,4)
            print(f"Generating synthetic LOP instance with seed={seed} n={n} m={m} phi={phi}")
            return cls.generate_synthetic(int(n), int(m), float(phi), int(seed))
        else:
            print(f"Reading instance from {filename}")
            with open(filename) as f:
                n = int(f.readline().strip())
                instance = np.loadtxt(f, max_rows=n)
            best_fitness = None
      
            if opt_filename is not None:
                with open(opt_filename) as f:
                    for line in f:
                        name, value = line.strip().split("\t")
                        if filename.find(name) >= 0:
                            # The instances in http://grafo.etsii.urjc.es/optsicom/lolib/#instances and best sols try to maximize the superdiagonal.
                            # We minimize the subdiagonal, which is equivalent, but the best sol needs to be updated
                            best_fitness = instance.sum() - int(value)
                            print(f"Reading best-known fitness {best_fitness} from {opt_filename}")
                            break
                    if best_fitness == None:
                        print(f"Instance {filename} not found in {opt_filename}")
            return LOP(n, instance, best_fitness = best_fitness, instance_name = filename)

    # Methods
    def __init__(self, n, instance, best_sol = None, worst_sol = None, instance_name = None, best_fitness = None):
        self.n = n
        # FIXME: Can we find a better name than instance?
        self.instance = instance
        super().__init__(best_sol = best_sol, worst_sol = worst_sol, instance_name = instance_name, best_fitness = best_fitness)

        print("identity, reverse and best-known fitnesses",self.fitness_nosave(np.arange(n)),self.fitness_nosave(np.arange(n)[::-1]), best_fitness)

    def fitness_nosave(self, x):
        # In case it is not numpy array.
        x = np.asarray(x, dtype=int)
        xinverse = np.argsort(x)#### OJOOOOO se evalua la inversa
        # Sum of the upper triangle. We have to minimize this.
        f = np.tril(self.instance[np.ix_(xinverse, xinverse)]).sum()
        return f
