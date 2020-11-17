from problem import Problem
import mallows_kendall as mk
import numpy as np
import re
import os

def generate_list_of_instances(opt_filename):
    with open(opt_filename) as f:
        for line in f:
            name, value = line.strip().split("\t")
            instance_name = find_in_lop_folder(name, "path")
            if instance_name:
                print(instance_name)


def find_in_lop_folder(instance_name, ret_value='instance'):
    for fol in os.listdir("lop"):
        subdir = "lop/" + fol
        if os.path.isdir(subdir):
            for file in os.listdir(subdir):
                if instance_name == file:
                    path = subdir + "/" + file
                    with open(path) as f:
                        n = int(f.readline().strip())
                        instance = np.loadtxt(f, max_rows=n)
                    if ret_value == 'instance' : return instance
                    if ret_value == 'path': return path
    return None

def read_best_known(opt_filename, instance_name):
    with open(opt_filename) as f:
        for line in f:
            name, value = line.strip().split("\t")
            if instance_name.find(name) >= 0:
                # The instances in http://grafo.etsii.urjc.es/optsicom/lolib/#instances and best sols try to maximize the superdiagonal.
                # We minimize the subdiagonal, which is equivalent, but the best sol needs to be updated
                return int(value)

    print(f"Instance {instance_name} not found in {opt_filename}")
    return None

def synthetic_LOP(n, m, phi):
    instance = np.zeros((n,n))
    # FIXME: It should already return an array.
    s = np.asarray(mk.samplingMM(m, n, phi=phi, k=None))
    central = np.random.permutation(n)
    s = s[:, central] #compose
    for i in range(n):
        for j in range(i+1, n):
            instance[i, j] = (s[:, i] < s[:, j]).sum() / m
            instance[j, i] = 1 - instance[i, j]
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
        worst_sol = np.argsort(np.argsort(best_sol)[::-1])
        #OJO, si la fitness NO evalua la inversa, las bests y worst son estas de abajo
        worst_sol = np.argsort(worst_sol)
        best_sol = np.argsort(best_sol)
        return cls(n, instance, best_sol = best_sol, worst_sol = worst_sol,
                   instance_name = f"LOP-synthetic,seed={seed},n={n},m={m},phi={phi}")

    @classmethod
    def read_instance(cls, filename, opt_filename = "./lop/best_knowns.csv"):
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
                best_fitness = read_best_known(opt_filename, filename)
                if best_fitness != None:
                    best_fitness = instance.sum() - best_fitness
                    print(f"Reading best-known fitness {best_fitness} from {opt_filename}")

            return LOP(n, instance, best_fitness = best_fitness, instance_name = filename)

    # Methods
    def __init__(self, n, instance, best_sol = None, worst_sol = None, instance_name = None, best_fitness = None):
        self.n = n
        # FIXME: Can we find a better name than instance?
        self.instance = instance
        super().__init__(best_sol = best_sol, worst_sol = worst_sol, instance_name = instance_name, best_fitness = best_fitness)
        identity = np.arange(n)
        print("identity, reverse, best-known and worst-known fitnesses",
              self.fitness_nosave(identity),
              self.fitness_nosave(identity[::-1]),
              self.best_fitness,
              self.worst_fitness)

    def fitness_nosave(self, x):
        # In case it is not numpy array.
        x = np.asarray(x, dtype=int)
        # Sum of the lower triangle. We have to minimize this.
        f = np.tril(self.instance[np.ix_(x, x)]).sum()
        return f
