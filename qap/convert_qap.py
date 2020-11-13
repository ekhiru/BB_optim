import numpy as np

def convert_qaplib(filename):
    with open(filename) as f:
        n = int(f.readline().strip())
        X = np.loadtxt(f, dtype=int)
        X = X.reshape((2*n,n))
    fnew = open(filename + "_new", "w")
    fnew.write(str(n) + "\n")
    np.savetxt(fnew, X, fmt="%d")
    fnew.close()
    print("Write " + filename + "_new")    
