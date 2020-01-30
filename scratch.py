import numpy as np
import numpy.random as rng
import matplotlib.pylab as plt

def parse(rdr):
    data = []
    for ln in rdr:
        if "ERROR" in ln:
            continue
        data.append([float(w) for w in ln.split()])

    data = np.array(data)
    x    = data[:,0] - np.sqrt(data[:,1]*data[:,2])
    y    = data[:,3]

    return x, y

def simulate(N=1000):
    f = lambda x: .25*(1-(1/(1+np.exp(-(x-3))))) + .25
    x = np.linspace(0, 12, N)
    F = f(x)

    return x, np.abs(F + .05*rng.randn(*x.shape))

