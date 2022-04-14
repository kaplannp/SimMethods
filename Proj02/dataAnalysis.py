import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
N_RUNS=1000

#THIS IS HOW TSV's ARE DONE C++!!!!
test = pd.read_csv("Out/out1.tsv", sep='\t', names=['R1','R2','R3']) 

X = -1*np.ones((N_RUNS, test.shape[0], test.shape[1])).astype(float)
for i in range(N_RUNS):
    filepath = "Out/out{}.tsv".format(i)
    X[i] = np.loadtxt(filepath,delimiter='\t')
means = X.mean(axis=0)
stds = X.std(axis=0)

#plot the nice curves
def plotCurve(ind, c):
    plt.plot(means[:,ind], color=c, label="species{}".format(ind))
    plt.plot(means[:,ind]+stds[:, ind], color=c, linestyle='dashed')
    plt.plot(means[:,ind]-stds[:, ind], color=c, linestyle='dashed')
plotCurve(0, 'c')
plotCurve(1, 'g')
plotCurve(2, 'b')
plt.legend()
plt.title("Decomposition of Species")
plt.ylabel("Population of Species")
plt.xlabel("Time")
plt.show()
plt.savefig("population_graph.png")
