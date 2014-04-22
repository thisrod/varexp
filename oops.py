import numpy as np
import matplotlib.pyplot as plt

t = np.arange(19)
plt.ylabel("Military spending, trillions")
plt.plot(1971+t, 5*np.exp(t/10.), "bp", label="USA")
plt.plot(1971+t, 19-t, "r*", label="USSR")
plt.legend( ("USSR", "USA"), 'upper left')
plt.savefig("oops.png")