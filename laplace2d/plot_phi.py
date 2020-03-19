import numpy as np
import matplotlib.pyplot as plt

phi = np.loadtxt("phi.csv", delimiter=",").T

plt.imshow(phi, interpolation='bicubic', origin='lower')
plt.show()
