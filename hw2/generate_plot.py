import struct
import numpy as np
from matplotlib import pyplot as plt

# first open the file
f = open("udata.txt", "r")

iterations = f.readline().rstrip('\n')

y = np.loadtxt(f)
f.close()

# we need to calculate the x values
x = np.linspace(0, 1, num=y.size, endpoint=True)

# now plot
fig, ax = plt.subplots()
ax.plot(x,y)

ax.set(title='u(x) at {} points with {} iterations.'.format(x.size, iterations))
ax.grid()

fig.savefig("udata.png")
