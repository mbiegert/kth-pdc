import struct
import numpy as np
from matplotlib import pyplot as plt

# first open the file
f = open("udata.txt", "r")

iterations = f.readline().rstrip('\n')

data = np.loadtxt(f)
f.close()

# split the data into values and the error
y = data[:,0]
err = data[:,1]

# we need to calculate the x values
x = np.linspace(0, 1, num=y.size, endpoint=True)

# now plot
fig = plt.figure()
ax = fig.add_subplot(211)
ax.plot(x,y)
ax.set(title='u(x) at {} points with {} iterations.'.format(x.size, iterations))
ax.grid()

ax = fig.add_subplot(212)
ax.plot(x, err)
ax.set(title='error')


fig.savefig("udata.png")
