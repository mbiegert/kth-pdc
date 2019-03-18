import sys
import struct
import numpy as np
from matplotlib import pyplot as plt

# first open the file
f = open(sys.argv[1], "r")

data = np.loadtxt(f)
f.close()

# check order of elements
prior_number = 0.0

for cur_number in data:
    if cur_number < prior_number:
        print("wrong order!")
        exit()
    prior_number = cur_number

print("Correct order")

