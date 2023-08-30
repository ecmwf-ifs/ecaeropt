import numpy as np

fname = "data/non_sphere_scaling/Kandler_nonsphere_scaling"
data  = np.loadtxt(fname).T
bins  = data[4]
print(bins)
print(len(bins))
print(len(bins)%3)
