

import numpy as np

dataset = np.loadtxt("data/refr_idx/du_ri_table_composite", skiprows=2).T

newdataset_row1 = [None]*26*12
newdataset_row2 = [None]*26*12 
newdataset_row3 = [None]*26*12
newdataset_row4 = [None]*26*12
newdataset_row5 = [None]*26*12

for ii in range(12):
    print(ii, range(ii*26,(ii+1)*26))
    ini = ii*26
    end = (ii+1)*26
    newdataset_row1[ini:end] = dataset[0]
    newdataset_row2[ini:end] = dataset[1]
    newdataset_row3[ini:end] = dataset[2]
    newdataset_row4[ini:end] = np.arange(1,27)
    newdataset_row5[ini:end] = np.ones(26)*(ii+1)

ndata = np.vstack([newdataset_row1, newdataset_row2, newdataset_row3, newdataset_row4, newdataset_row5])
print(ndata)
np.savetxt("testimg",ndata.T, fmt='%4.3e  %4.3e  %4.3e  %d  %d')   
