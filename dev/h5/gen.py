import numpy as np
import pandas as pd

filename = 'test.h5'

df = pd.DataFrame(np.arange(10).reshape((5,2)), columns=['A', 'B'])
print(df)
#    A  B
# 0  0  1
# 1  2  3
# 2  4  5
# 3  6  7
# 4  8  9

# Save to HDF5
df.to_hdf(filename, 'data', mode='w', format='table')
del df    # allow df to be garbage collected

# Append more data
df2 = pd.DataFrame(np.arange(10).reshape((5,2))*10, columns=['A', 'B'])
df2.to_hdf(filename, 'data', append=True)

print(pd.read_hdf(filename, 'data'))
