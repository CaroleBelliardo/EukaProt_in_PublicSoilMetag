
#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from time import time
import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())


# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[200000, 5])
data = arr.tolist()
data[:5]

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(results[:10])
#> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]




#nn='test1'
#str1 = "echo "+nn+" > temp.t" # duplique kraken file with header
#content = os.popen(str1).read()



#f=open('workfile', 'rw')
#f.write(b'0123456789abcdef')
