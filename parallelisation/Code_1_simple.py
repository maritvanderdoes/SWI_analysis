#%%
import numpy as np
from multiprocessing import Pool
import multiprocessing


#%%
def f1(x):
    print(x)
    return x*x


#%%
vals = np.arange(0,10)

if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    p = multiprocessing.Pool(5)
    results = p.map(f1, vals)
    print(results)
    p.close()