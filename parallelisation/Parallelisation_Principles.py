#%%
import numpy as np
from multiprocessing import Pool
import multiprocessing


#%%
def f1(x):
    return x*x


#%%
vals = np.arange(0,10)

if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    p = multiprocessing.Pool(5)
    results = p.map(f1, vals)
    print(results)
    p.close()

# if __name__ == '__main__':
#     with Pool(5) as p:
#         print(p.map(f, [1, 2, 3]))
# %%

def f2(x,y):
    print(x*y)
    return x*y


#%%
vals = np.arange(0,8)
vals2 = np.arange(0,8)

if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    p = multiprocessing.Pool(5)
    results = p.starmap(f2, zip(vals,vals2))
    print(results)
    p.close()
# %%
from utils import image_lists

def fL(l1,l2):
    return np.shape(l1)

(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%%
if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    p = multiprocessing.Pool(5)
    results = p.starmap(fL, zip(vals,vals2))
    print(results)
    p.close()

# %% multiprocessing_simpleargs.py
import multiprocessing


def worker(num):
    """thread worker function"""
    print('Worker:', num)


if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = multiprocessing.Process(target=worker, args=(i,))
        jobs.append(p)
        p.start()

# %%
def sum_up_to(number):
    return sum(range(1, number + 1))

if __name__ == '__main__':
    a_pool = multiprocessing.Pool()


    result = a_pool.map(sum_up_to, range(10))

    print(result)