import timeit

def tic():
    start = timeit.default_timer()
    return start

def toc(start):
    stop = timeit.default_timer()
    print('Time elapsed of ' +"{:.1f}".format((stop - start))+' seconds')
    return stop