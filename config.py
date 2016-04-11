MAX_EXP = 300
GRID_DEPTH = 3
INF = float('inf')
VERBOSE = True
USE_BIGFLOAT = False
STEP = 1.1
OPTIMIZATION_METHOD = 'L-BFGS-B'
# OPTIMIZATION_METHOD  = 'TNC'
INITIAL_GRID_COUNT = 16
INITIAL_GRID_STEP = 3
ERR_SCALE = 100
try:
    from multiprocessing import cpu_count

    DEFAULT_THREAD_COUNT = cpu_count()
except NotImplementedError:
    DEFAULT_THREAD_COUNT = 2
