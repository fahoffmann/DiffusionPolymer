import pandas as pd
import time
from concurrent.futures import ProcessPoolExecutor
from my_walker import random_walk, random_walk_fill, random_walk_fill_and_leave, random_walk_diffusion
import sys


def monte_carlo(constants):
    """ Monte Carlo program
        Performs simultanously as many iterations as available CPUs and returns concatenated results."""
    time_all = pd.DataFrame()
    length_all = pd.DataFrame()
    time_pa_all = pd.DataFrame()
    length_pa_all = pd.DataFrame()
    time_pe_all = pd.DataFrame()
    length_pe_all = pd.DataFrame()
    distance_walker_all = pd.DataFrame()
    time_walker_all = pd.DataFrame()
    procs = constants.number_of_iterations  # Number of processes to create
    n_cpus = 8 # multiprocessing.cpu_count()
    pool = ProcessPoolExecutor(max_workers=n_cpus)
    res = []
    for i in range(0, procs):
        print(constants.radius, i + 1)
        a = pool.submit(iteration, constants=constants)
        res.append(a)
        time.sleep(1)
    #sys.stdout.flush()
    for i in range(0, procs):
        if constants.diffusion:
            time_all=pd.concat([time_all,list(res[i].result())[0]],axis=1)
            length_all = pd.concat([length_all, list(res[i].result())[3]], axis=1)
            time_length_all = pd.concat([time_all, length_all], axis=1)
            time_pa_all = pd.concat([time_pa_all, list(res[i].result())[1]], axis=1)
            length_pa_all = pd.concat([length_pa_all, list(res[i].result())[4]], axis=1)
            time_length_pa_all = pd.concat([time_pa_all, length_pa_all], axis=1)
            time_pe_all = pd.concat([time_pe_all, list(res[i].result())[2]], axis=1)
            length_pe_all = pd.concat([length_pe_all, list(res[i].result())[5]], axis=1)
            time_length_pe_all = pd.concat([time_pe_all, length_pe_all], axis=1)
            distance_walker_all = pd.concat([distance_walker_all, list(res[i].result())[6]], axis=1)
            time_walker_all = pd.concat([time_walker_all, list(res[i].result())[7]], axis=1)
            distance_time_walker_all = pd.concat([distance_walker_all, time_walker_all], axis=1)
        else:
            time_all = pd.concat([time_all, res[i].result()], axis=1)
    if constants.diffusion:
        return time_length_all, time_length_pa_all, time_length_pe_all, distance_time_walker_all
    else:
        return time_all

def iteration(constants):
    """ MC iteration.
        Performs random walk and converts result into DataFrame."""
    if constants.diffusion:
        result = random_walk_diffusion(constants)
        [degradation, degradation_parallel, degradation_perpendicular, motions] = [zip(*result[i]) for i in range(
            4)]
        mc, time, length = degradation
        mc_pa, time_pa, length_pa = degradation_parallel
        mc_pe, time_pe, length_pe = degradation_perpendicular
        distance_walker, time_walker = motions
        return pd.DataFrame(list(time)), pd.DataFrame(list(time_pa)), pd.DataFrame(list(time_pe)), \
               pd.DataFrame(list(length)), pd.DataFrame(list(length_pa)), pd.DataFrame(list(length_pe)), \
               pd.DataFrame(list(distance_walker)), pd.DataFrame(list(time_walker))
    else:
        degradation_formated = zip(*random_walk_diffusion(constants))
        mc, time, length = degradation_formated
        return pd.DataFrame(list(time))