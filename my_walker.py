import random
from numpy import pi
from my_math import spherical_to_cartesian, distance
import pandas as pd
import sys

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

class walker(object):
    def __init__(self, x, y, z, t):
        """ Sets coordinates to walker.
            Input: coordinates.
        """
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        """ Checks if two walkers are equal.
            True if objects are walker and have the same coordinates"""
        if not isinstance(other, walker):
            return NotImplemented
        return self.x == other.x and self.y == other.y and self.z == other.z

    def start(self, x, y, z, t):
        """ Sets the initial coordinates.
            Input: Carthesian coordinates and MC time."""
        self.startx = x
        self.starty = y
        self.startz = z
        self.starttime = t

    def final(self, x, y, z, t):
        """ Sets the final coordinates.
            Input: Carthesian coordinates and MC time."""
        self.finalx = x
        self.finaly = y
        self.finalz = z
        self.finaltime = t

    def move(self, steps, direction, entry_direction, time, number, delay):
        """ Handles a move of a walker.
            Input: steps - number of steps per move
                   direction - direction lamellae (0-parallel to water entry, 1-perpendicular to water entry)
                   entry_direction - direction of water entry (0-x, 1-y, 2-z)
                   time - MC time before move + 1/number
                   number - number of walker in sphere
                   delay - additional delay, e.g. 10/N corresponds to delay=9
            Output: updated MC time"""
        for i in range(steps):
            (dx, dy, dz) = random.choice([(0, 1, 0), (0, -1, 0), (1, 0, 0), (-1, 0, 0), (0, 0, 1), (0, 0, -1)])
            self.x += dx
            self.y += dy
            self.z += dz
            if direction==0:
                if (dx!=0 and entry_direction == 0) or (dy!=0 and entry_direction == 1) or (dz!=0 and entry_direction == 2):
                    time=time + delay/number
            elif direction==1:
                if ((dy!=0 or dz!=0) and entry_direction == 0) or ((dx!=0 or dz!=0) and entry_direction == 1) or \
                        ((dx!=0 or dy!=0) and entry_direction == 2):
                    time=time + delay/number
            return time

def initial_state(number_of_walkers, radius):
    """ Create starting configuration.
        Input: number_of_walkers, radius
        Output: Updated walkers list"""
    walkers = []
    random.seed
    for i in range(1, number_of_walkers + 1):
        nw = len(walkers)
        while True:
            rn1 = random.randrange(100) / 100
            rn2 = random.randrange(100) / 100
            if (rn1 ** 2 >= rn2):
                r = rn1 * radius
                theta = random.randrange(100) / 100 * pi
                phi = random.randrange(100) / 100 * 2 * pi
                (x, y, z) = spherical_to_cartesian(r, theta, phi)
                new_walker = walker(int(round(x)), int(round(y)), int(round(z)))
                if new_walker not in walkers:
                    walkers.append(new_walker)
                if len(walkers) > nw:
                    break
    return walkers

def random_walk(constants):
    """ Handles a random walk
        Input: constants
        Output: list with degraded walkers"""
    walkers = initial_state(constants.number_of_walkers, constants.radius)
    time = 0.0
    degradation = [(0, time, len(walkers))]
    q = 1 - constants.D_r
    print(q)
    for i in range(1, constants.MC_steps + 1):
        if len(walkers) > 0:
            rn = random.randrange(100) / 100
            time += 1.0 / len(walkers)
            if (rn >= q):
                W_select = walkers[random.randrange(len(walkers))]
                W = walker(W_select.x, W_select.y, W_select.z)
                for walk_length in range(constants.min_walk_length, constants.max_walk_length):
                    W_old = walker(W.x, W.y, W.z)
                    W.move(walk_length)
                    if W in walkers:
                        W = W_old
                        # print("rejected move")
                        # degradation.append((i,time,len(walkers)))
                    else:
                        d = distance(W.x, W.y, W.z)
                        if d < constants.radius:
                            walkers.append(W)
                            walkers.remove(W_old)
                            # print("accepted move", len(walkers), d)
                            # degradation.append((i,time,len(walkers)))
                        else:
                            walkers.remove(W_select)
                            # print("polymer degraded", len(walkers), d)
                            # degradation.append((i,time,len(walkers)))
                            degradation.append((i, time, len(walkers)))
    return degradation

def initial_state_empty_complete(radius):
    """ Creates empty walkers list"""
    walkers = []
    return walkers

def initial_state_empty_oneside(radius):
    """ Creates walkers list with 1 walker at (-radius, 0, 0)"""
    walkers = []
    random.seed
    r = radius
    theta = pi/2
    phi = pi
    (x, y, z) = spherical_to_cartesian(r, theta, phi)
    walkers.append(walker(int(round(x)), int(round(y)), int(round(z)), 0))
    return walkers

def initial_state_empty_allsides(radius):
    """ Creates walkers list with 6 walkers at (+/-radius, 0, 0), (0, +/-radius, 0) and (0, 0, +/-radius)."""
    walkers = []
    random.seed
    r = radius
    for (theta,phi) in [(0,pi),(pi/2,pi),(pi/2,0),(pi/2,pi/2),(pi/2,-pi/2),(pi,pi)]:
        (x, y, z) = spherical_to_cartesian(r, theta, phi)
        walkers.append(walker(int(round(x)), int(round(y)), int(round(z)), 0))
    return walkers

def initial_state_empty_opposite(radius):
    """ Creates walkers list with walker at (-radius, 0, 0)"""
    walkers = []
    random.seed
    r = radius
    for (theta,phi) in [(pi/2,0)]:
        (x, y, z) = spherical_to_cartesian(r, theta, phi)
        walkers.append(walker(int(round(x)), int(round(y)), int(round(z)), 0))
    return walkers

def initial_state_empty_perpendicular(radius):
    """ Creates walkers list with walkers at (0, 0, +/-radius) and (0, +/-radius, 0) """
    walkers = []
    random.seed
    r = radius
    for (theta,phi) in [(0,pi),(pi/2,pi/2),(pi/2,-pi/2),(pi,pi)]:
        (x, y, z) = spherical_to_cartesian(r, theta, phi)
        walkers.append(walker(int(round(x)), int(round(y)), int(round(z)), 0))
    return walkers

def initial_state_load_from_file(file):
    """ Creates walkers list with walkers loaded from file"""
    walkers = []
    walkers_load = pd.read_csv(file, header=None, names=["x","y","z"])
    for row in walkers_load.itertuples(index=True):
        walkers.append(walker(getattr(row, "x"), getattr(row, "y"), getattr(row, "z"), 0))
    return walkers

def initial_state_load_from_file_crystalline(file):
    """ Creates walkers list with walkers loaded from file - lamellae positions"""
    walkers = []
    walkers_load = pd.read_csv(file, header=None, names=["x","y","z"])
    number = 1
    direction = 0
    distance = 1
    length=3
    crystalline = initial_state_crystalline(number, direction, distance, length)
    for row in walkers_load.itertuples(index=True):
        W=walker(getattr(row, "x"), getattr(row, "y"), getattr(row, "z"), 0)
        if W not in crystalline:
            walkers.append(W)
    return walkers

def initial_state_crystalline(number, direction, distance, length):
    """ Cretaes list blocked with all elements of lamellae"""
    blocked = []
    lower=-int(length/2.)
    upper=int(length/2.)
    if length%2==1:
        upper=upper+1
    for i in np.arange(lower,upper):
        for k in np.arange(lower, upper):
            lower2 = -int(distance/2.*(number-1))
            upper2=lower2+(distance*(number-1))
            for j in np.arange(lower2,upper2+1,distance):
                if direction==0: #parallel
                    blocked.append(walker(j,i,k, 0))
                else: #perpendicular
                    blocked.append(walker(i,j,k, 0))
    return blocked

def random_walk_fill(constants):
    """Random walk of water uptake."""
    walkers = initial_state_empty_allsides(constants.radius)
    time = 0.0
    degradation = [(0, time, len(walkers))]
    degradation.append((0., 0., len(walkers)))
    q = 1 - constants.D_r
    i=walkers[0]
    for i in range(1, constants.MC_steps + 1):
        if len(walkers) < constants.number_of_walkers:
            rn = random.randrange(100) / 100
            time += 1.0 / len(walkers)
            if (rn >= q):
                W_select = walkers[random.randrange(len(walkers))]
                W = walker(W_select.x, W_select.y, W_select.z)
                for walk_length in range(constants.min_walk_length, constants.max_walk_length):
                    W_old = walker(W.x, W.y, W.z)
                    W.move(walk_length)
                    if W in walkers:
                        W = W_old
                        # print("rejected move")
                        # degradation.append((i,time,len(walkers)))
                    else:
                        d = distance(W.x, W.y, W.z)
                        if d < constants.radius:
                            walkers.append(W)
                            if not(W_old == walkers[0]):
                                walkers.remove(W_old)
                            else:
                                degradation.append((i, time, len(walkers)))
                                #print(i, time, len(walkers))
                            #print("accepted move", len(walkers), d)
                            #print(i,time,len(walkers))
                        #else:
                        #    walkers.remove(W_select)
                        #    # print("polymer degraded", len(walkers), d)
                        #    # degradation.append((i,time,len(walkers)))
                        #    degradation.append((i, time, len(walkers)))
    return degradation

def random_walk_fill_and_leave(constants):
    """Random walk of water uptake and water exit
       Old implementation"""
    file = "initial_distribution_R10c05.csv"
    #walkers = initial_state_empty_allsides(constants.radius)
    walkers = initial_state_load_from_file(file)
    walkers_leave = initial_state_empty_allsides(constants.radius)
    time = 0.0
    degradation = [(0, time, len(walkers))]
    #degradation.append((0., 0., len(walkers)))
    q = 1 - constants.D_r
    #i=walkers[0]
    xyz = []
    x=[]
    y=[]
    z=[]
    for i in range(1, constants.MC_steps + 1):
        if len(walkers) < constants.number_of_walkers:
            rn = random.randrange(100) / 100
            time += 1.0 / len(walkers)
            if (rn >= q):
                W_select = walkers[random.randrange(len(walkers))]
                W = walker(W_select.x, W_select.y, W_select.z)
                for walk_length in range(constants.min_walk_length, constants.max_walk_length):
                    W_old = walker(W.x, W.y, W.z)
                    W.move(walk_length)
                    if W in walkers:
                        W = W_old
                        # print("rejected move")
                        # degradation.append((i,time,len(walkers)))
                    else:
                        d = distance(W.x, W.y, W.z)
                        if d < constants.radius:
                            walkers.append(W)
                            xyz.append([W.x,W.y,W.z])
                            x = [j[0] for j in xyz]
                            y = [j[1] for j in xyz]
                            z = [j[2] for j in xyz]
                            if not(W_old in walkers_leave):
                                walkers.remove(W_old)
                                xyz.remove([W_old.x,W_old.y,W_old.z])
                            else:
                                degradation.append((i, time, len(walkers)))
                            #print("accepted move", len(walkers), d)
    #pd.DataFrame(xyz).to_csv("initial_distribution_R10c05.csv", header=None, index=None)
    return degradation

def random_walk_diffusion(constants):
    """Random walk of diffusion through volume element."""
    walkers=[]
    walkers_degraded=[]
    file = "initial_distribution_R10c02.csv"
    #walkers = initial_state_empty_oneside(constants.radius)
    walkers_original = initial_state_load_from_file(file)
    direction=0 #0-parallel, 1-perpendicular
    number=1 # number of lamellae in sphere
    dist=2 # distance between two parallel lamellae
    length=11 # =2a+1, length and width of quadratic lamellae, height=1
    delay=9.0 #9.0 # how much slower is motion of MC move in parallel direction of the lamellae in comparison to
              # perpendicular direction
    walkers_crystalline = initial_state_crystalline(number,direction,dist,length)
    for wal in walkers_original:
        if not(wal in walkers_crystalline):
            walkers.append(wal)
    for wal in walkers:
        wal.start(wal.x,wal.y,wal.z, 0)
    walkers_leave = initial_state_empty_oneside(constants.radius)
    if walkers[0].x != 0:
        entry_direction = 0 #x
    elif walkers[0].y != 0:
        entry_direction = 1 #y
    elif walkers[0].z != 0:
        entry_direction = 2 #z
    walkers_opposite = initial_state_empty_opposite(constants.radius)
    walkers_perpendicular = initial_state_empty_perpendicular(constants.radius)
    time = 0.0
    degradation = [(0, time, len(walkers))]
    degradation_parrallel = [(0, time, 0)]
    degradation_perpendicular = [(0, time, 0)]
    parallel_degraded = 0
    perpendicular_degraded = 0
    #degradation.append((0., 0., len(walkers)))
    q = 1 - constants.D_r
    #i=walkers[0]
    xyz = []
    x=[]
    y=[]
    z=[]
    for i in range(1, constants.MC_steps + 1):
        #if len(walkers) < constants.number_of_walkers:
            rn = random.randrange(100) / 100
            time_old = time
            time += 1.0 / len(walkers)
            if (rn >= q):
                W_select = walkers[random.randrange(len(walkers))]
                W = walker(W_select.x, W_select.y, W_select.z, W_select.starttime)
                W.start(W_select.startx, W_select.starty, W_select.startz, W_select.starttime)
                for walk_length in range(constants.min_walk_length, constants.max_walk_length):
                    W_old = walker(W.x, W.y, W.z, W.starttime)
                    time = W.move(walk_length, direction, entry_direction, time, len(walkers), delay)
                    if ((W in walkers) or (W in walkers_crystalline)):
                        W = W_old
                        # print("rejected move")
                        # degradation.append((i,time,len(walkers)))
                    else:
                        d = distance(W.x, W.y, W.z)
                        if d <= constants.radius:
                            #ax = plt.axes(projection='3d')
                            #xyz.append([W.x,W.y,W.z])
                            #x = [j[0] for j in xyz]
                            #y = [j[1] for j in xyz]
                            #z = [j[2] for j in xyz]
                            #ax.scatter3D(x, y, z)
                            #print(x,y,z)
                            if not(W_old in walkers_leave):
                                walkers.append(W)
                                walkers.remove(W_old)
                                #print("normal move",W_old.x,W_old.y,W_old.z)
                                #xyz.remove([W_old.x,W_old.y,W_old.z])
                            else:
                                W_new=walker(W.x,W.y,W.z,time)
                                W_new.start(W_new.x, W_new.y, W_new.z, time)
                                walkers.append(W_new)
                                degradation.append((i, time, len(walkers)))
                                #print("new molecule",W_old.x,W_old.y,W_old.z, len(walkers))
                        else:
                            W.final(W.x, W.y, W.z, time)
                            if (W_old in walkers_opposite):
                                parallel_degraded += 1
                                degradation_parrallel.append((i, time, parallel_degraded))
                                print(i, "parallel", W_old.x, W_old.y, W_old.z,len(walkers))
                                walkers_degraded.append(W)
                                walkers.remove(W_old)
                            elif (W_old in walkers_perpendicular):
                                perpendicular_degraded += 1
                                degradation_perpendicular.append((i, time, perpendicular_degraded))
                                print(i, "perpendicular", W_old.x, W_old.y, W_old.z,len(walkers))
                                walkers_degraded.append(W)
                                walkers.remove(W_old)
    #for wal in walkers:
    #    xyz.append([wal.x, wal.y, wal.z])
    #pd.DataFrame(xyz).to_csv("initial_distribution_R10c02.csv", header=None, index=None)
    for wal in walkers:
        wal.final(wal.x, wal.y, wal.z, time)
        walkers_degraded.append(wal)
    motions = []
    for wal in walkers_degraded:
        x=distance(wal.startx-wal.finalx, wal.starty-wal.finaly, wal.startz-wal.finalz)
        t=wal.finaltime-wal.starttime
        motions.append((x, t))
    print("pa:",parallel_degraded)
    print("pe:",perpendicular_degraded)
    print(len(walkers))
    return [degradation, degradation_parrallel, degradation_perpendicular, motions]