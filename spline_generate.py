import os, itertools, math, logging
from collections import defaultdict, OrderedDict

from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plot
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from utils import *

indir = "."
outdir = "."
    
infile = "data_exp_tracked.txt"
infilevel = "data_vel_tracked.txt"
infileids = "data_exp_tracked_ids.txt"
outfile = "data_exp_tracked_splines.txt"
outfilevel = "data_vel_tracked_splines.txt"

L = 512
stop = 4000
smoothing = 2000.0

"""trajectories = defaultdict(list)

with open(os.path.join(indir, infile), "r") as fr:
    for time, (ids, positions) in enumerate(tracking_ids(fr, stop=50)):
        for id, pos in zip(ids, positions):
            trajectories[id].append([time, pos[0], pos[1]])"""

trajectories = defaultdict(list)
trajectories_angles = defaultdict(list)
trajectories_velocities = defaultdict(list)
ids_list = []

with open(os.path.join(indir, infile), "r") as fr1, open(os.path.join(indir, infileids), "r") as fr2, open(os.path.join(indir, infilevel), "r") as frvel:
    for time, ((_, positions, angles), ids, (_, vels)) in enumerate(itertools.izip(readvicsekdata(fr1, stop=stop, step=2), readiddata(fr2, stop=stop), readlendata(frvel, stop=stop))):
        assert len(ids) == len(set(ids)), ids
        assert len(positions) == len(ids), ids
        for id, pos, angle, vel in zip(ids, positions, angles, vels):
            trajectories[id].append([time, pos[0], pos[1]])
            trajectories_angles[id].append(angle)
            trajectories_velocities[id].append([vel[0], vel[1]])
            #if id == 65:
            #    print time
        ids_list.append(ids)

#trajectories = {id:poss for id, poss in trajectories.iteritems() if len(poss) >= 500}

def fix_angles(angles, x_i_der, y_i_der):
    anglevectors = [Vector2D(math.cos(angle), math.sin(angle)) for angle in angles]
    last_trajv = anglevectors[0]
    
    for index, (ang, angv, dx, dy) in enumerate(itertools.izip(angles, anglevectors, x_i_der, y_i_der)):
        #translated from "ang_distance > PI/2" from tracking.cpp
        #angv*trajv follows a Plot[{1,0}.{Cos[a],Sin[a]},{a,-Pi,Pi}] function
        try:
            trajv = Vector2D(dx, dy).Normalized()
        except RuntimeError:
            logging.warning("Vector (from spline fit) has length 0, falling back to detection: {} {}".format(last_trajv.x, last_trajv.y))
            trajv = last_trajv
        if angv*trajv < 0:
            if (ang > 0):
                angles[index] -= math.pi;
            else:
                angles[index] += math.pi;
            last_trajv = trajv
    return angles


for (id, positions), (id2, angles) in itertools.izip(trajectories.iteritems(), trajectories_angles.iteritems()):
    #assert id==id2 and id2==id3
    length = len(positions)
    t_start = positions[0][0]
    if length > 3:
        #print(positions)
        """knots = np.arange(t_start, t_start+length, 1)
        t,x,y = zip(*positions)
        tck= interpolate.splrep(t, x, s=smoothing, k=3) #, u=np.linspace(0, 1, length),  u=range(t_start, t_start+length)
        x_i = interpolate.splev(t, tck, ext=2) #(0, 1, length)
        #if id == 65:
        #    print t
        #    for i in zip(t, map(float, x_i)):
        #        print i
        trajectories[id] = zip(t, x_i, x_i)"""
        tck, u = interpolate.splprep(zip(*positions), u=range(t_start, t_start+length), s=smoothing, k=3) #, u=np.linspace(0, 1, length),  u=range(t_start, t_start+length)
        t_i, x_i, y_i = interpolate.splev(u, tck, ext=2) #(0, 1, length)
        t_i_der, x_i_der, y_i_der = interpolate.splev(u, tck, der=1, ext=2) #(0, 1, length)
        #if id == 65:
        #    for i in zip(map(float, t_i), map(float, x_i),  map(float, y_i)):
        #        print i
        assert len(t_i) == len(set(np.rint(t_i).astype(int))), u
        trajectories[id] = zip(np.rint(t_i).astype(int), x_i, y_i)
        trajectories_angles[id] = fix_angles(angles, x_i_der, y_i_der)
        trajectories_velocities[id] = zip(x_i_der, y_i_der)
        #print np.rint(t_i).astype(int), length, t_start

        """if length > 400:
            fig = plot.figure()
            ax = fig.gca(projection='3d')
            ax.plot(t_i, x_i, y_i, label='parametric curve')
            t, x, y = zip(*positions)
            ax.scatter(t, x, y)
            #ax.quiver(t_i, x_i, y_i, t_i_der, x_i_der, y_i_der)
            ax.legend()
            plot.show()"""

        #break #TODO DELETE
    elif length == 3:
        tck, u = interpolate.splprep(zip(*positions), u=range(t_start, t_start+length), s=smoothing, k=2)
        t_i, x_i, y_i = interpolate.splev(u, tck, ext=2)
        t_i_der, x_i_der, y_i_der = interpolate.splev(u, tck, der=1, ext=2)
        assert len(t_i) == len(set(np.rint(t_i).astype(int))), u
        trajectories[id] = zip(np.rint(t_i).astype(int), x_i, y_i)
        trajectories_angles[id] = fix_angles(angles, x_i_der, y_i_der)
        trajectories_velocities[id] = zip(x_i_der, y_i_der)
    else:
        pass


times = []
timesang = []
timesvel = []
for ids in ids_list:
    times.append(OrderedDict(zip(ids, itertools.repeat(None))))
    timesang.append(OrderedDict(zip(ids, itertools.repeat(None))))
    timesvel.append(OrderedDict(zip(ids, itertools.repeat(None))))

for (id, poslist), (id2, anglist), (id3, vellist) in itertools.izip(trajectories.iteritems(), trajectories_angles.iteritems(), trajectories_velocities.iteritems()):
    #print id
    for (t, x, y), angle, (dx, dy) in itertools.izip(poslist, anglist, vellist):
        #print(id, t)
        assert times[t][id] is None, "id={}, t={}, x={}, y={}".format(id, t, x, y)
        times[t][id] = (x, y)
        timesang[t][id] = angle
        timesvel[t][id] = (dx, dy)

with open(os.path.join(indir, infileids), "r") as fr2:
    for t, (timelist, ids) in enumerate(itertools.izip(times, readiddata(fr2, stop=stop))):
        assert len(timelist) == len(ids), t

with open(os.path.join(outdir, outfile), "w") as fw:
    fw.write("W: {} H: {} smoothing: {}\n".format(L, L, smoothing))
    for poslist, anglist in itertools.izip(times, timesang):
        fw.write(" ".join("{} {}".format(*coord) for coord in poslist.itervalues()))
        fw.write("\n")
        fw.write(" ".join(str(angle) for angle in anglist.itervalues()))
        fw.write("\n")

with open(os.path.join(outdir, outfilevel), "w") as fw:
    fw.write("W: {} H: {} smoothing: {}\n".format(L, L, smoothing))
    for vellist in timesvel:
        fw.write(" ".join("{} {}".format(*coord) for coord in vellist.itervalues()))
        fw.write("\n")

