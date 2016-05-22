from __future__ import absolute_import, division, print_function
import os, math, itertools
from collections import defaultdict
import numpy as np

from utils import readvicsekdata, writevicsekdata, readlendata, writelendata, IsEqual



infile = "data_exp.txt"
infile2 = "data_len.txt"
outfile = "data_exp_noduplicates.txt"
outfile2 = "data_len_noduplicates.txt"
start = 0
stop = 2000
step = 1


indir = ""
outdir = ""



with open(os.path.join(indir, infile), "r") as fr, open(os.path.join(indir, infile2), "r") as fr2, open(os.path.join(outdir, outfile), "w") as fw, open(os.path.join(outdir, outfile2), "w") as fw2:
    fw.write("W: 512 H: 512\n") #WRITE MANUALLY
    fw2.write("W: 512 H: 512\n") #WRITE MANUALLY
    for (indices, positions, angles), (_, lengths) in itertools.izip(readvicsekdata(fr, start, stop, step), readlendata(fr2, start, stop, step)):
        unique_positions = []
        unique_angles = []
        unique_leng = []
        for pos, angle, leng in zip(positions, angles, lengths):
            duplicate = False
            for upos in unique_positions:
                if IsEqual(pos, upos, 0.001):
                    duplicate = True
                    print("line", 2*indices+2, "x,y", pos)
            if not duplicate:
                unique_positions.append(pos)
                unique_angles.append(angle)
                unique_leng.append(leng)
        writevicsekdata(fw, unique_positions, unique_angles)
        writelendata(fw2, unique_leng)

