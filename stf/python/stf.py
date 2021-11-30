#!/usr/bin/env python3

import gzip
import numpy as np


class traj:
    def __init__(self,filename,compresslevel=5):
        self.frames_read=0
        self.traj=gzip.open(filename,compresslevel=compresslevel)
        self.header=""
        for i in self.traj:
            if "**" in str(i):
                self.natoms=int(i.split()[-1])
                self.frame=np.zeros((self.natoms,3))
                return
            self.header+=" "+i.rstrip("\n")
    #each call returns the next frame
    #when there are no more frames left, it raises a GeneratorExit
    def next(self,skip=False):
        r=0
        for i in self.traj:
            i=i.decode("utf-8")
            if "*" in i:
                self.frames_read+=1
                return self.frame
            if not skip:
                n=i.rstrip("\n").split()
                self.frame[r][0]=float(n[0])
                self.frame[r][1]=float(n[1])
                self.frame[r][2]=float(n[2])
                r+=1
        raise GeneratorExit #I guess not the best one.
    def get_header(self):
        return self.header
    def get_natoms(self):
        return self.natoms
    def get_frames_read(self):
        return self.frames_read


