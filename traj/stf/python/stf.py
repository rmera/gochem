#!/usr/bin/env python3

import io
import numpy as np
import zstandard as zstd #you need to have this one installed. It's not part of the std library.

class rtraj:
    def __init__(self,filename):
        self.prec=pow(10,precision)
        self.frames_read=0
        self.file=open(filename,'rb')
        self.dctx = zstd.ZstdDecompressor()
        self.stream = self.dctx.stream_reader(self.file)
        self.traj = io.TextIOWrapper(self.stream, encoding='utf-8')
        self.header=""
        self.headerdict={}
        self.readable=True
        for i in self.traj:
#            i=i.decode("utf-8")
            if "**" in i:
                self.natoms=int(i.split()[-1])
                self.frame=np.zeros((self.natoms,3))
                return
            self.header+=" "+i.rstrip("\n")
        for i in self.header.split():
            k,v=i.split("=")
            self.headerdict[k]=v
        try:
            p=self.headerdic["prec"]
        except KeyError:
            p="2"
        self.prec=pow(10,int(p))
    #each call returns the next frame
    #when there are no more frames left, it raises an EOFError
    #This reads every frame to the same array, to save memory. This is normally 
    def next(self,skip=False, box=[]):
        if not self.readable:
            raise EOFError #I guess not the best one, but I hate exceptions.
        r=0
        for i in self.traj:
  #          i=i.decode("utf-8")
            if "*" in i:
                if len(box)>=9:
                    line=i.rstrip().split()
                    if len(line)>=10:
                        for i,v in enumerate(line[1:]):
                            box[i]=float(v)
                self.frames_read+=1
                return self.frame
            if not skip:
                n=i.rstrip("\n").split()
                self.frame[r][0]=float(n[0])/self.prec
                self.frame[r][1]=float(n[1])/self.prec
                self.frame[r][2]=float(n[2])/self.prec
                r+=1
        raise EOFError #I guess not the best one.
    def next_list(self,skip=False):
        f=self.next(skip)
        return f.tolist()
    def get_header(self):
        return self.header
    def get_natoms(self):
        return self.natoms
    def get_frames_read(self):
        return self.frames_read
    def close(self):
        if self.readable:
            self.traj.close()
            self.readable=False


#wtraj is a write-mode trajectory. It requires the name of the file
#And natoms, the number of atoms per frame
#d is a dictionary of strings.
class wtraj:
    def __init__(self, filename,natoms,compressionlevel=9,d=None,precision=2):
        self.prec=pow(10,precision)
        self.natoms=natoms
        self.file=open(filename,'wb')
        self.cctx=zstd.ZstdCompressor(level=compressionlevel)
        self.traj=self.cctx.stream_writer(self.file)
        if precision!=2:
            if not d:
                d={"prec":str(precision)}
            else:
                d["prec"]=str(precision)
        if d:
            for k,v in d:
                self.traj.write(b"%s=%s\n"%(str(k),str(v))) #I'll attempt converting things to string, but it's your responsibility to ensure that is possible.
            self.traj.write(b"** %d\n"%natoms)
    #Writes the next frame from data, which needs to be an Nx3 list of floats or numpy array.
    def wnext(data, box=[]):
        if len(data)<self.natoms or len(data[0])<3:
            raise ValueError
        for i in range(self.natoms):
            d=data[i]
            self.traj.write(b"%d %d %d\n"%(round(d[0]*self.prec),round(d[1]*self.prec),round(d[2]*self.prec)))
        if len(box)>=9:
            b=box
            self.traj.write(b"%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n"%(b[0],b[1],
                b[2],b[3],b[4],b[5],b[6],b[7],b[8]))
        self.traj.write(b"*\n")
    def close(self):
        if self.readable:
            self.traj.close()
            selr.readable=False


