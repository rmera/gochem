#!/usr/bin/env python3


"""
Copyright (c) 2022 Raul Mera A.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


import io
import numpy as np
import zstandard as zstd #you need to have this one installed. 
#It's not part of the std library. Nevertheless, it is
#Available in the Fedora repos, and, I assume, also in those
#from Ubuntu, Arch, etc., as well as via PIP.


#This reader is quite slow, which is hardly surprising from a
#pure-Python implementation (othern than z-standard and a couple 
#of Numpy funcions). It is more of a proof of concept, although
#it can be used if you don't mind waiting a bit for the trajectory
#to load.
class rtraj:
    def __init__(self,filename):
        self.prec=100
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
    #This reads every frame to the same array
    #if box is given, it should be an array of 9 or more floats, or compabile type, (only the first 9 will be considered)
    #with the 3 vectors of the box, in A.
    def next(self,skip=False, box=[0]):
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
    def __init__(self, filename,natoms,compressionlevel=9,d=None):
        def_prec=1
        # Precision is always written to the file, so if no dictionary is given,
        # a dictionary containing only the precision is made. If a dictionary
        # _is_ given but the "prec"
        # key is not present, it's added. Same if the key is present, but not 
        # a positive integer. In all these cases, the default value def_prec is
        # used.
        if not d:
            d={"prec":str(def_prec)}
        elif not "prec" in d or not d["prec"].isdigit():
            d["prec"]=str(def_prec)
        precision=int(d["prec"])
        self.prec=pow(10,precision)
        self.natoms=natoms
        self.file=open(filename,'wb')
        self.cctx=zstd.ZstdCompressor(level=compressionlevel)
        self.traj=self.cctx.stream_writer(self.file)
        for k,v in d.items():
            self.traj.write(b"%s=%s\n"%(str(k).encode("utf-8"),str(v).encode("utf-8"))) #I'll attempt converting things to string, but it's your responsibility to ensure that is possible.
        self.traj.write(b"** %d\n"%natoms)
    #Writes the next frame from data, which needs to be an Nx3 numpy array.
    def wnext(self, data, box=[0]):
        if len(data)<self.natoms or len(data[0])<3:
            raise ValueError
        #we center the data on the mean. This reduces a bit the size of the file.
       # mean=data.mean(0)
       # data=data-mean
        for i in range(self.natoms):
            d=data[i]
            self.traj.write(b"%d %d %d\n"%(round(d[0]*self.prec),round(d[1]*self.prec),round(d[2]*self.prec)))
        bstr=b"*\n"
        if len(box)>=9:
            #We center the box on the mean, also
            #for i in (0,1,2):
            #    box[3*i+0]-=mean[0]
            #    box[3*i+1]-=mean[1]
            #    box[3*i+2]-=mean[2]

            bstr=b"* %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n"%(box[0],box[1],
                box[2],box[3],box[4],box[5],box[6],box[7],box[8])
        self.traj.write(bstr)
    def close(self):
            self.traj.close()


