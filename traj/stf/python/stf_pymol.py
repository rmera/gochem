#!/usr/bin/env python
from pymol import cmd
import stf
import copy

res2code={"SER":"S", "THR":"T", "ASN":"N", "GLN":"Q", "SEC":"U", "CYS":"C", "GLY":"G", "PRO":"P", "ALA":"A", "VAL": "V", "ILE": "I", "LEU":"L", "MET":"M", "PHE":"F", "TYR":"Y", "TRP":"W", "ARG":"R", "HIS":"H", "LYS":"K", "ASP":"D", "GLU":"E", "HIP":"H", "HID":"H", "HIE":"H"}

def loadstf(filename,objname="",skip=0,begin=0):
    if objname=="":
         objname=filename.replace(".stf","")
    t=stf.rtraj(filename)
    fr=0
    st=0
    skip=int(skip)
    d={}
    while True:
        read=True
        if fr%(skip+1)!=0 or fr+1<int(begin):
            read=False
        fr+=1
        try:
            f=t.next(read)
        except EOFError: #not an actual error
            t.close()
            break
        if not read:
            continue
        b=cmd.get_model(objname)
        #For _some_ reason, PyMOL changes the order of a system when it reads it
        #Here we map the new order to the original one (which, of course, is the one
        #used by the trajectory.
        #We also add back the one-letter code for aminoacids, which PyMOL seem to remove
        #when producing a chempy object from a selection.
        for i,v in enumerate(b.atom):
            d[i]=v.id-1
            try:
                v.resn_code=res2code[v.resn]
            except KeyError:
                v.resn_code=v.resn
        cmd.load_coordset(f, objname,st)
        st+=1
    cmd.center(objname)




cmd.extend("loadstf",loadstf) 
