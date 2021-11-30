#!/usr/bin/env python
from pymol import cmd
import stf
import copy

res2code={"SER":"S", "THR":"T", "ASN":"N", "GLN":"Q", "SEC":"U", "CYS":"C", "GLY":"G", "PRO":"P", "ALA":"A", "VAL": "V", "ILE": "I", "LEU":"L", "MET":"M", "PHE":"F", "TYR":"Y", "TRP":"W", "ARG":"R", "HIS":"H", "LYS":"K", "ASP":"D", "GLU":"E", "HIP":"H", "HID":"H", "HIE":"H"}

def loadstf(filename,sele="all",objname="",skip=0,begin=0):
    if objname=="":
        objname=filename.replace(".stz","")
    t=stf.traj(filename)
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
        except: #not an actual error
            break
        if not read:
            continue
        b=cmd.get_model(sele)
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
        for i in range(len(b.atom)):
            b.atom[i].coord[0]=f[d[i]][0]
            b.atom[i].coord[1]=f[d[i]][1]
            b.atom[i].coord[2]=f[d[i]][2]
        cmd.load_model(b,objname,finish=1,discrete=0)
        st+=1
    cmd.center(objname)
    cmd.hide("wire")
    cmd.show("cartoon")
    cmd.set("seq_view_discrete_by_state",0)




cmd.extend("loadstf",loadstf) 
