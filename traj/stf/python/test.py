#!/usr/bin/env python3

import stf

name="out.stf"
oname="testout.stf"

r=stf.rtraj(name)

w=stf.wtraj(oname,r.get_natoms())

while True:
    try:
        d=r.next()
    except EOFError:
        break
    w.wnext(d)

r.close()
w.close()





