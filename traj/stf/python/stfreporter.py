# Copyright (c) 2021, Raul Mera A.

#This program is based on the DCD reporter dcdreporte.py
#From DCD. The copyright note and license for that program is included below.
#The included license applies to this program, also.

"""
This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

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


from __future__ import absolute_import
import numpy as np
import stf
from openmm.unit import nanometer, nanometers

class STFReporter(object):
    """STFReporter outputs a series of frames from a Simulation to an STF file.

    To use it, create a STFReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, natoms, append=False, enforcePeriodicBox=None, d={}):
        """Create a STFReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
                Only for compatibility. Appending not supported (yet!).    
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        natoms: The number of atoms per frame of the trajectory.
        d: a string:string dictionary that will be added to the trajectory header
        """
        self._reportInterval = reportInterval
        self._append = False #not actually supported
        self._enforcePeriodicBox = enforcePeriodicBox
        d["description"]="openmm trajectory"
        self._stf=stf.wtraj(file,natoms,d=d)

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        fb=state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometers)*10
        #np.array(cq[0].Quantity(),)
      #  npbox=[np.array([fb[0].x,fb[0].y,fb[0].z]),np.array([fb[1].x,fb[1].y,fb[1].z]),np.array([fb[2].x,fb[2].y,fb[2].z])]
       # print(fb,npbox) ##############   
        fbox=[]
        for i in fb:
           fbox.append(i[0])
           fbox.append(i[1])
           fbox.append(i[2])
        #this will be somewhat slow, and the np conversion might be unnecessary, but yeah, it will have to do for now.
        self._stf.wnext(state.getPositions(asNumpy=True).value_in_unit(nanometers)*10, box=fbox)

    def __del__(self):
        self._stf.close()
