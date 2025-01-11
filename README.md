To the long life of the Ven. Khenpo Phuntsok Tenzin Rinpoche


## goChem: a Library for Computational Chemistry written in the Go programming language. 

Check out www.gochem.org for more information.



Design goals:

* Simplicity.
* Readability.
* Fast and light
* Concurrent when possible and necessary
* Easy to extend
* Useful for Computational Chemistry/Biochemistry 
at a classical and QM levels.


See the Wiki por goChem's current capabilities
(https://github.com/rmera/gochem/wiki)


If you use goChem, or a program based on goChem (such as goMD and goAnalize) 
we ask you to support the library by citing goChem's paper:

Domínguez, M., A. Jiménez, V., Savasci, G., Araya-Osorio, R., Pesonen, J. and Mera-Adasme, R. (2025), goChem: A Composable Library for Multi-Scale Computational Chemistry Data Analysis. J Comput Chem, 46: e70004. https://doi.org/10.1002/jcc.70004


Currently, gochem is licensed under LGPL2.1. This might change 
towards BSD in the future. Meanwhile, if you want to use some 
of this for a BSD-licensed project, contact the developer.



goChem uses gonum (https://www.gonum.org) for matrix
operations. The user can choose to have gonum backed by a pure-go 
(and some assembly) implementation of BLAS (which allows not having 
runtime dependencies) or by the somewhat more efficient cBLAS. 

Reading xtc files requires the xdrfile
library from Gromacs (www.gromacs.org) 

All dependencies of goChem are open source.


LICENSE

Copyright 2012 Raul Mera rmeraaatacademicosdotutadotcl


This program, including its documentation, 
is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation; either version 2.1 of the 
License, or (at your option) any later version.
	  
This program and its documentation is distributed in the hope that 
it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the GNU General Public License for more details.
		    
You should have received a copy of the GNU Lesser General 
Public License along with this program. If not, see 
<http://www.gnu.org/licenses/>. 


The mascot is a modification by Sebastian Franchini of the
Go language mascot by Renee French.


