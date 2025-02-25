﻿Big Data Extension for the MEDA Toolbox for its use in MATLAB.

Coded by: José Camacho (josecamacho@ugr.es), Jesús García (gsus@ugr.es)
Last modification of this document: 22/Nov/2024


Please, acknowledge the use of this software by refercing it: "Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E. Multivariate 
Exploratory Data Analysis (MEDA) Toolbox. Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49-57, available at 
https://github.com/josecamachop/MEDA-Toolbox" Also, please check the documentation of the routines used for more related references. 


Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have 
suggestions, please contact josecamacho@ugr.es or gsus@ugr.es


Items in the folder:

- Big Data toolbox routines:

	- projection models: Lpca.m (PCA), Lpls.m (PLS), Lgpca.m (GPCA)

	- exploratory tools: scoresLpca.m & scoresLpls.m (Score plots), loadingsLpca.m & loadingsLpls.m (Loading plots), 
		medaLpca.m & medaLpls.m (MEDA), omedaLpca.m & omedaLpls.m (oMEDA) mspcLpca.m & mspcLpls.m & mspcAdicov (MSPC), 
		leveragesLpca.m & leveragesLpls.m (leverages of variables)

	- tools to select the number of LVs: varLpca.m & varLpls.m (Variance plots) 

	- auxiliary routines: preprocess2Di.m, psc.m 

	- Lmodel management routines: iniLmodel.m, checkLmodel.m

	- Lmodel update routines: updateEwma.m, updateIterative.m

	- file management routines: addData.m, addDndices.m, cfilesys.m, readData.m, readIndices.m, VCfile.m 
  

You may find an example of use in the Examples\Networkmetric\KDD directory.

To start using the routines, you need your data in a specific format. Data should be splitted in a number of ".mat" files stored in a
specific directory. These files are Matlab storage files where you store partitions of the data, with 10000 observations at most. The data 
corresponding to each file should be in variable 'x', output/quality variables for the same observations should be named 'y' and a numbered 
class for the observations should be stored in variable 'class'. Variables 'x' and 'class' should be in the mat file, while 'y' is optional. 
For an example of this, please have a look at the data directory in the KDD example. 

Copyright (C) 2024  Universidad de Granada
Copyright (C) 2024  Jose Camacho
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.