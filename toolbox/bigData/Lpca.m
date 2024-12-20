function Lmodel = Lpca(Lmodel)
  
% PCA for large data. 
%
% Lmodel = Lpca(Lmodel) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.lvs: [1x1] number of PCs A.
%       Lmodel.centr: (LxM) centroids of the clusters of observations
%
%
% OUTPUTS:
%
% Lmodel: (struct Lmodel) output model
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 0:10;
% Lmodel = Lpca(Lmodel)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
[ok,Lmodel] = checkLmodel(Lmodel);


%% Main code

[P,d2] = eig(Lmodel.XX);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
P = P(:,inddd(Lmodel.lvs));
T = Lmodel.centr*P;
Lmodel.loads = P;
Lmodel.scores = T;
Lmodel.var = sum(dd);
Lmodel.sdT = dd(inddd(Lmodel.lvs));
Lmodel.type = 'PCA';
