function L = leveragesLpca(Lmodel,varargin)

% Compute and plot the leverages of variables in the PCA model for large data.
%
% L = leveragesLpca(Lmodel) % minimum call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lvs: (1x1) number of PCs.
%       Lmodel.vclass: [Mx1] class associated to each variable.
%       Lmodel.varl: {ncx1} label of each variable.
%
% Optional INPUTS (parameters):
%
% 'Option': (str or num) options for data plotting.
%       0: no plots.
%       1: plot bar plot of leverages (default)
%
%
% OUTPUTS:
%
% L: [Mx1] leverages of the variables
%
%
% EXAMPLE OF USE: Random leverages
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 1:3;
% L = leveragesLpca(Lmodel);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 21/Nov/2024
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
[ok, Lmodel] = checkLmodel(Lmodel);
% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Option','1');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
opt = p.Results.Option;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Validate dimensions of input data
assert (isequal(size(opt), [1 1]), 'Dimension Error: 2nd argument must be a string or num of maximum 2 bits. Type ''help %s'' for more info.', routine(1).name);
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel = Lpca(Lmodel);
P = Lmodel.loads;

L = diag(P*P');

%% Show results

if opt == '1' 
    plotVec(L,'EleLabel',Lmodel.varl,'ObsClass',Lmodel.vclass,'XYLabel',{'Variables','Leverages'});
end