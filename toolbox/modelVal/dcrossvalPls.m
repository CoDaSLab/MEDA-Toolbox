function [Qm,Q,lvso,MSE] = dcrossvalPls(x,y,varargin)

% Row-wise k-fold (rkf) double cross-validation for PLS. The algorithm uses
% repetitions of the dCV loop to estimate the stability: see Szymanska, E., 
% Saccenti, E., Smilde, A.K., Weterhuis, J. Metabolomics (2012) 8: 3.
%
% Qm = dcrossvalPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional INPUTS:
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% 'MaxBlock': [1x1] maximum number of blocks of samples (N by default)
%
% 'PreprocessingX': [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% 'PreprocessingY': [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% 'Repetitions': [1x1] number of repetitions for stability.
%
% 'Option': [1x1] options for data plotting
%       0: no plots
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% Qm: [1x1] Mean Goodness of Prediction
%
% Q: [rep x 1] Goodness of Prediction
%
% lvso: [rep x blocksr] optimum number of LVs in the inner loop
%
% MSE: [1x1] Mean Square Error
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% Q = dcrossvalPls(X,Y,'LVs',lvs,'Repetitions',5);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 20/Nov/2024
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LVs',lat'); 
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Repetitions',10);
addParameter(p,'Option',1);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

lvs = p.Results.LVs;
blocksr = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
rep = p.Results.Repetitions;
opt = p.Results.Option;

% Extract LVs length
A = length(lvs);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(rep), [1 1]), 'Dimension Error: parameter ''Repetitions'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: parameter ''Option'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>3, 'Value Error: parameter ''MaxBlock'' must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

for j=1:rep
    % Cross-validation
    
    rows = rand(1,N);
    [a,rind]=sort(rows);
    elemr=N/blocksr;
    
    for i=1:blocksr
        indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
        i2 = ones(N,1);
        i2(indi)=0;
        val = x(indi,:);
        rest = x(find(i2),:);
        valy = y(indi,:);
        resty = y(find(i2),:);
        
        cumpress = crossvalPls(rest,resty,'LVs',lvs,'MaxBlock',blocksr-1,'PreprocessingX',prepx,'PreprocessingX',prepy,'Option',0);
        
        lvso(j,i) = lvs(find(cumpress==min(cumpress),1));
        
        [ccs,av,st] = preprocess2D(rest,'Preprocessing',prepx);
        [ccsy,avy,sty] = preprocess2D(resty,'Preprocessing',prepy);
        
        vcs = preprocess2Dapp(val,av,'Scale',st);
        vcsy = preprocess2Dapp(valy,avy,'Scale',sty);
        
        model = simpls(ccs,ccsy,'LVs',1:lvso(i));
        srec = vcs*model.beta;
        
        Qu(i) = sum(sum((vcsy-srec).^2));
        Qd(i) = sum(sum(vcsy.^2));
    end
    
    Q(j) = 1-sum(Qu)/sum(Qd);
    
    SSE(j) = sum(Qu);
    
end

Qm = mean(Q);
MSE = mean(SSE);

%% Show results

if opt == 1
   figh = plotVec(Q,'XYLabel',{'#Repetition','Goodness of Prediction'},'Option','11'); 
end

