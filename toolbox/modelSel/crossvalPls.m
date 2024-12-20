function [cumpress,press] = crossvalPls(x,y,varargin)

% Row-wise k-fold (rkf) cross-validation for square-prediction-errors computing in PLS.
%
% [cumpress,press] = crossvalPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional INPUTS (parameter):
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
% 'Option': (str or num) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% cumpress: [Ax1] Cumulative PRESS
%
% press: [AxO] PRESS per variable.
%
%
% EXAMPLE OF USE: Random data with structural relationship, two examples
% of plot.
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% cumpress = crossvalPls(X,Y,'LVs',lvs);
% 
% % Mean centering example
% cumpress = crossvalPls(X,Y,'LVs',lvs,'PreprocessingX',1,'PreprocessingY',1);
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
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LVs',lat'); 
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Option',1);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

lvs = p.Results.LVs;
blocksr = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
opt = p.Results.Option;

% Extract LVs length
A = length(lvs);

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: parameter ''Option'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>2, 'Value Error: parameter ''MaxBlock'' must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain a binary value. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Initialization
cumpress = zeros(length(lvs),1);
press = zeros(length(lvs),O);

rows = rand(1,N);
[a,rind]=sort(rows);
elemr=N/blocksr;

% Cross-validation
        
for i=1:blocksr
    
    indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
    i2 = ones(N,1);
    i2(indi)=0;
    sample = x(indi,:);
    calibr = x(find(i2),:); 
    sampley = y(indi,:);
    calibry = y(find(i2),:); 

    [ccs,av,st] = preprocess2D(calibr,'Preprocessing',prepx);
    [ccsy,avy,sty] = preprocess2D(calibry,'Preprocessing',prepy);
        
    scs = preprocess2Dapp(sample,av,'Scale',st);
    scsy = preprocess2Dapp(sampley,avy,'Scale',sty);
    
    model = simpls(ccs,ccsy,'LVs',0:max(lvs));
    Q = model.yloads;
    R = model.altweights;
    
    for lv=1:length(lvs)
    
        if lvs(lv) > 0
            beta = R(:,1:min(lvs(lv),end))*Q(:,1:min(lvs(lv),end))';
            srec = scs*beta;
            
            pem = scsy-srec;
            
        else % Modelling with the average
            pem = scsy;
        end
        
        press(lv,:) = press(lv,:) + sum(pem.^2,1);
        
    end
end

cumpress = sum(press,2);

%% Show results

if opt == '1' 
    figh = plotVec(cumpress,'EleLabel',lvs,'XYLabel',{'#LVs','PRESS'},'Option','01'); 
end

