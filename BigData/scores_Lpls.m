function [T,TT] = scores_Lpls(Lmodel,test,opt,label,classes)

% Compute and plot scores in PLS for large data. The original 
% paper is Camacho J. Visualizing Big data with Compressed Score Plots: 
% Approach and Research Challenges. Chemometrics and Intelligent Laboratory
% Systems, 2014, 135: 110-125.
%
% scores_Lpca(Lmodel) % minimum call
% [T,TT] = scores_Lpca(Lmodel,test,opt,label,classes) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lvs: [1x1] number of LVs. 
%       Lmodel.centr: (LxM) centroids of the clusters of observations
%       Lmodel.multr: (Lx1) multiplicity of each cluster.
%       Lmodel.class: (Lx1) class associated to each cluster.
%       Lmodel.av: [1xM] sample average according to the preprocessing method.
%       Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%       Lmodel.weight: [1xM] weight applied after the preprocessing method.
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% opt: (str or num) options for data plotting: binary code of the form 'abcd' for:
%       a:
%           0: no plots
%           1: plot scores
%       b:
%           0: scatter plot of pairs of PCs 
%           1: bar plot of each single PC
%       c:
%           0: plot calibration and test data
%           1: plot only test data 
%       d:
%           00: plot multiplicity info in the size of the markers.
%           01: plot multiplicity info in the form of the markers.
%           10: plot multiplicity information in the Z axis.
%           11: plot multiplicity info in the size of the markers and 
%               classes in Z-axis
%           
%   By deafult, opt = '10000'. If less than 5 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0, c=0, d=00. 
%   If a=0, then b, c and d are ignored. If b=1, then d is ignored.
%
% label: [Lx1] name of the test observations (numbers are used by default)
%
% classes: [Lx1] groups in test for different visualization (a single group 
%   by default)
%
%
% OUTPUTS:
%
% T: [LxA] calibration scores.
%
% TT: [NxA] test scores.
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Lmodel = Lmodel_ini(X,Y);
% Lmodel.lvs = 1:3;
% T = scores_Lpca(Lmodel);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% X = simuleMV(n_obs,n_vars,8);
% Y = 0.1*randn(n_obs,2) + X(:,1:2);
% Lmodel = Lmodel_ini(X,Y);
%
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(X)*(n_obst-1)/(n_obs-1));
%
% Lmodel.lvs = 1;
% scores_Lpca(Lmodel,test);
% Lmodel.lvs = 1:2;
% [T,TT] = scores_Lpca(Lmodel,test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/May/17.
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
% 
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

check_Lmodel(Lmodel);

N = Lmodel.nc;
M = size(Lmodel.XX, 2);

if nargin < 2, test = []; end;
L = size(test, 1);

if nargin < 3 || isempty(opt), opt = '100'; end; 

A = length(Lmodel.lvs);

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
while length(opt)<5, opt = strcat(opt,'0'); end
if opt(3) == '1',
    K = L;
else
    K = N+L;
end

if nargin < 4 || isempty(label), 
    if  opt(3) == '1',
        label = cellstr(num2str((1:L)'));
    elseif isempty(Lmodel.obs_l),
        label = cellstr(num2str([1:N 1:L]'));
    else
        if L
            lb1 = cellstr(num2str((1:L)'));
            label = {Lmodel.obs_l{:} lb1{:}};
        else
            label = Lmodel.obs_l;
        end
    end
else
    if  opt(3) == '0',
        if isempty(Lmodel.obs_l),
            lb1 = cellstr(num2str((1:N)'));
            label = {lb1{:} label{:}};
        else
            label = {Lmodel.obs_l{:} label{:}};
        end
    end
end

if nargin < 5 || isempty(classes),
    if opt(3) == '1', 
        classes = ones(L,1); 
    else
        classes = [ones(N,1);2*ones(L,1)];  
    end
elseif opt(3) == '0' && length(classes)==L,
        classes = [Lmodel.class;2*classes];
end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 1sr argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)==5, 'Dimension Error: 3rd argument must be a string or num of maximum 5 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 4th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[beta,W,P,Q,R] = Lpls(Lmodel);
T = Lmodel.centr*R;

if ~isempty(test)
    testcs = preprocess2Dapp(test,Lmodel.av,Lmodel.sc,Lmodel.weight);
    TT = testcs*R;
else
    TT = [];
end

%% Show results

if opt(1) == '1',
     
    if opt(3) == '0'
        ttt = [T;TT];
        mult = [Lmodel.multr;ones(size(TT,1),1)];
    else
        ttt = TT;
        mult = ones(size(TT,1));
    end
    
    if length(Lmodel.lvs) == 1 || opt(2) == '1',
        for i=1:length(Lmodel.lvs)
            plot_vec(ttt(:,i), label, classes, {'',sprintf('Compressed Scores PC %d',Lmodel.lvs(i))}, [], [], [], mult);
        end
    else
        for i=1:length(Lmodel.lvs)-1,
            for j=i+1:length(Lmodel.lvs),
                plot_scatter([ttt(:,i),ttt(:,j)], label, classes, {sprintf('Scores PC %d',Lmodel.lvs(i)),sprintf('Scores PC %d',Lmodel.lvs(j))}, [], strcat('1',opt(4:5)), mult);
            end
        end
    end
end



