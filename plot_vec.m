
function fig_h = plot_vec(vec,elabel,classes,xylabel,lcont,opt,vlabel)

% Bar or line plot.
%
% fig_h = plot_vec(vec) % minimum call
% fig_h = plot_vec(vec,elabel,classes,xylabel,lcont,opt,vlabel) % complete call
%
%
% INPUTS:
%
% vec: [NxM] vector/s to plot. 
%
% elabel: [Nx1] name of the vector elements (numbers are used by default)
%
% classes: [Nx1, str(N), {N}] groups for different visualization (a single 
%   group by default)
%
% xylabel: {2} xlabel and ylabel (nothing by default)
%
% lcont: [NxL or Lx1] L control limits (nothing by default)
%
% opt: (str or num) options for data plotting.
%       0: line plot
%       1: bar plot (default)
%
% vlabel: [Mx1] name of the vectors (numbers are used by default)
%
%
% OUTPUTS:
%
% fig_h: [1x1] figure handle
%
%
% EXAMPLE OF USE: To plot three lines with constant control limits:
%
% fig_h = plot_vec(randn(100,3),[],[],{'Functions','Time'},[1, -1, 3], 0);
%
%
% EXAMPLE OF USE: with labels and classes in observations and variable limit:
%
% fig_h = plot_vec(randn(5,3),{'one','two','three','four','five'},[1 1 1 2 2],{[],'Functions'},randn(5,1));
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 19/Apr/2016
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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


%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if size(vec,1) == 1,     vec = vec'; end;
N = size(vec, 1);
M = size(vec, 2);
if nargin < 2 || isempty(elabel), elabel = 1:N; end;
if nargin < 3 || isempty(classes), classes = []; end;
if nargin < 4 || isempty(xylabel), xylabel = {'',''}; end;
if nargin < 5 || isempty(lcont),  lcont = []; end;
if nargin < 6 || isempty(opt),  opt = '1'; end;
if nargin < 7 || isempty(vlabel),  vlabel = 1:M; end;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(vlabel,1)  == 1, vlabel = vlabel'; end;

% Convert int arrays to str
if ~isempty(elabel) && isnumeric(elabel), 
    vecn = elabel; 
    veci = 1:length(vecn);
    if length(elabel)>2,  
        max_lab = 30; % limit the number of labels displayed
        ini = 2;
        stepN = [];
        while isempty(stepN),
            lenv = length(vecn(ini:end));
            div = 1:(lenv-1);
            div = div(rem(lenv,div)==0);
            stepN = div(find(div>lenv/max_lab,1));
            ini = ini+1;
        end
        veci = 1:(lenv+ini-2);
        veci = veci(round([1 (ini-2+stepN):stepN:end]));
    end
    for i=veci,
        labele{i} = num2str(vecn(i));
    end
    elabel=labele'; 
end
if ~isempty(vlabel) && isnumeric(vlabel), vlabel=num2str(vlabel); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;
if ischar(vlabel),  vlabel = cellstr(vlabel); end;

% Validate dimensions of input data
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: 3rd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: 4th argument must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (isequal(size(lcont,1), N) || isequal(size(lcont,2), 1), 'Dimension Error: 5th argument must be N-by-L or L-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (isequal(size(opt), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(vlabel), assert (isequal(size(vlabel), [M 1]), 'Dimension Error: 7th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); end;
 
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 6th argument must contain a binary value. Type ''help %s'' for more info.', routine(1).name);
   
% Convert constant limits in vectors
if ~isempty(lcont) && ~isequal(size(lcont,1), N), lcont = (lcont*ones(1,N))'; end;
    
% Exception: bar plot with multivariate vec and one-observation class  
if ~opt && ~isempty(classes) && size(vec, 2)>1,
    unique_classes = unique(classes);
    assert (min(hist(classes,unique(classes)))>1, 'Exception: Cannot visualize a multivariate bar plot with one-observation classes. Try setting the 6th argument to 1.'); 
end;
    

%% Main code

% Create figure window
fig_h = figure;
hold on;

% Preprocess classes to force them start with 1, 2...n,
unique_classes = unique(classes);
if iscell(classes)
    classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
    classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end

% Plot vectors

if ~isempty(classes)
    unique_classes = unique(classes);
    color_list = hsv(length(unique_classes));
    if opt == '0',
        plot(vec,'k','HandleVisibility', 'off');
    end  
    for i=1:length(unique_classes)
        ind = classes == unique_classes(i);
        if opt == '0',
            plot(find(ind), vec(ind,:), 'Color', 'none', 'Marker','O', 'MarkerFaceColor', color_list(i,:), 'DisplayName', num2str(unique_classes(i)));
        else 
            bar(find(ind), vec(ind,:), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', num2str(unique_classes(i)));
        end
    end
else
    color_list = hsv(M);
    for i=1:M,
        if opt == '0',
            plot(vec(:,i), 'LineWidth', 2, 'Color', color_list(i,:), 'DisplayName', vlabel{i});
        else
            bar(vec(:,i), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', vlabel{i});
        end
    end
end

% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    for i=1:size(lcont,2),
        a = [lcont(:,i)';lcont(:,i)'];
        plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
    end
end    

% Get axes handler
axes_h = get(fig_h,'Children');
if length(axes_h)>1, axes_h = axes_h(1); end;

% Set ticks and labels
label_length = max(cellfun('length', elabel));
label_size = 300/(length(find(~cellfun('isempty', elabel)))*label_length);
set(axes_h, 'FontSize', max(min(14,round(label_size)), 10));
if ~isempty(elabel)
    stepN = ceil(0.2*N/label_size);
    if stepN==1,
        vals = 1:N;
        set(axes_h,'XTick',vals);
        set(axes_h,'XTickLabel',elabel(vals));
    else
        set(axes_h,'XTickMode','auto');
        set(axes_h, 'FontSize', 14);
    end
end

if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 16);
    ylabel(xylabel{2}, 'FontSize', 16);
end

% Set axis
axis tight
ax = axis;
axis auto
ax2 = axis;
axis([ax(1:2) ax2(3:4)])

legend off
box on
hold off

       