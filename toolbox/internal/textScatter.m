function ax = textScatter(fig_h,bdata,varargin)

% Print text in a Scatter plot.
%
% textScatter(fig_h,bdata) % minimum call
%
%
% INPUTS:
%
% fig_h: (1x1) figure handle
%
% bdata: (Nx2) bidimensional data 
%
%
% Optional INPUTS (parameter):
%
% 'EleLabel': [Nx1] name of the elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single
%   group by default)
%
% 'Multiplicity': [Nx1] multiplicity of each row (1s by default)
%
% 'PlotMult': str
%      'none': do not plot multiplicity (by default)
%      'zaxis': plot multiplicity information in the Z axis.
%      'zsize': plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance with other points where a label is allowed to be 
%   visualized. For a value of 0, all labels are printed, while for a 
%   large value only uncluttered labels are printed. By default Inf is 
%   chosen, where only indices is visualized. 
%
%
% OUTPUTS:
%
% 'ax': [1x4] axis enclosing the text.
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es) and Jesús García (gsus@ugr.es)
% Last modification: 17/Feb/2026
% Dependencies: Matlab R2017b, MEDA v1.9
%
% Copyright (C) 2026  University of Granada, Granada
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

figure(fig_h);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',ones(N,1));  
addParameter(p,'PlotMult','none'); 
addParameter(p,'Multiplicity',ones(N,1)); 
addParameter(p,'BlurIndex',0.3);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
plottype = p.Results.PlotMult;
classes = p.Results.ObsClass;
mult = p.Results.Multiplicity;
blur = p.Results.BlurIndex;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel  = elabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(mult,1) == 1, mult = mult'; end;

% Convert num arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if ~isempty(classes) && isnumeric(classes), classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: paramter ''bdata'' must be N-by-2. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: paramter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

% Get ordering of classes
unique_classes = unique(classes,'stable');
if iscell(classes)
     ord_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
     ord_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end
unique_ord_classes = unique(ord_classes);

ax = axis;
deltax = (ax(2)-ax(1)); 
deltay = (ax(4)-ax(3)); 

c = 0.01; % bias with text
if ~isempty(elabel)
    for i=1:N
        ind = [1:(i-1) (i+1):size(bdata,1)];
        
        diffs = (bdata(ind,:)-bdata(i,:))./[deltax/sqrt(4) deltay/sqrt(3)]; % This creates a set of neighbours closer to the blur index
        dist_vec = sqrt(sum(diffs.^2, 2));
        knn = find(dist_vec < blur); 
        extr = [length(find(diffs(knn,1)>= 0)) length(find(diffs(knn,1)<= 0)) length(find(diffs(knn,2)>= 0)) length(find(diffs(knn,2)<= 0))]; % none right, none left, none up, none down 
        if blur==0 | length(knn)==0 | any(extr==0) % either no blur control activated, no close neighbour or a extreme point wiht neigjbours

            posx = bdata(i,1) + c*deltax;
            posy = bdata(i,2) + c*deltay;
            halign = 'left';
            valign = 'bottom';
            bias = -mean(diffs(knn,:),1);
    
            % Label position based on cuadrant & interference
            if ~isempty(knn) % if there are neighbours (interference), locate from them
                if bias(1) < 0 
                    posx = bdata(i,1) - c*deltax;
                    halign = 'right';
                end
                if bias(2) < 0
                    posy = bdata(i,2) - c*deltay;
                    valign = 'top';
                end
            else % otherwise locate from que quadrant
                if bdata(i,1) < 0 
                    posx = bdata(i,1) - c*deltax;
                    halign = 'right';
                end
                if bdata(i,2) < 0
                    posy = bdata(i,2) - c*deltay;
                    valign = 'top';
                end
            end
    
            if any(strcmp(plottype,{'zaxis','zshape'}))
                t = text(posx, posy, mult(i), strtrim(elabel(i,1)),'VerticalAlignment',valign, 'HorizontalAlignment',halign,'FontSize', 12, 'Interpreter', 'none');
            elseif strcmp(plottype,'zsize')
                t = text(posx, posy, ord_classes(i), strtrim(elabel(i,1)),'VerticalAlignment',valign, 'HorizontalAlignment',halign,'FontSize', 12, 'Interpreter', 'none');
            else
                t = text(posx, posy, strtrim(elabel(i,1)),'VerticalAlignment',valign, 'HorizontalAlignment',halign,'FontSize', 12, 'Interpreter', 'none');
            end
            
            ext = t.Extent; % entend axis if necessary
            text_right_edge = ext(1) + ext(3);
            text_top_edge   = ext(2) + ext(4);
            ax(1) = min(ax(1), ext(1)/1.1);
            ax(2) = max(ax(2), 1.1*text_right_edge);
            ax(3) = min(ax(3), ext(2)/1.1);
            ax(4) = max(ax(4), 1.1*text_top_edge);
            axis(ax);
        end
    end

    txt_objs = findobj(gca, 'Type', 'text');
    n = length(txt_objs);
    padding = 0; % Extra pixels to keep them apart
    iterations = 2; % Repeat a few times to resolve secondary collisions

    for iter = 1:iterations % Do some additional iterations to separate labels
        for i = 1:n
            for j = i+1:n
                % Get current positions
                ext1 = txt_objs(i).Extent;
                ext2 = txt_objs(j).Extent;

                % Calculate overlap in X and Y
                overlapX = min(ext1(1)+ext1(3), ext2(1)+ext2(3)) - max(ext1(1), ext2(1));
                overlapY = min(ext1(2)+ext1(4), ext2(2)+ext2(4)) - max(ext1(2), ext2(2));

                % If both are positive, they overlap
                if overlapX > 0 && overlapY > 0
                    % Move along the smaller overlap to minimize distance
                    if overlapX < overlapY
                        moveX = (overlapX + padding) / 2;
                        signX = sign(ext1(1) - ext2(1));
                        txt_objs(i).Position(1) = txt_objs(i).Position(1) + signX * moveX;
                        txt_objs(j).Position(1) = txt_objs(j).Position(1) - signX * moveX;
                    else
                        moveY = (overlapY + padding) / 2;
                        signY = sign(ext1(2) - ext2(2));
                        txt_objs(i).Position(2) = txt_objs(i).Position(2) + signY * moveY;
                        txt_objs(j).Position(2) = txt_objs(j).Position(2) - signY * moveY;
                    end
                end
            end
        end
    end
end
