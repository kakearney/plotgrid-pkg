function h = plotgrid(fun, varargin) %x, y, rowname, colname, varargin)
%PLOTGRID Plots data from cell array in grid pattern of axes
%
% h = plotgrid(fun, in1, in2, in3, ...)
% h = plotgrid(fun, in1, in2, in3, ..., rowname, colname)
% h = plotgrid(fun, in1, in2, in3, ..., rowname, colname, prop, val, ...)
%
% Input arguments:
%
%   fun:        function handle of plotting function.  Must have an
%               explicit number of input arguments, so if using something
%               that can take a variable number, use an anonymous function
%               to specify (for example, @(y) plot(y) or @(x,y) plot(x,y)
%               instead of @plot).
%
%               or this can be the string 'setup', which just sets up the
%               axes without plotting any data, 
%
%   in#:        Equally-sized 2D cell arrays (n x m), with values
%               corresponding to the input for the function handle. (If
%               'setup', then a single cell array with the desired number
%               of rows and columns should be passed).  The figure will
%               have the same number of rows and columns as these cell
%               arrays.  
%
%   rowname:    1 x n cell array of strings,  labels for each row.  If none
%               desired, pass empty array
%
%   colname:    1 x m cell array of strings, labels for each column.  If none
%               desired, pass empty array.
%
% Optional input variables:
%
%   'staggery': if included, left side of axes in every other row will be
%               indented, with some overlap in vertical extent (see
%               plotses.m)  
%
%   Also, any of the positioning parameters accepted by subaxis can be
%   passed as parameter.value pairs.
%
% Output variables:
%
%   h:          1 x 1 structure array of handles:
%
%               fig:    figure
%               
%               ax:     n x m array, axis handles
%
%               rlab:   n x 1 array, row label text handles
%
%               clab:   m x 1 array, column label text handles

% Copyright 2013 Kelly Kearney

% TODO row and column labels don't work if you include x and y labels on
% each axis


%--------------------------
% Parse input
%--------------------------

if ischar(fun) && strcmp(fun, 'setup')
    setup = true;
    nin = 1;
else
    setup = false;

    if ~isa(fun, 'function_handle')
        error('First input must be a function handle');
    end
        
    % Determine number of inputs passed to plotting function

    nin = nargin(fun);
    if nin == -1
        error('Plotting function must accept an explicit number of input arguments');
    end 
    
end

funin = varargin(1:nin);

% Expand any plotting function inputs that aren't cell arrays

iscl = cellfun(@iscell, funin);
sz = cellfun(@size, funin(iscl), 'uni', 0);
sz = unique(cat(1, sz{:}), 'rows');
if size(sz,1) ~= 1
    error('Function inputs must be equal-sized cell arrays or arrays');
end

for iin = 1:nin
    if ~iscl(iin)
        temp = funin{iin};
        funin{iin} = cell(sz);
        [funin{iin}{:}] = deal(temp);
    end
end

% Get other arguments

if nargin == (nin + 1)
    rowname = [];
    colname = [];
    axprops = cell(0);
elseif nargin == (nin + 2)
    rowname = varargin{nin + 1};
    colname = [];
    axprops = cell(0);
elseif nargin == (nin + 3)
    rowname = varargin{nin + 1};
    colname = varargin{nin + 2};
    axprops = cell(0);
elseif nargin > (nin + 3)
    rowname = varargin{nin + 1};
    colname = varargin{nin + 2};
    axprops = varargin(nin+3:end);
end

if ~isempty(rowname) && (~isvector(rowname) || length(rowname)~=sz(1))
    error('rowname wrong size');
end

if ~isempty(colname) && (~isvector(colname) || length(colname)~=sz(2))
    error('colname wrong size');
end

% Check for staggered axis indicators

isstr = cellfun(@ischar, axprops);
staggerx = cellfun(@(x) ischar(x) && strcmp(x, 'staggerx'), axprops);
axprops(staggerx) = [];
staggerx = any(staggerx);

staggery = cellfun(@(x) ischar(x) && strcmp(x, 'staggery'), axprops);
axprops(staggery) = [];
staggery = any(staggery);

if staggery & setup
    error('Cannot set up staggered axes without a plotting function');
end

% Check for sparse indicator

usesparse =  cellfun(@(x) ischar(x) && strcmp(x, 'sparse'), axprops);
axprops(usesparse) = [];
usesparse = any(usesparse);

if usesparse
    for ii = 1:nin
        isemp{ii} = cellfun('isempty', funin{ii});
    end
    isemp = cat(3, isemp{:});
    isemp = any(isemp, 3);
else
    isemp = false(size(funin{1}));
end

% Check for figure properties

isstr = cellfun(@ischar, axprops);
isfigprop = cellfun(@(x) ischar(x) && strcmp(x, 'figprop'), axprops);
if any(isfigprop)
    idx = find(isfigprop);
    figprop = axprops{idx+1};
    axprops([idx idx+1]) = [];
else
    figprop = cell(0);
end

% Check for output cell

isout = cellfun(@iscell, axprops);
if any(isout)
    outname = axprops{isout};
    axprops(isout) = [];
    nout = length(outname);
else
    nout = 0;
end


%--------------------------
% Create axes
%--------------------------

h.fig = figure(figprop{:});

nrow = sz(1);
ncol = sz(2);

if ~staggerx & ~staggery

    for irow = 1:nrow
        for icol = 1:ncol
            if isemp(irow,icol)
%                 h.ax(irow,icol) = subaxis(nrow, ncol, icol, irow, axprops{:});
%                 set(h.ax(irow,icol), 'visible', 'off');
            else
                h.ax(irow,icol) = subaxis(nrow, ncol, icol, irow, axprops{:});
            end
        end
    end
    
end

if staggery
    
    htemp = plotgrid(fun, funin{:}, [], [], axprops{:});
    axis(htemp.ax(:), 'tight'); % Trying to eliminate space on x axis, since can't adjust after staggered
    xlim = get(htemp.ax, 'xlim');
    close(htemp.fig);
    xlim = minmax(cat(1, xlim{:}));
    
    for icol = 1:ncol
        axtemp(icol) = subaxis(1, ncol, icol, 1, axprops{:});
        [hl, h.ax(:,icol)] = plotses(xlim, rand(2,nrow));
        delete(hl);
    end
    h.ax = flipud(h.ax);
    
end

if staggerx
    error('Can''t stagger x axes yet... look into plotses bug');
end

%--------------------------
% Plot data
%--------------------------

for iout = 1:nout
    h.(outname{iout}) = cell(nrow, ncol);
end


if ~setup
    
    for irow = 1:nrow
        for icol = 1:ncol

            if ~isemp(irow,icol)


%                 axes(h.ax(irow,icol));
%                 hold on;
%                 
                set(h.fig, 'currentaxes', h.ax(irow,icol));
                set(h.ax(irow,icol), 'nextplot', 'add');
                
                invar = cellfun(@(x) x{irow,icol}, funin, 'uni', 0);
                
                tmp = cell(1,nout);
                
                if nout > 0
                    [tmp{:}] = fun(invar{:});
                    for iout = 1:nout
                        h.(outname{iout}){irow,icol} = tmp{iout};
                    end
                else
                    fun(invar{:});
                end
                

            end
        end
    end
    
end

%--------------------------
% Add row and column labels
%--------------------------

Subax = get(h.fig, 'UserData');

xpos = linspace(Subax.MarginLeft, 1-Subax.MarginRight, Subax.cols+1);
ypos = linspace(1-Subax.MarginTop, Subax.MarginBottom, Subax.rows+1);
left = xpos(1)/2;
top = (ypos(1)+1)/2;
xpos = (xpos(1:end-1) + xpos(2:end))./2;
ypos = (ypos(1:end-1) + ypos(2:end))./2;

labelax = axes('position', [0 0 1 1], 'visible', 'off');
set(labelax, 'xlim', [0 1], 'ylim', [0 1], 'handlevisibility', 'off');

if ~isempty(rowname)
%     for irow = 1:nrow
    h.rlab = text(left*ones(size(ypos)), ypos, rowname, 'parent', labelax, 'horiz', 'center', 'rotation', 90);
%         h.rlab(irow) = suplabel('axes', h.ax(irow,:), 'ylabel', rowname{irow});
%     end
end
   
if ~isempty(colname)
    h.clab = text(xpos, top*ones(size(xpos)), colname, 'parent', labelax, 'horiz', 'center');
    
%     for icol = 1:ncol
% %         h.clab(icol) = suplabel('axes', h.ax(:, icol), 'xlabel', colname{icol});
%         h.clab(icol) = suplabel('axes', h.ax(:, icol), 'title', colname{icol});
%     end
end


    
