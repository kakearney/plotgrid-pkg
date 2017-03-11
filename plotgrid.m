function h = plotgrid(varargin) %x, y, rowname, colname, varargin)
%PLOTGRID Plots data from cell array in grid pattern of axes
%
% h = plotgrid(...)
%
% This function sets up a grid of axes, returning handles in arrays with
% the same geometry as the axes themselves, making it easier to reference
% plotted objects.  It allows for several customization options, including
% the margins around and spacing between axes, optional row and column
% labels, and options to offset axes from each other.
%
% Optional input arguments (passed as parameter/value pairs):
%
%   size:           1 x 2 vector, indicating the number of rows (n) and
%                   columns (m), respectively, to add to the figure.
%                   Default: [1 1]
%
%   function:       cell array, where first element is a function handle to
%                   a plotting function, and remaining elements are n x m
%                   cell arrays of input data to that function.  A n x m
%                   grid of axes will be created, with data plotted to
%                   those axes corresponding to the respective input in the
%                   cell arrays. For example, an input of {@(y) plot(y),
%                   {rand(10,1) 1:10}} will plot the two line plots in
%                   side-by-side axes.  If included, the size of the input
%                   data cell arrays overrides the 'size' input to set the
%                   axes geometry.  
%                   Default: {}
%
%   staggery:       scalar. If non-zero, the y-axis of every other row
%                   (starting second from the bottom) will be offset by the
%                   specified fraction (as a fraction of the axis width,
%                   see offsetaxis.m).  This can be combined with a
%                   negative vertical spacing value to create overlapping
%                   axes. 
%                   Default: 0
%
%   staggerx:       scalar.  If non-zero, the x-axis of every other column
%                   (starting second from the left) will be offset by the
%                   specific fraction (as a fraction of the axis height).
%                   This can be combined with a negative horizontal spacing
%                   to create overlapping axes.
%                   Default: 0
%
%   rowlabel:       n x 1 cell array of strings, text labels to be applied
%                   like a y-axis label to the left of each row in the
%                   grid. If empty, no labels are added.  Rows are numbered
%                   from top to bottom.    
%                   Default: {}
%
%   collabel:       m x 1 cell array of strings, text labels to be applied
%                   like titles to each column in the grid. If empty, no
%                   labels are added.
%                   Default: {} 
%
%   rowlabeloffset: scalar, distance row labels should be offset from the
%                   grid of axes, expressed as a fraction of the entire
%                   axis grid width (right of column m - left of column 1)
%                   Default: 0.05
%
%   collabeloffset: scalar, distance column labels should be offset from the
%                   grid of axes, expressed as a fraction of the entire
%                   axis grid height (top of column 1 - botoom of column n)
%                   Default: 0.05
%
%   figprop:        cell array holding parameter/value pairs of figure
%                   properties to apply to the newly-created figure. 
%
%   outputs:        cell array of strings.  If the 'function' input is
%                   supplied, this can be used to save any outputs
%                   returned by the called function.  The length of this
%                   cell array should match the number of outputs of the
%                   returned by the supplied function.  Each string will be
%                   used as the name of a field in the output structure,
%                   and will hold an n x m cell array holding the output
%                   associated with each axis.
%
%   In addition to the above parameters, any parameter accepted by the
%   subaxis function (spacing, padding, margin), or their abbreviations,
%   with the exception of 'Holdaxis', can be included. 
%
% Output variables:
%
%   h:          1 x 1 structure array of handles:
%
%               fig:    figure
%               
%               ax:     n x m array, axis handles.  Geometry of this array
%                       matches the geometry of the axes on the plot, e.g.
%                       h.ax(1,:) refers to the top row of axes.
%
%               rlab:   n x 1 array, row label text handles
%
%               clab:   m x 1 array, column label text handles
%
%               yax:    floor(n/2) x m array, axis handles to offset y
%                       axes.  These are linked to the axes in
%                       h.ax(end-1:-2:1,:).  These are for decoration only;
%                       see offsetaxis.m for details.
%
%               xax:    n x floor(m/2) array, axis handles to offset x
%                       axes.  These are linked to the axes in
%                       h.ax(:,2:2:end).  These are for decoration only;
%                       see offsetaxis.m for details.

% Copyright 2013-2017 Kelly Kearney


%--------------------------
% Parse input
%--------------------------

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('size', [1 1], @(x) validateattributes(x, {'numeric'}, {'integer', 'positive', 'size', [1 2]}));
p.addParameter('function', {}, @(x) validateattributes(x, {'cell'},{}));
p.addParameter('staggery', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('staggerx', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('rowlabel', {}, @(x) validateattributes(x, {'cell'}, {'vector'}));
p.addParameter('collabel', {}, @(x) validateattributes(x, {'cell'}, {'vector'}));
p.addParameter('rowlabeloffset', 0.05, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('collabeloffset', 0.05, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('figprop', {}, @(x) validateattributes(x, {'cell'},{}));
p.addParameter('outputs', {}, @(x) validateattributes(x, {'cell'},{}));
p.parse(varargin{:});
Opt = p.Results;
SubaxOpt = p.Unmatched;

% If function input is provided, parse this

setup = isempty(Opt.function);
if ~setup
    if ~isa(Opt.function{1}, 'function_handle')
        error('First element of function input must be a function handle');
    end
    
    nin = nargin(Opt.function{1});
    fun = Opt.function{1};
    funin = Opt.function(2:end);
    
    % Determine size of axis grid based on cell array inputs to function.
    % All input cell arrays should be the same size; any non-cell array
    % input is assumed to apply to all axes.
    
    iscl = cellfun(@(x) iscell(x) && ndims(x) <= 2, funin);    
    
    sz = cellfun(@size, funin(iscl), 'uni', 0);
    sz = unique(cat(1, sz{:}), 'rows');
    if size(sz,1) ~= 1
        error('Function inputs must be equal-sized cell arrays');
    end
    for iin = 1:nin
        if ~iscl(iin)
            temp = funin{iin};
            funin{iin} = cell(sz);
            [funin{iin}{:}] = deal(temp);
        end
    end
    Opt.size = sz;
end

% Check that row/column labels match

if ~isempty(Opt.rowlabel) && length(Opt.rowlabel)~=Opt.size(1)
    error('rowlabel length must match number of rows');
end
if ~isempty(Opt.collabel) && length(Opt.collabel)~=Opt.size(2)
    error('collabel length must match number of columns');
end
    
% Check for stagger flags

if Opt.staggery > 0 && Opt.staggerx > 0
    error('Cannot stagger both sets of axes');
end

% Check for output cell

nout = length(Opt.outputs);

if ~setup
    isemp = cell(nin,1);
    for ii = 1:nin
        isemp{ii} = cellfun('isempty', funin{ii});
    end
    isemp = cat(3, isemp{:});
    isemp = any(isemp, 3);
end

%--------------------------
% Create axes
%--------------------------

h.fig = figure(Opt.figprop{:});

nrow = Opt.size(1);
ncol = Opt.size(2);

axprops = [fieldnames(SubaxOpt) struct2cell(SubaxOpt)]';

for irow = 1:nrow
    for icol = 1:ncol
        h.ax(irow,icol) = subaxis(nrow, ncol, icol, irow, axprops{:}, 'holdaxis');
    end
end

if Opt.staggery > 0
    h.yax = offsetaxis(h.ax(end-1:-2:1,:), 'y', Opt.staggery);
    set(h.ax, 'color', 'none', 'box', 'off');
    set(h.ax(1:end-1,:), 'xcolor', 'none');
end
if Opt.staggerx > 0
    h.xax = offsetaxis(h.ax(:,2:2:end), 'x', Opt.staggerx);
    set(h.ax, 'color', 'none', 'box', 'off');
    set(h.ax(:,2:end), 'ycolor', 'none');
end

%--------------------------
% Plot data
%--------------------------

for iout = 1:nout
    h.(Opt.outputs{iout}) = cell(nrow, ncol);
end

if ~setup
    for irow = 1:nrow
        for icol = 1:ncol

            if ~isemp(irow,icol)
             
                set(h.fig, 'currentaxes', h.ax(irow,icol));
                set(h.ax(irow,icol), 'nextplot', 'add');
                
                invar = cellfun(@(x) x{irow,icol}, funin, 'uni', 0);
                
                tmp = cell(1,nout);
                
                if nout > 0
                    [tmp{:}] = fun(invar{:});
                    for iout = 1:nout
                        h.(Opt.outputs{iout}){irow,icol} = tmp{iout};
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

pos = reshape(get(h.ax, 'Position'), size(h.ax));
xmid = cellfun(@(x) x(1)+x(3)/2, pos(1,:));
ymid = cellfun(@(x) x(2)+x(4)/2, pos(:,1));
xleft = pos{1,1}(1) - (pos{1,end}(1)+pos{1,end}(3) - pos{1,1}(1))*Opt.rowlabeloffset;
ytop = pos{1,1}(4)+pos{1,1}(2) + (pos{1,1}(2)+pos{1,1}(4) - pos{end,1}(2))*Opt.collabeloffset;


labelax = axes('position', [0 0 1 1], 'visible', 'off');
set(labelax, 'xlim', [0 1], 'ylim', [0 1], 'handlevisibility', 'off');

if ~isempty(Opt.rowlabel)
    h.rlab = text(xleft*ones(Opt.size(1),1), ymid, Opt.rowlabel, 'parent', labelax, 'horiz', 'center', 'rotation', 90);
end
   
if ~isempty(Opt.collabel)
    h.clab = text(xmid, ytop*ones(1,Opt.size(2)), Opt.collabel, 'parent', labelax, 'horiz', 'center');
end


    
