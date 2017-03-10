function h = plotgrid(varargin) %x, y, rowname, colname, varargin)
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


    
