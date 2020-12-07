function h = pcolor3(varargin)
%PCOLOR Pseudocolor (checkerboard) plot.
%   PCOLOR(C) is a pseudocolor or "checkerboard" plot of matrix C.
%   The values of the elements of C specify the color in each
%   cell of the plot. In the default shading mode, 'faceted',
%   each cell has a constant color and the last row and column of
%   C are not used. With shading('interp'), each cell has color
%   resulting from bilinear interpolation of the color at its 
%   four vertices and all elements of C are used. 
%   The smallest and largest elements of C are assigned the first and
%   last colors given in the color table; colors for the remainder of the 
%   elements in C are determined by table-lookup within the remainder of 
%   the color table.
%
%   PCOLOR(X,Y,C), where X and Y are vectors or matrices, makes a
%   pseudocolor plot on the grid defined by X and Y.  X and Y could 
%   define the grid for a "disk", for example.
%
%   PCOLOR(AX,..) plots into AX instead of GCA.
%
%   H = PCOLOR(...) returns a handle to a SURFACE object.
%
%   PCOLOR is really a SURF with its view set to directly above.
%
%   See also CAXIS, SURF, MESH, IMAGE, SHADING.

%-------------------------------
%   Additional details:
%
%
%   PCOLOR sets the View property of the SURFACE object to directly 
%   overhead.
%
%   If the NextPlot axis property is REPLACE (HOLD is off), PCOLOR resets 
%   all axis properties, except Position, to their default values
%   and deletes all axis children (line, patch, surf, image, and 
%   text objects).  View is set to [0 90].

%   Copyright 1984-2015 MathWorks, Inc. 

%   J.N. Little 1-5-92

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 4
    error(message('MATLAB:narginchk:tooManyInputs'));
end
% do error checking before calling newplot. This argument checking should
% match the surface(x,y,z) or surface(z) argument checking.
if nargs == 2
  error(message('MATLAB:pcolor:InvalidNumberOfInputs'))
end
if isvector(args{end})
  error(message('MATLAB:pcolor:NonMatrixColorInput'));
end
if nargs == 3 && LdimMismatch(args{1:3})
  error(message('MATLAB:pcolor:InputSizeMismatch'));
end
for k = 1:nargs
  if ~isreal(args{k})
    error(message('MATLAB:pcolor:NonRealInputs'));
  end
end

cax = newplot(cax);
hold_state = ishold(cax);

if nargs == 1
    x = args{1};
    hh = surface(zeros(size(x)),x,'Parent',cax);
    [m,n] = size(x);
    lims = [ 1 n 1 m];
elseif nargs == 3
    [x,y,c] = deal(args{1:3});
    hh = surface(x,y,zeros(size(c)),c,'Parent',cax);
    lims = [min(min(x)) max(max(x)) min(min(y)) max(max(y))];
elseif nargs==4
    [x,y,z,c] = deal(args{1:4});
    hh = surface(x,y,z,c,'Parent',cax);
    lims = [min(min(x)) max(max(x)) min(min(y)) max(max(y)) min(min(z)) max(max(z))];
end
set(hh,'AlignVertexCenters','on');
if ~hold_state
    set(cax,'View',[0 90]);
    set(cax,'Box','on');
    if lims(2) <= lims(1)
        lims(2) = lims(1)+1;
    end
    if lims(4) <= lims(3)
        lims(4) = lims(3)+1;
    end
    if lims(6) <= lims(5)
        lims(6) = lims(5)+1;
    end
    axis(cax,lims);
end
if nargout == 1
    h = hh;
end

function ok = LdimMismatch(x,y,z)
[xm,xn] = size(x);
[ym,yn] = size(y);
[zm,zn] = size(z);
ok = (xm == 1 && xn ~= zn) || ...
     (xn == 1 && xm ~= zn) || ...
     (xm ~= 1 && xn ~= 1 && (xm ~= zm || xn ~= zn)) || ...
     (ym == 1 && yn ~= zm) || ...
     (yn == 1 && ym ~= zm) || ...
     (ym ~= 1 && yn ~= 1 && (ym ~= zm || yn ~= zn));



