function mask = connectdots(xp,yp,H,W)
%Fills in a mask delineated by points
%
%function mask = connectdots(xp,yp,H,W)
%
% This function takes "sparse" points xp,yp (usually on the perimeter of an
% embryo or of an object) and fills in the holes to make a solid mask.
%
% Inputs:
% "xp,yp": the coordinate points around the perimeter of some object that
%	you would like to make a mask for.  These must be the same length and
%	should be simply-connected.  Typically the final point is the same as
%	the first point to make a connected polygon, but this isn't necessary
%	as the program will check for that.  They must be vectors and not 2D or
%	higher-D arrays.
% "H,W": the height and width of your original image.
%
% Outputs:
% "mask": a logical image of dimension H-by-W that is true inside the
%	polygon created by (xp,yp), and false outside.

%
% Size-checking
%
sizex = size(xp);
sizey = size(yp);
if ~isequal(sizex,sizey)
	error('xp and yp must be vectors of the same length')
end
if length(sizex) ~= 2 || length(sizey) ~= 2
	error('xp and yp must be vectors and not arrays')
end
[~,k] = min(sizex);
if sizex(k) ~= 1 
	error('xp and yp must be vectors and not arrays')
end

%
% Ensure polygon wrap-around
%
if xp(1) ~= xp(end) || yp(1) ~= yp(end)
	xp(end+1) = xp(1);
	yp(end+1) = yp(1);
end
nt = length(xp) - 1;

ds = sqrt(diff(xp).^2 + diff(yp).^2);
s = round(sum(ds)); % estimate of the number of points to expect

%
% Loop to fill-in points
%
X = zeros(2*s,1);
Y = X;
count = 0;
for i = 1:nt
	x1 = xp(i); x2 = xp(i+1);
	y1 = yp(i); y2 = yp(i+1);
	
	xm1 = (x1:sign(x2-x1):x2)';
	ym1 = (y2 - y1)/(x2 - x1)*(xm1 - x1) + y1;
	
	ym2 = (y1:sign(y2-y1):y2)';
	xm2 = (x2 - x1)/(y2 - y1)*(ym2 - y1) + x1;
	
	xm = [xm1;xm2];
	ym = [ym1;ym2];
	m = length(xm);
	
	X(count+1:count+m) = xm;
	Y(count+1:count+m) = ym;
	count = count + m;
end
X(count+1:end) = [];
Y(count+1:end) = [];

x = (1:W)';  y = (1:H)';
X1 = roundx(X,x); Y1 = roundx(Y,y);
E = sparse(Y1,X1,1,H,W);
E = full(logical(E));
mask = imfill(E,'holes');

