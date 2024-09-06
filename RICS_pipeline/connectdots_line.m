function [mask,xp,yp] = connectdots_line(xp,yp,H,W)
%Fills in a mask delineated by points; extrapolates ends to edge of image
%
%function mask = connectdots(xp,yp,H,W)
%
% This function takes "sparse" points xp,yp (usually on the border of an
% embryo or of an object) and connects the dots to make a solid mask. The
% points are assumed to divide the image into two separate pieces. The
% endpoints are not assumed to be on the edge of the image, so the curve
% dividing the image in two is extrapolated using a line from the last two
% endpoints on both sides.
%
% Inputs:
% "xp,yp": the coordinate points on the border of some object that divides
%	the image in two parts. These must be the same length and should be
%	simply-connected.  They must be vectors and not 2D or higher-D arrays.
% "H,W": the height and width of your original image.
%
% Outputs:
% "mask": a logical image of dimension H-by-W. The "true" pixels are
%	determined by the curve (xp,yp), which divides the image into two
%	parts. The part of the image that has the smaller area is considered the
%	mask (i.e., pixels are set to true). The true pixels include the curve
%	(xp,yp).
% "xp,yp": The original inputs with the endpoints extrapolated and the
%	corner(s) added.

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
% Extrapolate endpoints to the edge of the image. To do this, we first must
% figure out when edge of the image the line will hit. We start with the
% first two points.
%
dx = xp(2) - xp(1); dy = yp(2) - yp(1);
if dx > 0
	nx = xp(1) - 1; 
	xextrap = 1;
elseif dx < 0
	nx = W - xp(1);
	xextrap = W;
else
	nx = Inf;
	xextrap = Inf;
end
if dy > 0
	ny = yp(1) - 1; 
	yextrap = 1;
elseif dy < 0
	ny = H - yp(1); 
	yextrap = H;
else
	ny = Inf;
	yextrap = Inf;
end
if nx*abs(dy) < ny*abs(dx)
	xp0 = xextrap;
	yp0 = yp(1) - nx*dy/abs(dx);
else
	xp0 = xp(1) - ny*dx/abs(dy);
	yp0 = yextrap;
end

%
% Now do the same thing for the final two points.
%
dx = xp(end) - xp(end-1); dy = yp(end) - yp(end-1);
if dx < 0
	nx = xp(end) - 1; 
	xextrap = 1;
elseif dx > 0
	nx = W - xp(end);
	xextrap = W;
else
	nx = Inf;
	xextrap = Inf;
end
if dy < 0
	ny = yp(end) - 1; 
	yextrap = 1;
elseif dy > 0
	ny = H - yp(end); 
	yextrap = H;
else
	ny = Inf;
	yextrap = Inf;
end
if nx*abs(dy) < ny*abs(dx)
	xpend = xextrap;
	ypend = yp(end) + nx*dy/abs(dx);
else
	xpend = xp(end) + ny*dx/abs(dy);
	ypend = yextrap;
end

xp = round([xp0;xp;xpend]);
yp = round([yp0;yp;ypend]);

%
% Now we look at different cases to figure out which corner points to put
% in. The corner points are so that the (xp,yp) vectors that define the
% boundary make a closed polygon.
%
if (xp(1) == 1 && yp(1) == 1) || (xp(1) == 1 && yp(1) == H) || ...
		(xp(1) == W && yp(1) == 1) || (xp(1) == W && yp(1) == H) || ...
		(xp(end) == 1 && yp(end) == 1) || (xp(end) == 1 && yp(end) == H) || ...
		(xp(end) == W && yp(end) == 1) || (xp(end) == W && yp(end) == H)
		% In this case, the boundary contains a point that is a corner of
		% the image. When this file was originally written, this case was
		% not worked out, so it simply returns a NaN mask
	mask = NaN;
	return
	
elseif (xp(1) == 1 || xp(end) == 1) && (xp(1) == W || xp(end) == W)
	A1 = abs(trapz(xp,yp-1));
	A2 = abs(trapz(xp,H-yp));
		
	xp = [xp(1); xp; xp(end); xp(1)];
	if A1 < A2
		yp = [1; yp; 1; 1];
	else
		yp = [H; yp; H; H];
	end
	
elseif (yp(1) == 1 || yp(end) == 1) && (yp(1) == H || yp(end) == H)
	
	A1 = abs(trapz(yp,xp-1));
	A2 = abs(trapz(yp,W-xp));
		
	yp = [yp(1); yp; yp(end); yp(1)];
	if A1 < A2
		xp = [1; xp; 1; 1];
	else
		xp = [W; xp; W; W];
	end
	
elseif xp(1) == 1 || xp(end) == 1
	xp = [1; xp; 1];
	
	if yp(1) == 1 || yp(end) == 1
		yp = [1; yp; 1];
	elseif yp(1) == H || yp(end) == H
		yp = [H; yp; H;];
	end
		
elseif xp(1) == W || xp(end) == W
	xp = [W; xp; W];
	
	if yp(1) == 1 || yp(end) == 1
		yp = [1; yp; 1];
	elseif yp(1) == H || yp(end) == H
		yp = [H; yp; H];
	end
	
end

%
% Now that we have completed our xp,yp set with the corner points, we have
% an enclosed region, so we can call on "connectdots."
%
mask = connectdots(xp,yp,H,W);


