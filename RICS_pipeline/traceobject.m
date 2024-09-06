function [B,X,Y,Idx] = traceobject(mask,m,n)
%Finds boundary of bw object.
%
%function [B,X,Y,Idx] = traceobject(mask)
%
% This function takes either a bw image "mask", which has exactly one
% object in it, or an n-by-1 structure output of regionprops with field
% "PixelIDxList" (where "n" is the number of objects in the original image
% that regionprops was performed on) and finds the pixels that perfectly
% trace the boundary.  This function is based on the image processing
% toolbox function "bwtraceboundary", but it makes the implementation
% easier for my purposes.
%
% The input "mask" can be one of two things:
% (1) A bw image with exactly one object. In this case, the interpretation
%	of this function is straightforward: the output is a set of pixels that
%	trace that object.  The inputs "m" and "n" are not used.  In this case,
%	the outputs B,X,Y are double arrays, with X being the x-coordinates of
%	the boundary, and Y the y-coords, and B = [Y X].  (Thus, "B" is more in
%	the form of array format, with rows coming first, then columns).  The
%	fourth output, "Idx" is the indices in B but in long-column format.
% (2) A structure which is the output of "regionprops" performed on a bw or
%	label image with "N" objects.  This structure must have the field
%	"PixelIdxList".  Also, you must pass the size of the original image in
%	"m" and "n" (with "m" the number of rows, "n" the number of columns).
%	Alternatively, you could put [nrows ncols] into the input "m", in which
%	case the third input is not needed.  In this case, the outputs change.
%	B becomes the same as the structure array input with an additional
%	field: PerimXYList (see below). X,Y,Idx become cell arrays with N rows
%	and one column, and the i-th element of these corresponds to the
%	"X","Y","Idx" of the i-th object. 
%
% Outputs:
% "B": if the first input is a logical image, then B is an n-by-2 array,
%	where "n" is the number of pixels in the boundary. The first column
%	represents the Y values of the pixels and the second column the X
%	values. If the first input is a structure, then "B" is the same
%	structure, but with the "PerimXYList" field added, and the value of
%	that field is the fliplr of what "B" would have been if the first input
%	were the logical image.
% "X,Y": the "X" and "Y" values described within "B" above
% "Idx": the indices of the "X" and "Y" values in long column format

if isstruct(mask)
	N = length(mask);
	if ~isfield(mask,'PixelIdxList')
		error('You need to have a PixelIdxList in your structure')
	end
	if ~exist('m','var') || (length(m) < 2 && ~exist('n','var'))
		error('If you pass a regionprops struct as 1st input, you also need array size.')
	end
	if length(m) == 2
		n = m(2);
		m = m(1);
	end
else
	N = 1;
end

B = cell(N,1); X = B; Y = B; Idx = B;
for i = 1:N
	
	%
	% Create binary image of just the object
	%
	if isstruct(mask)
		P = mask(i).PixelIdxList;
		bw = false(m,n);
		bw(P) = true;
	else
		[m,n] = size(mask);
		bw = ~~mask;
	end
	
	
	%
	% Need a starting point that is on the boundary of the object. The
	% below lines robustly find a starting point by finding the upper-most
	% point in the left-most column that has any "true" pixels.
	%
	maxbw = max(bw);
	j0 = find(maxbw); j0 = j0(1);
	bwj0 = bw(:,j0); i0 = find(bwj0); i0 = i0(1); % now (i0,j0) is true.
		
	B{i} = bwtraceboundary(bw,[i0,j0],'N');
	x = B{i}(:,2); 
	y = B{i}(:,1);
	X{i} = x; Y{i} = y;
	Idx{i} = Y{i} + m*(X{i} - 1);
	
	%
	% Check to see if there are any holes. To do that, we will calculate
	% the approximate number of pixels contained in the perimeter and
	% compare that to the  number of  pixels in the object. If that differs
	% by more than some small percent (like, 5%?) then we will assume there
	% is a hole in the image. 
	%
	a1 = 0.5*(-sum(y(1:end-1).*diff(x)) + sum(x(1:end-1).*diff(y)) + length(x));
	a2 = sum(bw(:));
	e = (a1 - a2)/a1;
	if e > 0.05
		
		%
		% If there are holes, we will find them as objects using
		% regionprops and then recursively call traceobject 
		%
		holes = imfill(bw,'holes') & ~bw;
		holesstats = regionprops(holes,'Area','PixelIdxList');
		b = traceobject(holesstats,m,n);
		
		xh = []; yh = [];
		for j = 1:length(b)
			xh = [xh;b(j).PerimXYList(:,1);NaN];
			yh = [yh;b(j).PerimXYList(:,2);NaN];
		end
		xh(end) = []; yh(end) = [];
		X{i} = [x;NaN;xh]; Y{i} = [y;NaN;yh];
		B{i} = [Y{i} X{i}];
		Idx{i} = Y{i} + m*(X{i} - 1);
	end
	
	
end

if N == 1 && ~isstruct(mask)
	B = B{1}; X = X{1}; Y = Y{1}; Idx = Idx{1};
elseif isstruct(mask)
	
	for j = 1:length(Idx)
		mask(j).PerimXYList = fliplr(B{j});
	end	
	B = mask;
end




