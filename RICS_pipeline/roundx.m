function [Y,idx] = roundx(y,x)
%Rounds element of vector "y" to closest element found in vector "x".
%
%function [Y,idx] = roundx(y,x)
%
% This function takes each element of vector "y" and finds the closest
% value to that element within the vector "x", outputting the corresponding
% element in vector "Y".  For example, if y(i) is closer to x(j) than to
% x(k) for any k ~= j, then Y(i) = x(j).
%
% "y": vector that you wish to round.  Can be scalar, vector, or 2D array
% "x": vector containing values you wish to round to.  Must be scalar or
%	vector.
%
% "Y": output, same size as y, and Y(i) is "close to" y(i).
% "idx": indices that we use to round: Y = x(idx).  Same size as "y".

[my,ny] = size(y);
[mx,nx] = size(x);

if (nx > 1 && mx > 1)
	error('Inputs "x" must be a vector or scalar')
end

%
% First, we make "y" into a column vector, and "x" into a row vector
%
y = y(:);
x = x(:)';

Ly = length(y);
Lx = length(x);

Y = repmat(y,1,Lx);
X = repmat(x,Ly,1);

D = abs(Y - X);
[~,idx] = min(D,[],2);
Y = x(idx)'; % since we made x into a row vec, trnsps so that Y a col vec.

%
% Now we transform "Y" and "idx" into the original size of "y".
%
count1 = (1:my)'; A = repmat(count1,1,ny);
count2 = 0:(ny-1); B = repmat(count2,my,1);
t = A + my*B;
Y = Y(t);
idx = idx(t);





