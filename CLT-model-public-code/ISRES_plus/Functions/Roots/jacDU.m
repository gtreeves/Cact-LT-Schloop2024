function f_x = jacDU(fHandle,x0,f,varargin)
%Numerical calculation of the Jacobian.
%
%function f_x = jacDU(fHandle,x0,f,varargin)
%
%'fHandle' is a function handle that points to the 
%   multivariate function of which the jacobian is to be taken.
%'x0' is the point at which the derivative should be taken.
%Optional argument 'varargin' can consist of two things, in this order:
%   (1) 'dxj': The finite difference size.  Default is 0.01.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.
%   (2) 'p': Other parameters that may change the evaluation of the 
%       function, but the function is not actually differentiated wrt them.

N = size(x0,1);
nArg = size(varargin,2);
% f_x = [];

if nArg > 0 && isnumeric(varargin{1})
    dxj = varargin{1};
else
    dxj = 0.01;
end

if nArg > 1
	p = {varargin{2:nArg}};
% 	existOtherParameters = 1;
else
	p = {};
%     existOtherParameters = 0;
end

f_x = zeros(N,N);
% if existOtherParameters
if dxj < 0
	for j = 1:N
		dx = zeros(N,1);
		dx(j) = dxj;
		f_x(:,j) = (fHandle(x0*(1+dx),p{:}) - f)/dxj/x0(j);
	end
else
	for j = 1:N
		dx = zeros(N,1);
		dx(j) = dxj;
		f_x(:,j) = (fHandle(x0+dx,p{:}) - f)/dxj;
	end
end
% else
%     for j = 1:N
% 		dx = zeros(N,1);
% 		dx(j) = dxj;
%         f_x(:,j) = (fHandle(x0+dx) - f)/dxj;
%     end
% end

% if existOtherParameters
%     for j = 1:N
% 		dx = zeros(N,1);
% 		dx(j) = dxj;
%         f_x = [f_x (feval(fHandle,x0 + dx,p{:}) - feval(fHandle,x0 - dx,p{:}))/(2*dx(j))];
%     end
% else
%     for j = 1:N
% 		dx = zeros(N,1);
% 		dx(j) = dxj;
%         f_x = [f_x (feval(fHandle,x0 + dx) - feval(fHandle,x0 - dx))/(2*dx(j))];
%     end
% end