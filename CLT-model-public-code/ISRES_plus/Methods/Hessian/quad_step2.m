function [H,lam,xquad,X0,a] = quad_step2(p,f,xcen,d,beta,lb,ub)
%
%
%
%
% This function calculates the best-fit hyper-paraboloid (dimensionality
% n, where n = number of parameters) using a group of "np" parameter
% sets, their corresponding objective function values, "f", and their
% gradients "fdot".
%
% Actually, this function is "2" because we are not using the gradients
% now.

xquad = [];
X0 = [];
a = [];


[np,n] = size(p);

if np < 0.5*(n^2+n)+n+1
	error('need more param sets')
end

% p = p - repmat(xcen,np,1);

O = true(n); O = triu(O);
Z = zeros(n);

%
% Create matrix that defines the f-portion of the fit
%
P = zeros(np,0.5*(n^2+n));
for i = 1:np
	p1 = p(i,:);
	P1 = p1'*p1;
	P1 = 2*P1 - diag(diag(P1));
	P(i,:) = P1(O)';
end
A = [P p ones(np,1)];

% lastwarn('')
lam = A \ f;
% msg = lastwarn;
% if ~isempty(msg)
% 	1;
% end

a = lam(1:0.5*(n^2+n));
b = lam(0.5*(n^2+n)+1:0.5*(n^2+n)+n);

%
% Convert "a" back into the approximate Hessian: make a symmetric
% n-by-n matrix out of it.
%
H = Z;
H(O) = a;
H = H + H' - diag(diag(H));
lam = eig(H);

% 
% if ~all(lam > 0)
% 	xquad = NaN;
% 	X0 = NaN(n,1);	
% else
% 	X0 = (-(H \ b))';
% 	
% 	if exist('xcen','var')
% 	
% 		
% 		Now the quadratic explanatory model suggests that X0 is your set
% 		of optimum parameter values.  This can be discarded if your
% 		predicted X0 is too far from the centroid of your current x cloud.
% 		
% 		Dsquared = sum((X0 - xcen).^2,2);
% 		D = sqrt(Dsquared/n);
% 		if D < beta*d;
% 			dXm = d.*randn(mu,n); % mutation
% 			xquad = repmat(X01,mu,1) + dXm;
% 			LB = lb(1:mu,:);
% 			V = xquad < LB;
% 			xquad(V) = LB(V);
% 			UB = ub(1:mu,:);
% 			V = xquad > UB;
% 			xquad(V) = UB(V);
% 		end
% 	else
% 		xquad = 'Not asked for';
% 	end
% end
