function [x,J,f,nSteps] = advNtnDU(fhandle,x0,varargin) 
%Newton's method for scalar functions.  Uses bounds.
%
%function [x,J,f] = advNtnDU(fHandle,x0,varargin) 
%
%'fhandle' is a function handle that points to the 
%   multivariate function for which the zero is to be found.
%'x0' is an intial guess.
%Optional argument varargin can consist of six things, in this order:
%   * 'bounds': The bounds on x and the evaluation of the function at those
%		bounds.  This array looks like: [xL xU; fL fU];
%		You may specify the x-bounds without the f-bounds, in which case,
%		this function will evaluate the f-bounds for you.  If the value of
%		f is not real and of different signs at the two bounds, then an
%		error occurs.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.  
%   * 'convergenceCriterion': The convergence tolerance.  Default is 
%       set to 1e-6.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.  
%	* 'nStepsMax': The maximum number of steps before the program quits.
% 		Default is set to 100. 
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument. 
%	* 'jacHandle': the Jacobian matrix handle, if the user supplies
%		one.  It is important that the parameters passed to this
%		function are identical to those passed to the fHandle.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument. 
%	*	'dxj': The finite difference size for the numerical Jacobian evaluation in
%		calling jacDU(fHandle,x0,varargin).  Of course, you would leave this out if
%		you specified a 'jacHandle'.  Default is 1e-4.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.
%	*	'alpha': a relaxation constant, that premultiplies the Jacobian inverse.
%		Must be less than one.  Default is 1.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.
%	*	'p': A list (which becomes a cell array) of any other parameters that
%		may change the evaluation of the function, but are not actually 
%		varied to find the root.
%       
% x is the final answer, and J is the last jacobian.

nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	bounds = varargin{iArg}; else
	bounds = [];
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	convergenceCriterion = varargin{iArg}; else
	convergenceCriterion = 1e-6;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nStepsMax = varargin{iArg}; else
	nStepsMax = 100;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	jacHandle = varargin{iArg}; else
	jacHandle = @jacDU;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	dxj = varargin{iArg}; else
	dxj = 1e-4;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Alpha1 = varargin{iArg}; else
	Alpha1 = 1;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
% 	p = {varargin{iArg:nArg}}; else
	p = varargin(iArg:nArg); else
    p = {};
end%, iArg = iArg + 1;

%
% Evaluating the bounds
%
[m,n] = size(bounds);
if m < 2 && n == 2
	x1 = bounds(1,1);
	x2 = bounds(1,2);
	f1 = fhandle(x1,p{:});
	f2 = fhandle(x2,p{:});
elseif m == 2 && n == 2
	x1 = bounds(1,1);
	x2 = bounds(1,2);
	f1 = bounds(2,1);
	f2 = bounds(2,2);
	
else
	x1 = NaN;
end

%
% Running the loop
%
f = 1; delta_x = 1; x = x0; nSteps = 1;
while norm([f;delta_x]) > convergenceCriterion && nSteps < nStepsMax
	f = fhandle(x,p{:});
	J = jacHandle(fhandle,x,f,dxj,p{:});
	delta_x = - J \ f;
	
	if ~isnan(x1)
		[x,x1,x2,f1,f2,delta_x,stepflag] = advNtnStep(x,x1,x2,f,f1,f2,delta_x); %#ok<ASGLU>
	else
		x = x + Alpha1*delta_x;
	end
	F(nSteps) = f;
	X(nSteps) = x;
	F1(nSteps) = f1;
	F2(nSteps) = f2;
	Deltax(nSteps) = delta_x;
	X1(nSteps) = x1;
	X2(nSteps) = x2;
	Stepflag(nSteps) = stepflag;
	nSteps = nSteps + 1;
end

if delta_x > convergenceCriterion
    x = 'Did not converge';
end