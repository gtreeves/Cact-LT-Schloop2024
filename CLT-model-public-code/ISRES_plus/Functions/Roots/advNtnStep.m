function [x,x1,x2,f1,f2,deltax,stepflag] = advNtnStep(x0,x1,x2,f,f1,f2,deltax)
%Advanced newton on monotone, scalar ftn, with no dir and perhaps unbdd
%
%function [x,x1,x2,f1,f2,deltax] = advNtnStep(x0,x1,x2,f,f1,f2,deltax)
%
% This function works only on helping to take the newton step on scalar
% nonlinear equations.  It needs as input, in particular, your starting
% guess (x0), the value of the function (for which you are trying to find
% the root) at x0 (f), and the newton step deltax = -J \ f.  It is mostly
% useful when you have natural bounds on your function, x1 and x2 (it
% doesn't matter whether x1 < x2 or not).  If you have (either of) these
% bounds, you also should pass f1 and f2, the values of f at the respective
% x's.  Once given these parameters, this function will give you a smart
% estimate of your next x (and in addition update your upper or lower
% bounds, and return the deltax that the function actually took).  For
% example, if your original deltax is such that x = x0 + deltax would be
% past one of your natural boundaries, this function will instead modify
% deltax so that it will not jump past your boundary.

if isempty(f1) && isempty(f2)
	x1 = x0; f1 = f;
	stepflag = 1;
	
elseif isempty(f1)
	if f*f2 < 0 % means root is in [x0,x2]
		[x1,f1,deltax] = rootByX2(x0,x2,f,f2,deltax);
		stepflag = 2;
	elseif f*f2 > 0 % means root is in [x1,x0]
		stepflag = 3;
		if deltax*(x2 - x0) > 0 % but if step is contrary
			deltax = -deltax; % switch directions and hope for the best.
			stepflag = 4;
		end
		x2 = x0; f2 = f;
	else % means root is on x0
		 % do nothing.
		 stepflag = 5;
	end

elseif isempty(f2)
	if f*f1 > 0 % means root is in [x0,x2]
		stepflag = 6;
		if deltax*(x1 - x0) > 0 % but if step is contrary
			deltax = -deltax; % switch directions and hope for the best.
			stepflag = 7;
		end
		x1 = x0; f1 = f;
	elseif f*f1 < 0 % means root is in [x1,x0]
		[x2,f2,deltax] = rootByX1(x0,x1,f,f1,deltax);
		stepflag = 8;
	else % means root is on x0
		 % do nothing.
		 stepflag = 9;
	end
	

	
elseif f1*f2 < 0 % means we can use our bounds
	if f*f1 > 0 % means root is in [x0,x2]
		[x1,f1,deltax] = rootByX2(x0,x2,f,f2,deltax);
		stepflag = 10;
	elseif f*f1 < 0% means root is in [x1,x0]
		[x2,f2,deltax] = rootByX1(x0,x1,f,f1,deltax);
		stepflag = 11;
	else % means root is on x0
		 % do nothing.
		 stepflag = 12;
	end
	
else % means we cannot use our bounds
	% do nothing
	stepflag = 13;
end

x = x0 + deltax; % now we have a new iterate


%----------- subfunction for if root is in [x0,x2] ---------------
function [x1,f1,deltax] = rootByX2(x0,x2,f,f2,deltax)

x1 = x0; f1 = f; % define new lower bounds
if deltax*(x2 - x0) < 0 % but if step is contrary
	deltax = -f/(f2 - f)*(x2 - x0); % redefine as secant step
elseif abs(deltax) >= abs(x2 - x0); % but if step is too big
	deltax = (x2 - x0)*deltax/((x2 - x0) + deltax); % hypblc step
end

%----------- subfunction for if root is in [x1,x0] ---------------
function [x2,f2,deltax] = rootByX1(x0,x1,f,f1,deltax)

x2 = x0; f2 = f; % define new upper bounds
if deltax*(x1 - x0) < 0 % but if step is contrary
	deltax = -f/(f1 - f)*(x1 - x0); % redefine as secant step
elseif abs(deltax) >= abs(x1 - x0); % but if step is too big
	deltax = (x1 - x0)*deltax/((x1 - x0) + deltax); % hypblc step
end






