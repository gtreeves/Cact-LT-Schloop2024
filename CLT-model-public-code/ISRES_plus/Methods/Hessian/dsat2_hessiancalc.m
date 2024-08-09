function H = dsat2_hessiancalc(params,idx,soln0,fdot0)
%
% This function will calculate the Hessian by doing a for-loop.
%
% "params": the parameter set in log-space.  Passed to dsat_errorcalc,
%	which expects them to be in log-space.
%
% "H": the Hessian, which is in log-space.


if ~exist('idx','var') || isempty(idx)
	nparams = length(params);   
	idx = 1:nparams;
else
	nparams = length(idx);
end


if ~exist('soln0','var') 
	%[~,~,soln0] = dsat2_errorcalc_fake(params,false);
    [~,~,soln0] = calc_error(params);
end


if ~exist('fdot0','var')
	fdot0  = dsat2_gradcalc(soln0);
    fdot0  = fdot0(1:16);
end



delt = 1e-5;
H    = zeros(nparams);
for i = 1:nparams
	params1         = params; 
	params1(idx(i)) = params1(idx(i)) + log10(1 + delt);
	
	[~,~,solnwt1]   = calc_error(params1);
	fdot1           = dsat2_gradcalc(solnwt1);
	fdot1           = fdot1(1:16);
    
	H(i,:)          = (fdot1 - fdot0)/(delt*10^params(i)*log(10));
    %H1(i,:) = (fdot1 - fdot0)/log10(1 + delt);
    
    %H(i,:) = (fdot1 - fdot0)/(delt).*params(idx);
end

H = 0.5*(H + H');








