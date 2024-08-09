function grad = calc_grad_finitedifference(xb,F,pert)

if isstruct(xb)
   xb = log10(xb.params); 
end
if ~exist('F','var')
    F = calc_error(xb);
end
if ~exist('pert','var')
    pert = 1e-7;        % perturbation
end


%
% Loop over all params
%
nParams = length(xb);
grad    = zeros(length(xb),1);

for i=1:nParams
    x_pert     = 10.^xb;
    x_pert(i)  = x_pert(i) + pert*(x_pert(i));    
    Fpert      = calc_error(log10(x_pert));
    grad(i)    = (Fpert - F)/(pert*x_pert(i));
end

end