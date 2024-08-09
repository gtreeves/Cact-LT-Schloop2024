function [H,fdot0,f] = calc_hessian_grads(x,gradhandle)
  
if ~exist('gradhandle','var')
   gradhandle = @dsat2_gradcalc; 
end


np   = length(x);   
pert = 1e-2;

xlin            = x;
X               = repmat(xlin,np,1);
X               = X + diag(pert*(xlin)); 
X               = [X;xlin];
[F,~,~,Soln]    = calc_error(X);
Fdot            = gradhandle(Soln);
fdot            = Fdot(1:np,1:np);
fdot0           = Fdot(np+1,1:np);
H               = (fdot - fdot0)./(pert*(xlin'));
H               = 0.5*(H + H');
f               = F(end);

% LOOP
%{
for i=1:nParams
   x_pert         = x;
   x_pert(i)      = x_pert(i) + pert*(x_pert(i));    
   [~,~,~,soln]   = calc_error(log10(x_pert));
   fdot           = gradhandle(soln);
   fdot           = fdot(1:16);

   H(i,:)         = (fdot - fdot0)/(pert*(x_pert(i)));
end
toc
H = 0.5*(H + H');
%}
end  