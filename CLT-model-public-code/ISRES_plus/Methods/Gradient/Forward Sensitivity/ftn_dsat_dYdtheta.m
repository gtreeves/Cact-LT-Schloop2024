function dYdt = ftn_dsat_dYdtheta(t,dYdtheta,params,x,P,Vn,Vc,An,Am,stage,deval_soln)

% nparams = length(params);
% N = length(dYdtheta);
% m = size(P,1);
% ns = N/nparams/m;
% dYdtheta = reshape(dYdtheta,[m*ns,nparams]);

% n_params = length(params);
% if strcmp(modelname,'decon')
%     index = 1:(n_params-3);
% elseif strcmp(modelname,'dsat')
%     index = 1:(n_params-2);
% elseif strcmp(modelname,'negfb')
%     index = 1:(n_params);
% end



G     = G_dsat(t,dYdtheta,params,x,P,Vn,Vc,An,Am,stage,deval_soln);

% Multiply Jac and dYdtheta element-wise 
C     = jac_dydtheta(t,dYdtheta,params,Vn,Vc,An,Am,stage,deval_soln);
dYdt  = C + G(:);



% Old method
%{
% %
% % Method 1 - using sparse diagonals
% %
% J     = jac(t,dYdtheta,params,x,P,Vn,Vc,An,Am,stage,deval_soln);
% 
% 
% %
% % Compare Method 1 and 2
% %
% b = J*dYdtheta;
% c = [C,b];
% a = norm(C-b)
% 
% 
% dYdt  = J*dYdtheta + G(:);


%dYdt  = J*dYdtheta + G;
% dYdt = dYdt(:);
%}
