function C = jac_dydtheta(t,dYdtheta,params,Vn,Vc,An,Am,stage,deval_soln)
% This function calculates the Jacobian of the functions from 
% ftn_dsat with respect to the State Variables

%load('sample2.mat')

Y       = deval(deval_soln,t);
M       = length(Y)/6;
nparams = length(params);
N       = length(dYdtheta);
nblock  = N/nparams;


%
% Unpack state variables
%
un   = Y(1:M);
uc = Y(M+1:2*M);
% wn = Y(2*M+1:3*M);
wc = Y(3*M+1:4*M);
vn = Y(4*M+1:5*M);
vc = Y(5*M+1:6*M);


%
% Unpack params
%
lambdaU = params(1);
lambdaW = params(2);
lambdaV = params(3);
sigmaU  = params(4);
sigmaW  = params(5);
sigmaV  = params(6);
muU     = params(7);
muW     = params(8);
muV     = params(9);
gamma   = params(10);
psi     = params(11);
alpha   = params(12);
beta    = params(13);
phi     = params(14);
beta0   = params(15);
kappa   = params(16);
delta   = params(17);
K       = params(18);


% Groupings
if strcmp(stage,'interphase')
	an = An/Vn;
elseif strcmp(stage,'mitosis')
	an = 0;
end


% Other matrices
I      = ones(M,1);
Z      = zeros(M,1);
Pmdiag = -2*I;
% Prdiag = [0;2;I(1:M-2)];
% Pldiag = [I(1:M-2);2;0];

Prdiag = [2;I(1:M-2)];
Pldiag = [I(1:M-2);2];

%
% Defining other derivatives
%
% 1. Derivative of the Michaelis-Menten function
% ftn = g*wc./(kappa + wc)
x = linspace(0,1,M)';
g = exp(-x.^2/2/phi^2);
g = g.*kappa./(kappa + wc).^2;


% 2. Derivative of the function that produces Cact neg Fbk
% f = un.^nH./(K.^nH + un.^nH);
nH = 20;
f = nH*un.^(nH-1)*K.^nH./(K.^nH + un.^nH).^2;




%
% Jacobian elements
% 
% 1st row block, derivatives of f1, the eqn for un.
% f1 = (sigmaU*uc - muU*un) - gamma*un.*vn + beta0*wn;
f1_un = -muU*an*I - gamma*vn;
f1_uc = sigmaU*an*I;
f1_wn = beta0*I;
f1_wc = Z;
f1_vn = -gamma*un;
f1_vc = Z;


% Second row block, derivatives of f2, the eqn for uc.
% f2 = lambdaU*P*uc - (gamma*uc.*vc - beta*g*wc./(kappa + wc) - beta0*wc) ...
%	- (sigmaU*uc - muU*un);
f2_un = muU*An/Vc*I;
f2_uc = lambdaU*Am/Vc*Pmdiag - gamma*vc - sigmaU*An/Vc*I;
f2_wn = Z;
f2_wc = beta*g + beta0*I;
f2_vn = Z;
f2_vc = -gamma*uc;

f2_uc_l = lambdaU*Am/Vc*Pldiag;
f2_uc_r = lambdaU*Am/Vc*Prdiag;


% Third row block, derivatives of f3, the eqn for wn.
% f3 = (sigmaW*wc - muW*wn) + gamma*un.*vn - beta0*wn;
f3_un = gamma*vn;
f3_uc = Z;
f3_wn = -muW*an*I - beta0*I;
f3_wc = sigmaW*an*I;
f3_vn = gamma*un;
f3_vc = Z;


% Fourth row block, derivatives of f4, the eqn for wc.
% f4 = lambdaW*P*wc + (gamma*uc.*vc - beta*g*wc./(kappa + wc) - beta0*wc)...
%	- (sigmaW*wc - muW*wn);
f4_un = Z;
f4_uc = gamma*vc;
f4_wn = muW*An/Vc*I;
f4_wc = lambdaW*Am/Vc*Pmdiag - sigmaW*An/Vc*I - beta*g - beta0*I;
f4_vn = Z;
f4_vc = gamma*uc;

f4_wc_l = lambdaW*Am/Vc*Pldiag;
f4_wc_r = lambdaW*Am/Vc*Prdiag;


% Fifth row block, derivatives of f5, the eqn for zn.
% f5 = (sigmaV*vc - muV*vn) - psi*(gamma*un.*zn - beta0*wn);
f5_un = -psi*gamma*vn;
f5_uc = Z;
f5_wn = psi*beta0*I;
f5_wc = Z;
f5_vn = -muV*an*I - psi*gamma*un;
f5_vc = sigmaV*an*I;


% Sixth row block, derivatives of f6, the eqn for zc.
% f6 = lambdaV*P*vc - psi*(gamma*uc.*vc - beta*g*wc./(kappa + wc)) + ...
%	1 + delta*f - alpha*vc - (sigmaV*vc - muV*vn);
f6_un = delta*f;
f6_uc = -psi*gamma*vc;
f6_wn = Z;
f6_wc = psi*(beta*g + beta0*I);
f6_vn = muV*An/Vc*I;
f6_vc = lambdaV*Am/Vc*Pmdiag - psi*gamma*uc - alpha*I - sigmaV*An/Vc*I;

f6_vc_l = lambdaV*Am/Vc*Pldiag;
f6_vc_r = lambdaV*Am/Vc*Prdiag;


fdiag = [[f1_un; f1_uc; f1_wn; f1_wc; f1_vn; f1_vc], ...
         [f2_un; f2_uc; f2_wn; f2_wc; f2_vn; f2_vc], ...
         [f3_un; f3_uc; f3_wn; f3_wc; f3_vn; f3_vc], ...
         [f4_un; f4_uc; f4_wn; f4_wc; f4_vn; f4_vc], ...
         [f5_un; f5_uc; f5_wn; f5_wc; f5_vn; f5_vc], ...
         [f6_un; f6_uc; f6_wn; f6_wc; f6_vn; f6_vc]];

     
m  = M;
ii = [1:m, 1:m-1, 2:m];     %rows:    main diagonal, right diagonal, left diagonal
jj = [1:m, 2:m, 1:m-1];     %columns: main diagonal, right diagonal, left diagonal
vv = [-2*ones(1,m), 2, ones(1,2*m-4), 2]; 
P  = sparse(ii,jj,vv);

P2 = full(P);
     
%
% unpack lambda
%
C = zeros(N,1);
for i=1:nparams
   iblock       = (i-1)*nblock+1:i*nblock;
   delta      = dYdtheta(iblock);
   
   del1 = delta(1:M);
   del2 = delta(M+1:2*M);
   del3 = delta(2*M+1:3*M);
   del4 = delta(3*M+1:4*M);
   del5 = delta(4*M+1:5*M);
   del6 = delta(5*M+1:6*M);
   
   
   C1 = f1_un.*del1 + f1_uc.*del2 + f1_wn.*del3 + f1_wc.*del4 + f1_vn.*del5 + f1_vc.*del6;
   C2 = f2_un.*del1 + f2_uc.*del2 + f2_wn.*del3 + f2_wc.*del4 + f2_vn.*del5 + f2_vc.*del6; 
   C3 = f3_un.*del1 + f3_uc.*del2 + f3_wn.*del3 + f3_wc.*del4 + f3_vn.*del5 + f3_vc.*del6; 
   C4 = f4_un.*del1 + f4_uc.*del2 + f4_wn.*del3 + f4_wc.*del4 + f4_vn.*del5 + f4_vc.*del6;   
   C5 = f5_un.*del1 + f5_uc.*del2 + f5_wn.*del3 + f5_wc.*del4 + f5_vn.*del5 + f5_vc.*del6;
   C6 = f6_un.*del1 + f6_uc.*del2 + f6_wn.*del3 + f6_wc.*del4 + f6_vn.*del5 + f6_vc.*del6;

   
   
   % first add diagonal elements
   %C(iblock,1)  = sum(fdiag.*delta,2);
   
   % then add the off-diagonal elements from the sub-matrices
   mblock    = M+1:2*M;
   deltr     = delta(mblock(2:end));
   deltl     = delta(mblock(1:end-1));
   C2 = C2 + [f2_uc_r.*deltr;0] + [0;f2_uc_l.*deltl];
   %C2        = C2+ [0;f2_uc_r.*deltr] + [f2_uc_l.*deltl;0];
   
   mblock    = 3*M+1:4*M;
   deltr     = delta(mblock(2:end));
   deltl     = delta(mblock(1:end-1));
   %C4        = C4 + [f4_wc_r.*deltr;0] + [0;f4_wc_l.*deltl];
   C4        = C4 + [f4_wc_r.*deltr;0] + [0;f4_wc_l.*deltl];
   
   mblock    = 5*M+1:6*M;
   deltr     = delta(mblock(2:end));
   deltl     = delta(mblock(1:end-1));
   C6        = C6+ [f6_vc_r.*deltr;0] + [0;f6_vc_l.*deltl];
   
   C(iblock,1) = [C1; C2; C3; C4; C5; C6]';
end

end



