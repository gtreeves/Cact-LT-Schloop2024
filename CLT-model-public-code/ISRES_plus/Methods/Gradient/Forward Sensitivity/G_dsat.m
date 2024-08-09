function G = G_dsat(t, Y, params, x, P, Vn, Vc, An, Am, stage, deval_soln,index)

% This function calculates "G" the derivative of f with respect to 
% parameters.

%----------------
% Routine Checks
if exist('deval_soln','var')
	Y = deval(deval_soln,t);
end
n_params = length(params);
if ~exist('index','var')
	index = 1:n_params;
end

M = length(Y)/6;

%-----------------
% Unpack st vbls
un = Y(1:M);
uc = Y(1*M+1:2*M);
wn = Y(2*M+1:3*M);
wc = Y(3*M+1:4*M);
vn = Y(4*M+1:5*M);
vc = Y(5*M+1:6*M);

%-----------------
% Other matrices
Z = sparse(M,1);

%---------------
% Unpack params
params = num2cell(params);
[~,~,~,~,~,~,~,~,~,...
	gamma, psi,~, beta, phi, beta0, kappa, ~, K] = params{:};

%----------
% Groupings
if strcmp(stage,'interphase')
	an = An/Vn;
	ac = An/Vc;
    nH = 20;
    f = un.^nH./(K.^nH + un.^nH);
    df_by_dK = -nH*K^(nH-1)*un.^nH./(K^nH + un.^nH).^2;
elseif strcmp(stage,'mitosis')
	an = 0;
	ac = 0;
    f = sparse(M,1);
    df_by_dK = sparse(M,1);
end
am  = Am/Vc;

% =========================================================================
% G-elements
% =========================================================================

G = sparse(6*M, n_params);

if any(index == 1)
    % lambdaU
    G(:,1) = [Z; am*P*uc; repmat(Z,4,1)]; 
    
    %G(:,1) = [Z; am*P*uc; repmat(Z,6,1)]; 
    % Note this is 13+13+78 = 104. Or 1+1+6 = 8 (bcoz lenght = 8. Need to chanage)that remains constant
end

if any(index == 2)
    % lambdaW
    G(:,2) = [repmat(Z,3,1); am*P*wc; repmat(Z,2,1)]; 
    
	%G(:,2) = [repmat(Z,3,1); am*P*wc; repmat(Z,4,1)]; 
    % Note this is 39+13+52 = 104. Or 3+1+4 = 8 that remains constant
end

if any(index == 3)
    % lambdaV
	G(:,3) = [repmat(Z,5,1); am*P*vc]; 
    
    %G(:,3) = [repmat(Z,5,1); am*P*zc; repmat(Z,2,1)];
end

if any(index == 4)
    % sigmaU
	G(:,4) = [an*uc; -ac*uc; repmat(Z,4,1)]; 
    
    %G(:,4) = [an*uc; -ac*uc; repmat(Z,6,1)];
end

if any(index == 5)
    % sigmaW
	G(:,5) = [repmat(Z,2,1); an*wc; -ac*wc; repmat(Z,2,1)]; 
    
    %G(:,5) = [repmat(Z,2,1); an*wc; -ac*wc; repmat(Z,4,1)];
end

if any(index == 6)
    % sigmaV
	G(:,6) = [repmat(Z,4,1); an*vc; -ac*vc]; 
    
    %G(:,6) = [repmat(Z,4,1); an*zc; -ac*zc; repmat(Z,2,1)]; 
end

if any(index == 7)
    % muU
	G(:,7) = [-an*un; ac*un; repmat(Z,4,1)]; 
    
    %G(:,7) = [-an*un; ac*un; repmat(Z,6,1)]; 
end

if any(index == 8)
    % muW
	G(:,8) = [repmat(Z,2,1); -an*wn; ac*wn; repmat(Z,2,1)]; 
    
    %G(:,8) = [repmat(Z,2,1); -an*wn; ac*wn; repmat(Z,4,1)]; 
end

if any(index == 9)
    % muV
	G(:,9) = [repmat(Z,4,1); -an*vn; ac*vn];
    
    %G(:,9) = [repmat(Z,4,1); -an*zn; ac*zn; repmat(Z,2,1)];
end

if any(index == 10)
    % gamma
	G(:,10) = [-un.*vn; -uc.*vc; un.*vn; uc.*vc; -psi*un.*vn; -psi*uc.*vc];

   % G(:,10) = [-un.*zn; -uc.*zc; un.*zn; uc.*zc; -psi*un.*zn; -psi*uc.*zc; Z;Z]; % gamma
end

if any(index == 11)
    % psi
    g   = exp(-x.^2/2/phi^2);
	G(:,11) = [Z;Z;Z;Z; -(gamma*un.*vn - beta0*wn); -(gamma*uc.*vc - beta*g.*wc./(kappa + wc) - beta0*wc)]; % psi

   % G(:,11) = [Z;Z;Z;Z; -(gamma*un.*zn - beta0*wn); acs*omega*epsilon*xx-(gamma*uc.*zc - beta0*wc); Z;Z];
end

if any(index == 12)
	% alpha
    G(:,12) = [repmat(Z,5,1); -vc]; % alpha
    
  %  G(:,12) = [repmat(Z,5,1); -zc; repmat(Z,2,1)];
end

if any(index == 13)
    % beta
	G(:,13) = [Z; g.*wc./(kappa + wc); Z; -g.*wc./(kappa + wc); Z; psi*g.*wc./(kappa + wc) ]; % beta
    
%   G(:,13) = [repmat(Z,6,1); phi./(phi + x.^4); Z]; % beta
% 	G(:,13) = [repmat(Z,6,1); phi^4./(phi^4 + x.^4); Z]; % beta
end

if any(index == 14)
    % phi
    dBEc_by_dPhi = -beta*wc.*x.^2.*exp(-x.^2./2/phi^2)./((kappa + wc)*phi^3);
	G(:,14) = [Z; -dBEc_by_dPhi; Z; dBEc_by_dPhi; Z; -psi*dBEc_by_dPhi]; % phi
    
%   G(:,14) = [repmat(Z,6,1); beta*x.^4./(phi + x.^4).^2; Z]; % phi
% 	G(:,14) = [repmat(Z,6,1); beta*4*phi^3*x.^4./(phi^4 + x.^4).^2; Z]; % phi
% 	G(:,14) = [repmat(Z,6,1); beta*exp(-0.5*(x/phi).^2).*x.^2/phi.^3; Z]; % phi
% 	G(:,14) = [repmat(Z,6,1); beta*exp(-0.5*(x.^2/phi)).*0.5.*x.^2/phi.^2; Z]; % phi
end

if any(index == 15)
    % beta0
	G(:,15) = [wn; wc; -wn; -wc; psi*wn; psi*wc];  
    
   % G(:,15) = [wn; wc; -wn; -wc; psi*wn; psi*wc; repmat(Z,2,1)];
end


if any(index == 16)
    % kappa
    dBec_by_dkappa = (beta*wc.*exp(-x.^2./2/phi^2))./(kappa + wc).^2;
	G(:,16) = [Z; -dBec_by_dkappa; Z; dBec_by_dkappa; Z; -psi*dBec_by_dkappa]; 
    
    % nu - Make it kappa
	%G(:,16) = [Z;Z;Z; acs*epsilon*xx; Z;Z; xx; -xx];
end

% 
% if any(index == 17)
%     % delta
% 	G(:,17) = [repmat(Z,5,1); f]; 
%     
%     % omega - Make it delta
% 	%G(:,17) = [Z; acs*epsilon*xx; Z;Z;Z; acs*psi*epsilon*xx; xx; -xx];
% end
% 
% if any(index == 18)
%     % K
% 	G(:,18) = [repmat(Z,5,1); df_by_dK]; 
%     
%     % eta - Make it K
% 	% G(:,18) = [Z;Z;Z; -acs*epsilon*wc.*y; Z;Z; -y.*wc; y.*wc]; 
% end

%G = G(:,index);


%{
if any(index == 19)
    % rho - Take it off
	G(:,19) = [repmat(Z,6,1); -y; Z]; 
end

if any(index == 20)
     % epsilon - Take it off
	G(:,20) = [Z; acs*omega*xx; Z; -acs*(eta*wc.*y-nu*xx); Z; acs*psi*omega*xx; Z;Z];
end
%}









