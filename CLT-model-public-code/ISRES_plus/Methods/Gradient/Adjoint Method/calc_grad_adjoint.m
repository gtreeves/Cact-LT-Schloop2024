function [grad,Grad,Qall,L,lamsoln] = calc_grad_adjoint(protein)
% This function evaluates the gradeint of the decon model by solving the
% adjoint equation individually for every nuclear cycle. Between nuclear
% cycles, lambda is interpolated onto the nc's grid such that the
% dimensions of concentrations and lambda match. 
%
% Opinion: This should be the most accurate way to obtain the gradient via
% the adjoint method. 

%load('sample.mat')

% tolerances
reltol = 1e-2;
abstol = 1e-2;


% extract data from model soln
Tspan   = protein.T;
M       = protein.M;
params  = protein.params;
nparams = length(params);


% dlvenus data
[C ,~,Cstruct,~,Cmat]    = dlVenusData('nuclear','SingleMatrix',true);
%[~,~,~,dCstruct]  = dlVenusData('nuclear','IncludeMitosis',true);
Cstruct.('NC10')  = zeros(size(protein.dlNuc.('NC10')));
Cstruct           = orderfields(Cstruct);
C = C/1000;


% C  = [Cmat.NC11, Cmat.NC12, Cmat.NC13, Cmat.NC14]/1000;
% C = C(:);
Cmat.('NC10')  = zeros(size(protein.dlNuc.('NC10')));
Cmat           = orderfields(Cmat);



% names
names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt'};
num   = length(names);
ncs   = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
nncs  = length(ncs);


% Simulation data as a column vector
D      = [ protein.nuclearDorsal.NC11(:);
           protein.nuclearDorsal.NC12(:);
           protein.nuclearDorsal.NC13(:);
           protein.nuclearDorsal.NC14(:)];
sumD2  = sum(D.*D);
sumDC  = sum(D.*C);
nlen   = length(D);
beta   = mean(D.*C)/mean(D.^2);
         


%
% Solve adjoint equations backwards in time!
%
options = odeset('RelTol',reltol,'AbsTol',abstol); 
grad = 0;
for i = fliplr(1:nncs)
    
    % get nuclear cycle
    nc = ncs{i};
    
	% stage-dependent factors
	m     = M.(nc);
	tspan = Tspan.(nc);
    x     = linspace(0,1,m)';
    
        
    ii = [1:m, 1:m-1, 2:m];     %rows:    main diagonal, right diagonal, left diagonal
    jj = [1:m, 2:m, 1:m-1];     %columns: main diagonal, right diagonal, left diagonal
    vv = [-2*ones(1,m), 2, ones(1,2*m-4), 2]; 
    P  = sparse(ii,jj,vv);
    
    
	% Initial conditions
    % First, get state variables from the solved equations!
    ns = 6*m;
    nt = length(tspan);
    Y  = zeros(ns,nt);
    for j = 1:num
        Y((j-1)*m+1:j*m,:) = protein.(names{j}).(nc);       
    end  
    
    % check by plotting
    %{
    un = Y(1:m,:);
    uc = Y(1*m+1 :2*m,:);
    wn = Y(2*m+1 :3*m,:);
    wc = Y(3*m+1 :4*m,:);
    vn = Y(4*m+1 :5*m,:);
    vc = Y(5*m+1 :6*m,:); 
    
    plot(x,(un+wn));
    plot(x,uc);
    plot(x,wc);
    plot(x,vn);
    plot(x,vc);
    %}

    
    % stage-dependent factors
    
    if i==nncs
       lambda0 = zeros(ns,1);
       stage = 'interphase';
       [An,Am,~,Vn,Vc] = nuclearSize(1, 'static', m, stage);
       Yobs  = Cstruct.(nc)/1000;       % This is Dl-Venus data.    
     
    elseif mod(i,2) == 0  
        x_old   = linspace(0,1,m_old)'; 
        
        lam1 = interp1(x_old,lam10,x,'spline');
        lam2 = interp1(x_old,lam20,x,'spline');
        lam3 = interp1(x_old,lam30,x,'spline');
        lam4 = interp1(x_old,lam40,x,'spline');
        lam5 = interp1(x_old,lam50,x,'spline');
        lam6 = interp1(x_old,lam60,x,'spline');  
        
%         lam1    = ((Vn+Vc)*lam1)./Vc;
%         lam2    = ((Vn+Vc)*lam2)./Vc;
%         lam3    = ((Vn+Vc)*lam3)./Vc;
%         lam4    = ((Vn+Vc)*lam4)./Vc;
%         lam5    = ((Vn+Vc)*lam5)./Vc;
%         lam6    = ((Vn+Vc)*lam6)./Vc; 

%         lam10 = interp1(x_old,lam10,x);
%         lam20 = lam10;
%         lam30 = interp1(x_old,lam30,x);
%         lam40 = lam30;
%         lam50 = interp1(x_old,lam50,x);
%         lam60 = lam50;
        
%         lam20 = interp1(x_old,lam20,x,'spline');
%         lam10 = lam20;         
%         lam40 = interp1(x_old,lam40,x,'spline');
%         lam30 = lam40;
%         lam60 = interp1(x_old,lam60,x,'spline');    
%         lam50 = lam60;
%         
%         lam10 = zeros(m,1);
%         lam20 = interp1(x_old,lam20,x,'spline');
%         lam30 = lam10;
%         lam40 = interp1(x_old,lam40,x,'spline');
%         lam50 = lam10;
%         lam60 = interp1(x_old,lam60,x,'spline'); 
        
        
%         lam1 = zeros(m,1);
%         lam2 = (Vn*lam10 + Vc*lam20)/(Vn + Vc);
%         lam3 = zeros(m,1);
%         lam4 = (Vn*lam30 + Vc*lam40)/(Vn + Vc);
%         lam5 = zeros(m,1);
%         lam6 = (Vn*lam50 + Vc*lam60)/(Vn + Vc);        
        %

%         lam1    = ((Vn)*lam10)./(Vc);
%         lam2    = ((Vn)*lam20)./(Vc);
%         lam3    = ((Vn)*lam30)./(Vc);
%         lam4    = ((Vn)*lam40)./(Vc);
%         lam5    = ((Vn)*lam50)./(Vc);
%         lam6    = ((Vn)*lam60)./(Vc); 
        
        lambda0 = [lam1; lam2; lam3; lam4; lam5; lam6];      
     
        %lambda0 = [lam10; lam20; lam30; lam40; lam50; lam60];
		stage   = 'mitosis';
		[An,Am,~,Vn,Vc] = nuclearSize(1,'static',m,stage);

        Yobs    = zeros(m,nt);
    else  
        

%         
%         x_old   = linspace(0,1,m_old)'; 
%         lambda0 = reshape(lambda(:,1),m_old,6);
%         lambda0 = interp1(x_old,lambda0,x,'spline');
%         %lambda0 = interp1(x_old,lambda0,x);
%         lambda0 = lambda0(:);
        Vc2 = Vc;
        Vn2 = Vn;
        stage   = 'interphase';
		[An,Am,~,Vn,Vc] = nuclearSize(1, 'static', m, stage);
% %         
%         lam1    = lam10;
%         lam2    = lam20;
%         lam3    = lam30;
%         lam4    = lam40;
%         lam5    = lam50;
%         lam6    = lam60;  
%         
%         lam1    = ((Vn+Vc)*lam10)./Vc;
%         lam2    = ((Vn+Vc)*lam20)./Vc;
%         lam3    = ((Vn+Vc)*lam30)./Vc;
%         lam4    = ((Vn+Vc)*lam40)./Vc;
%         lam5    = ((Vn+Vc)*lam50)./Vc;
%         lam6    = ((Vn+Vc)*lam60)./Vc; 
% %         
        lam1    = (Vc2)*lam10./Vc;
        lam2    = (Vc2)*lam20./Vc;
        lam3    = (Vc2)*lam30./Vc;
        lam4    = (Vc2)*lam40./Vc;
        lam5    = (Vc2)*lam50./Vc;
        lam6    = (Vc2)*lam60./Vc; 
        
%         lam1    = ((Vn2+Vc2)*lam10)./(Vc+Vn);
%         lam2    = ((Vn2+Vc2)*lam20)./(Vc+Vn);
%         lam3    = ((Vn2+Vc2)*lam30)./(Vc+Vn);
%         lam4    = ((Vn2+Vc2)*lam40)./(Vc+Vn);
%         lam5    = ((Vn2+Vc2)*lam50)./(Vc+Vn);
%         lam6    = ((Vn2+Vc2)*lam60)./(Vc+Vn); 
        
        
%         lam1    = ((Vn)*lam10)./(Vc+Vn);
%         lam2    = ((Vn+Vc)*lam20)./Vc;
%         lam3    = ((Vn)*lam30)./(Vc+Vn);
%         lam4    = ((Vn+Vc)*lam40)./Vc;
%         lam5    = ((Vn)*lam50)./(Vc+Vn);
%         lam6    = ((Vn+Vc)*lam60)./Vc; 

%         lam1 = lam10;
%         lam2 = (Vn*lam10 + Vc*lam20)/(Vn + Vc);
%         lam3 = lam30;
%         lam4 = (Vn*lam30 + Vc*lam40)/(Vn + Vc);
%         lam5 = lam50;
%         lam6 = (Vn*lam50 + Vc*lam60)/(Vn + Vc); 
        
%         lam2    = ((Vn+Vc)*lam20 - Vn*zeros(m,1))./Vc;
%         lam1    = lam2;
%         lam4    = ((Vn+Vc)*lam40 - Vn*zeros(m,1))./Vc;
%         lam3    = lam4;
%         lam6    = ((Vn+Vc)*lam60 - Vn*zeros(m,1))./Vc;
%         lam5    = lam6;
        
%         lam1    = lam10;
%         lam2    = ((Vn+Vc)*lam20 - Vn*lam10)./Vc;
%         lam3    = lam30;
%         lam4    = ((Vn+Vc)*lam40 - Vn*lam30)./Vc;
%         lam5    = lam50;
%         lam6    = ((Vn+Vc)*lam60 - Vn*lam50)./Vc; 
       
       lambda0 = [lam1; lam2; lam3; lam4; lam5; lam6];
        
        Yobs    = Cstruct.(nc)/1000;       % This is Dl-Venus data. 
    end
   

   % check by plotting
   %{
   figure(1)
   subplot(2,4,i)
   lam_old = reshape(lambda(:,1),m_old,6);
   lam_new = reshape(lambda0,m,6);
   plot(xgrid,lam_new,'linewidth',2); hold on; plot(xgrid_old,lam_old,'*');hold off
   title(nc); legend(['un';'uc';'wn';'wc';'vn';'vc'])
   set(gca,'fontsize',14)
   %}
     

    %lambda0 = zeros(ns,1);
    % dudY matrix
    spDiag = spdiags(ones(m,1),0,m,m);
    spMat  = sparse(m,m);
    dudY   = [spDiag, spMat, spDiag, spMat, spMat, spMat];
    pp     = spline(tspan,Yobs);
    %dpp    = spline(tspan,dYobs);

    
    % check by plotting
    %{
    figure(2)
    tcheck    = linspace(tspan(1),tspan(end),1000);
    yobsCheck = ppval(pp,tcheck);
    plot(tspan,Yobs','.'); hold on; plot(tcheck,yobsCheck); hold off;
    %}
    
    
    % solve ode
    protein_soln_nc = protein.soln.(nc);
    soln    = ode15s(@dlambdaT_by_dt2,flipud(tspan),lambda0,options, ...
              m,dudY,pp,protein_soln_nc,@lambda_jac,params,Vn,Vc,An,Am,stage,sumD2,sumDC,...
              beta,nlen,tspan,Yobs,P);

   
    if soln.stats.nsteps == 0
		error('protein.stats.nsteps was zero')
    end
    

    % get lambda and store it
    [lambda,lamdadot]  = deval(soln, tspan);
    f = zeros(size(lamdadot));
    for j = 1:length(tspan)
        f(:,j) = dlambdaT_by_dt2(tspan(j),lambda(:,j),m,dudY,pp,protein_soln_nc,@lambda_jac,...
            params,Vn,Vc,An,Am,stage,sumD2,sumDC,beta,nlen,tspan,Yobs,P);
    end
    lamsoln.lamdadot.(nc)           = lamdadot;
    lamsoln.F.(nc)                  = f;  
    lamsoln.F_minus_lamdadot.(nc)   = f - lamdadot;
    
    
    
    
    lam10   = lambda(0*m+1:1*m,1);
    lam20   = lambda(1*m+1:2*m,1);
    lam30   = lambda(2*m+1:3*m,1);
    lam40   = lambda(3*m+1:4*m,1);
    lam50   = lambda(4*m+1:5*m,1);
    lam60   = lambda(5*m+1:6*m,1);
    %Lambda.(nc) = lambda;
  
    
    
    % check by plotting
    %{
    figure(3)
    subplot(3,4,i)
    %lambda_plot = permute(lambda,[2,1]);
    %lambda_plot = reshape(lambda_plot,length(tspan),m,[]);
    plot(tspan,lambda); hold on; xlabel('time (s)'); ylabel('lambda'); 
    %plot(lambda)
    title(nc); set(gca,'fontsize',14);
    %}  
    
    delt   = tspan(2) - tspan(1);
    tspan  = linspace(tspan(1),tspan(end),500);
    lambda = deval(soln, tspan);

    %
    % Integrate forward in time!
    %
    %if mod(i,2)~=0
        Q = zeros(nt,nparams);
        for j = 1:length(tspan)
            dfdp    = G_dsat(tspan(j), [], params, x, P, Vn, Vc, An, Am, stage,protein_soln_nc);
            Q(j,:) = -lambda(:,j)'*dfdp;
        end
        grad        = grad + trapz(tspan,Q)/delt;
        Grad(i,:)   = trapz(tspan,Q)/delt;
        Qall{i}     = Q;
        L{i}        = lambda;
   % end
    %delt = tspan(2) - tspan(1);
    %grad        = grad + trapz(tspan,Q)/delt;
    
    
    
    % prep for the next iteration
    m_old = m;
end
grad = grad';
Grad = Grad';
Grad = [Grad, sum(Grad,2)];

disp('done')



