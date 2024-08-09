% script_fake_quadstep
%
% This script is dumb in that it makes a bunch of fake parameter sets
% really, really close to a single one and does quad-step with them.
%
% (1) H_delt(k,k) (calculating two fdot's, differencing and div by delt)
% (2) H1(k,k) (fitting with fdot)
% (3) H2(k,k) (fitting without fdot...need more parameters)
% (4) d2fdtheta2 (calc f three times and take (f1 - 2*f + f2)/delt^2
% (5) fit f(theta_k) vs a parabola
% (6) Calc many fdots and fit them to a line
% (7) Calc fdot twice and do (fdot1 - fdot)/delt
% (8) Analytical


clear
close all
pth = 'C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\Mat';
matfilenames = readdir2(pth,'mat');
v = strfindDU(matfilenames,'2017-05-1');
matfilenames(~v) = [];
n = 20;
O = true(n); O = triu(O);
Z = zeros(n);
ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};


%
% Real data
%
[C1,dC1] = dlVenusData('nuclear');
[C2,dC2] = dlVenusData('total');
C = [C1;C2]; dC = [dC1;dC2];
Cmax = max(C./dC);

lu = [-4*ones(1,n); 4*ones(1,n)];
lu(1,14) = -1; lu(2,14) = 0; % lower and upper limits for phi
lambda = 105; mu = 21;
ub = repmat(lu(2,:),lambda,1);
lb = repmat(lu(1,:),lambda,1);


nfiles = length(matfilenames);
% nfiles = 2;
G = 50;
d50 = zeros(G-1,nfiles);
D_ = d50; fmin = d50; gradmin = d50;
beta = 1.1; % arb param
delt = 1e-6;
A_out = cell(nfiles,1);
for ii = 1:nfiles
	% 	ii1 = ii;
	load(matfilenames{ii})
	
	
% 	%
% 	% Quad step
% 	%
% 	np = 26;
% 	dx = randn(np,n);
% 	x = repmat(xb,np,1) + 1e-4*dx;
% 	Q = dist(x')/sqrt(n);
% 	d = mean(sum(Q)/(size(Q,1)-1)); % average pairwise distance.
% 	
% 	[f,~,Solnwt] = dsat_errorcalc(x,false);
% 	[fdot] = dsat_gradcalc(Solnwt);
% 	
% 	[xquad,X0,a] = quad_step(x,f,fdot,xb,d,beta,lb,ub);
% 	Xcen = 10.^xb'*10.^xb;
% 	
% 	H1 = Z;
% 	H1(O) = a;
% 	H1 = H1 + H1' - diag(diag(H1));
% 	H1 = H1./Xcen/(log(10)^2);  % Estimated Hessian
% 	
% 	
% 		
	%
	% Quad step 2
	%
	np = 300;
	dx = randn(np,n);
	x = repmat(xb,np,1) + 1e-4*dx;	
	[f2] = dsat_errorcalc(x,false);	
	[~,~,a2] = quad_step2(x,f2);

	H2 = Z;
	H2(O) = a2;
	H2 = H2 + H2' - diag(diag(H2));
	H2 = H2./Xcen/(log(10)^2); % Estimated Hessian 2
	
	
	
	%
	% Calculated Hessian
	%
	H = dsat_hessiancalc(xb);
	
	
	
% 	%
% 	% Single-param fits
% 	%
% 	kp = find(mean(fdot) > 1e-6);
% 	N = 11;
% 	Delt = linspace(-0.1,0.1,N)';
% 	f01 = zeros(length(Delt),1);
% 	theta = zeros(length(Delt),1);
% % 	Lam1 = NaN(n,1);
% 	Lam2 = NaN(n,1);
% 	Lam3 = NaN(n,1);
% 	d2fdtheta2 = NaN(n,1);	
% 	a_out = zeros(length(kp),6);
% 	count = 1;
% 	for k = kp
% 		
% 		%
% 		% Fit many values of f to a parabola, varying only one param.  Also
% 		% fit many values of fdot to a line, varying only that one param.
% 		%
% 		x1 = repmat(xb,length(Delt),1);
% 		for i = 1:length(Delt)
% 			delt = Delt(i);
% 			x1(i,k) = x1(i,k) + log10(1 + delt);
% 			theta(i) = 10.^x1(i,k);
% 		end
% 		[f01,~,Soln] = dsat_errorcalc(x1,false);
% 		[fdot01] = dsat_gradcalc(Soln);
% 		
% 		figger
% 		plot(theta,f01)
% 		
% 		A = [theta.^2 theta ones(length(theta),1)];
% 		lam2 = A \ f01;
% 		Lam2(k) = lam2(1);
% 		
% 		A = [theta ones(length(theta),1)];
% 		lam1 = A \ f01;
% 		
% 		lam3 = A \ fdot01(:,k);
% 		Lam3(k) = lam3(1);
% 		
% 		
% % 		%
% % 		% Calc 2nd derivative for one parameter by calc f three times
% % 		%
% % 		[f] = dsat_errorcalc(xb,false);
% % 		delt = 1e-4;
% % 		x1 = xb; x1(k) = x1(k) + log10(1 + delt);
% % 		[f1] = dsat_errorcalc(x1,false);
% % 		x2 = xb; x2(k) = x2(k) + log10(1 - delt);
% % 		[f2] = dsat_errorcalc(x2,false);
% % 		d2fdtheta2(k) = (f2 - 2*f + f1)/delt^2/(10^xb(k)).^2/(log(10)^2);
% 	
% 			
% 		
% 		
% 		a_out(count,:) = [Lam2(k) Lam3(k) H(k,k)];% d2fdtheta2(k) H1(k,k) H2(k,k)];
% 		count = count + 1;
% 	end
% 	A_out{ii} = a_out;
	disp(['difference in Hessians = ',num2str(norm(H2(:)-H(:)))])
	
	
end















