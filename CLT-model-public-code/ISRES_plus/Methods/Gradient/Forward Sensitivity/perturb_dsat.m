function DYdtheta = perturb_dsat(soln,reltol,abstol)

% This function takes a previous run of the deconvolution model with Toll
% saturation (dsat) and calculates the derivative dYdtheta, where Y is the
% set of state variables and theta is the list of 15 params, from NC 10
% to NC 14.
% params    : a vector of parameters in linear space.
% dl0,cact0 : the initial concentrations (normalized) of dl,cact (resp)
% T         : a vector that delineates time points between mitoses and 
%             interphases.
%--------------------------------------------------------------------------
% Original Description
% This function takes a previous run of the deconvolution model with Toll
% saturation (dsat) and calculates the derivative dYdtheta, where Y is the
% set of state variables and theta is the list of all 20 params, from NC 10
% -- NC 14.
%
% params: a vector of parameters in linear space.
% dl0,cact0: the initial concentrations (normalized) of dl,cact (resp)
% T: a vector that delineates time points between mitoses and interphases.
%
% options = odeset('Jacobian',@jac_dsat2);
% options = odeset('RelTol',1e-1,'AbsTol',1e-2,'Jacobian',@jac_dsat2);


%
% Check inputs
%
if ~exist('reltol','var')
	reltol = 1e-7;
end
if ~exist('abstol','var')
	abstol = 1e-7;
end
options = odeset('RelTol',reltol,'AbsTol',abstol,'Jacobian',@jac);


%
% unpack params and ensure that length(params) = 18
%
params = [soln.params, soln.addParams];


% Number of species
ns = 6;


% Unpacking structure fields
nParams  = length(params);
M        = soln.M;               % Numbers of nuclei in each NC
Tspan    = soln.T;               % Time duration of each NC
X        = soln.X;
names    = fieldnames(soln)'; 
names    = names(1:ns);



%
% Looping through each nuclear cycle - mitosis & interphase
%
% Just take one of the fields, here - dlNuc, to calculate how many times
% to loop.
ncs     = fieldnames(soln.dlNuc)';      
nncs    = length(ncs);
m       = 0;
for i = 1:nncs
	nc = ncs{i};


	% Stage-dependent factors
	m_old   = m;
	m       = M.(nc);
	tspan   = Tspan.(nc);
	nt      = length(tspan);
	x       = X.(nc);
	Y       = zeros(ns*m,nt);
	soln1   = soln.soln.(nc);
	

	% Transport matrix
	e       = ones(m,1);
	P       = spdiags([e -2*e e],[-1 0 1], m, m);
	P(1,2)  = 2; 
    P(m,m-1)= 2; %#ok<SPRIX>
	

	% Protein values (Y)	
    for j = 1:ns
		Y((j-1)*m+1:j*m,:) = soln.(names{j}).(nc);
    end


	% Initial conditions
    if i == 1                   % the ultimate initial conditions
		stage = 'interphase';
		[An, Am, ~, Vn, Vc] = nuclearSize(1, 'static', m, stage);		
		dun0 = zeros(m,1,nParams);
        duc0 = zeros(m,1,nParams);
        dwn0 = zeros(m,1,nParams);
        dwc0 = zeros(m,1,nParams);
        dvn0 = zeros(m,1,nParams);
        dvc0 = zeros(m,1,nParams);
        
	elseif mod(i,2) == 0        % transition from interphase to mitosis
		dun0  = zeros(m,1,nParams); 
        duc0  = (Vn*dun(:,end,:) + Vc*duc(:,end,:))/(Vn + Vc);
        dwn0  = dun0; 
		dwc0  = (Vn*dwn(:,end,:) + Vc*dwc(:,end,:))/(Vn + Vc);
		dvn0  = dun0;
		dvc0  = (Vn*dvn(:,end,:) + Vc*dvc(:,end))/(Vn + Vc);     
		stage = 'mitosis';
		[An, Am, ~, Vn, Vc] = nuclearSize(1,'static',m,stage);
		
	elseif mod(i,2) ~= 0        % transition from mitosis to interphase		
		stage = 'interphase';
		[An, Am, ~, Vn, Vc] = nuclearSize(1, 'static', m, stage);
    
		% Interpolating onto the new x-mesh.  Also, Nuc species = Cyt		
		x_old = linspace(0,1,m_old)';	
		duc0  = interp1(x_old,duc(:,end,:),x); 
        dun0  = duc0;
		dwc0  = interp1(x_old,dwc(:,end,:),x); 
        dwn0  = dwc0;
        dvc0  = interp1(x_old,dvc(:,end,:),x); 
        dvn0  = dvc0;	
    end
    
	dYdtheta0 = [dun0; duc0; dwn0; dwc0; dvn0; dvc0];
	

	% Solving and Evaluating
	fhandle         = @ftn_dsat_dYdtheta;
	[~,dYdtheta]    = ode15s(fhandle, tspan, dYdtheta0, options, ...
		params, x, P, Vn, Vc, An, Am, stage, soln1);
	dYdtheta        = reshape (dYdtheta, [nt,ns*m,nParams]);
	dYdtheta        = permute (dYdtheta, [2 1 3]);
 %  dYdtheta        = dYdtheta (:, :, nparams);
   

	% Misc accounting
	dun = dYdtheta(1:m,:,:);
	duc = dYdtheta(1*m+1:2*m, :, :);
	dwn = dYdtheta(2*m+1:3*m, :, :);
	dwc = dYdtheta(3*m+1:4*m, :, :);
    dvn = dYdtheta(4*m+1:5*m, :, :);
	dvc = dYdtheta(5*m+1:6*m, :, :);
    
    for j = 1:ns
		DYdtheta.(names{j}).(nc) = dYdtheta((j-1)*m+1:j*m, :, :);
    end

    
    % Keeping derivatives of Nuclear and Total Dorsal
	DYdtheta.('nuclearDorsal').(nc) = ...
		DYdtheta.dlNuc.(nc) + DYdtheta.dlCactNuc.(nc);
	DYdtheta.('totalDorsal').(nc)   = ...
		(Vn*(DYdtheta.dlNuc.(nc) + DYdtheta.dlCactNuc.(nc)) + ...
		Vc*(DYdtheta.dlCyt.(nc) + DYdtheta.dlCactCyt.(nc))) / (Vc + Vn);
end










 


 