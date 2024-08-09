function Fdot = dsat2_gradcalc(Soln,fhandle,reltol,abstol)

% This function calculates the gradient of the error 
% Input "fhandle" should be:
%	@perturb_dsat (default if not passed)
%	@perturb_dsat_IE2
%	@perturb_dsat_ode15s
%	@perturb_dsat_RK4


% 
% Check inputs
%
if ~exist('fhandle','var') || isempty(fhandle) || isnumeric(fhandle) || islogical(fhandle)
	fhandle = @perturb_dsat;
end
if ~exist('reltol','var') || isempty(reltol)
	reltol = 1e-7;
end
if ~exist('abstol','var') || isempty(abstol)
	abstol = 1e-7;
end


% Misc. initial procedures
ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};


% Experimental data from dl-Venus
% [C1,dC1]    = dlVenusData('nuclear');
% [C2,dC2]    = dlVenusData('total');
% C           = [C1; C2];
% dC          = [dC1; dC2];
% Cmax        = max(C./dC);
[C,dC]  = dlVenusData('nuclear');
C       = C/1000;


% 
% LOOP: Initialize and run the loop over every parameter set 
%
nSets      = length(Soln);
nParams    = length([Soln{1}.params,Soln{1}.addParams]);
Fdot       = zeros(nSets, nParams);

parfor ii = 1:nSets
	data        = Soln{ii};
    
    %DYdtheta = fhandle(data, modelname, reltol, abstol);
    try
	    DYdtheta = fhandle(data, reltol, abstol);
	catch ME
		disp(ME.message)
		Fdot(ii,:) = NaN(1,nParams);
		continue
    end


	% Simulation data as a column vector
    D = [data.nuclearDorsal.NC11(:);
         data.nuclearDorsal.NC12(:);
         data.nuclearDorsal.NC13(:);
         data.nuclearDorsal.NC14(:)];

 
	% derivatives of our simulation data
	Dprime  = cell(6,1);
	count   = 1;
	for i = [3 5 7 9]
		nc          = ncs{i};
		D1          = DYdtheta.nuclearDorsal.(nc);
		D2          = DYdtheta.totalDorsal.(nc);
		[m1,n1,o1]  = size(D1);
		D1          = reshape(D1,m1*n1,o1);
		D2          = reshape(D2,m1*n1,o1);
		
		Dprime{count}   = D1; 
		%Dprime{count+4} = D2;
		count           = count + 1;
	end
	Dprime  = cell2mat(Dprime);

    
	% Error calculation
	Beta    = mean(D.*C)/mean(D.^2);
	epsln   = (C-Beta*D);
	gamma   = sum(epsln.^2);

    
	% Gradient calculation
	%dBetadtheta     = (sum(repmat(C,1,nparams).*Dprime)*sum(D.^2) - ...
	%	2*sum(C.*D)*sum(repmat(D,1,nparams).*Dprime))/  sum(D.^2).^2;
	%depslndtheta    = -(D*dBetadtheta + Beta*Dprime)./repmat(dC,1,nparams);
	%gammadot        = 2*sum(repmat(epsln,1,nparams).*depslndtheta);
    %gammadot        = sum(repmat(epsln,1,nparams).*depslndtheta)./sqrt(sum(epsln));
	%Fdot(ii,:)      = 1/(sqrt(length(C))*Cmax*2*sqrt(gamma))*gammadot;	
    %Fdot(ii,:)      = 1/(sqrt(length(C))*Cmax*2*sqrt(gamma))*gammadot;    
    %Fdot(ii,:)      = 1/(sqrt(length(C))*2*sqrt(gamma))*gammadot;
    
    
    dBetadtheta     = (sum(repmat(C,1,nParams).*Dprime)*sum(D.^2) - ...
		2*sum(C.*D)*sum(repmat(D,1,nParams).*Dprime))/  sum(D.^2).^2;
    depslndtheta    = -(D*dBetadtheta + Beta*Dprime);
	gammadot        = 2*sum(repmat(epsln,1,nParams).*depslndtheta);

    
    Fdot(ii,:)      = gammadot;
end










