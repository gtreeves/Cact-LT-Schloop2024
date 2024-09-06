function data = fit_both_models_2step(Gs_data,B0,D20,W0)
%Runs two-step fitting procedures to fit two different models to data ACFs
%
%function data = fit_both_models_2step(Gs_data,B0,D20,W0)
% 
% This function takes the structure "Gs_data", which contains the
% data-derived ACFs for the total image, the nuclear fraction, the
% cytoplasmic fraction, and the nuclear channel (nuclear fraction), as well
% as the CCF, and fits the appropriate models to them. For the ACFs, it
% fits both the 1-comp diffusion model and the 2-comp (diffusion + binding)
% model, according to the 2-step procedure in which the fast direction cut
% is fit to a Gaussian first, which yields an accurate ACF amplitude, then
% fits the whole ACF to the two models described above. For the CCF case,
% there is a specific CCF model (for a zero-diffusion component) that
% includes the peak being off center and the width of the Gaussian being
% larger than expected for the green channel.
%
% The output is the structure "data", with several more constants filled
% in.
%
% "B0": if you want to force the background to be zero, then pass this
%	argument as B0 = 0.
% "D2": if you want to force the diffusivity for the two-component model
%	to be held constant during the fit, then pass this argument.
% "W0": if you want to hold the xy waist of the confocal volume to be held
%	fixed, then pass this argument.


if Gs_data.metadata.usemask
	append = {'tot','nuc','cyt','cc','nuc_ch'};
else
	append = {'tot'};
end
idx = 1:length(append);

%
% Copy over metadata from input variable "Gs_data" into output variable
% "data"
%
data.pth = Gs_data.pth;
data.filenameshort = Gs_data.filenameshort;
data.filename = Gs_data.filename;
if isfield(Gs_data,'basename')
	data.basename = Gs_data.basename;
end
data.genotype = Gs_data.genotype;
data.side = Gs_data.side;
data.metadata = Gs_data.metadata;
if isfield(Gs_data,'t')
	data.t = Gs_data.t;
end
data.nucsignal = Gs_data.nucsignal;
data.cytsignal = Gs_data.cytsignal;
if isfield(Gs_data,'totsignal')
	data.totsignal = Gs_data.totsignal;
	data.nuc_chsignal = Gs_data.nuc_chsignal;
end
if isfield(Gs_data,'nucvar')
	data.nucvar = Gs_data.nucvar;
	data.cytvar = Gs_data.cytvar;
end
if isfield(Gs_data,'totvar')
	data.totvar = Gs_data.totvar;
	data.nuc_chvar = Gs_data.nuc_chvar;
end

dr = data.metadata.dr;
wz = data.metadata.wz;
taup = data.metadata.taup*1e-6; % This assumes the pixel dwell time is in microseconds
tauL = data.metadata.tauL/1000; % Line time in milliseconds
H = data.metadata.H;
W = data.metadata.W;
xilimit = min(round(W/2),round(5*data.metadata.w0/sqrt(2)/data.metadata.dr)); % was: 6.4*...
etalimit = min(round(H/2),round(5*data.metadata.w0/sqrt(2)/data.metadata.dr));

for i = idx

	%
	% Unpack and validate data
	%
	GS = Gs_data.(['Gs',append{i}]);
	if strcmp(append{i},'cc')

		%
		% The CCF is a special case and has a completely different shape
		% (peak is not at zero, for example) and, as such, has completely
		% different fitting parameters.
		%
		data = fit_immobile_CCF(GS,data,B0,W0);
		continue

	else


		%
		% Preallocate
		%
		ngroups = size(GS,3);
		A = NaN(ngroups,1); B_PSF = NaN(ngroups,1); w01 = NaN(ngroups,1);
		B1 = NaN(ngroups,1); D = NaN(ngroups,1); 
		B2 = NaN(ngroups,1); D2 = NaN(ngroups,1); phi = NaN(ngroups,1); w02 = NaN(ngroups,1);
		R2_PSF = NaN(ngroups,1); R2_1 = NaN(ngroups,1); R2_2 = NaN(ngroups,1);
		delta68_PSF = NaN(3,ngroups); delta95_PSF = NaN(3,ngroups); J_PSF = NaN; 
		resnorm_PSF = NaN(ngroups,1); nu_PSF = NaN(ngroups,1); rmse_PSF = NaN(ngroups,1); 
		se_PSF = NaN(3,ngroups); exitflag_PSF = NaN(ngroups,1); %output_PSF = NaN(ngroups,1);
		delta68_1 = NaN(4,ngroups); delta95_1 = NaN(4,ngroups); J1 = NaN;
		resnorm1 = NaN(ngroups,1); nu1 = NaN(ngroups,1); rmse1 = NaN(ngroups,1); 
		se1 = NaN(4,ngroups); exitflag1 = NaN(ngroups,1); %output1 = NaN(ngroups,1);
		delta68_2 = NaN(5,ngroups); delta95_2 = NaN(5,ngroups); J2 = NaN; 
		resnorm2 = NaN(ngroups,1); nu2 = NaN(ngroups,1); rmse2 = NaN(ngroups,1);
		se2 = NaN(5,ngroups); exitflag2 = NaN(ngroups,1); %output2 = NaN(ngroups,1);

		%
		% Spatial coordinates
		%
		if any(isnan(GS(:)))
			continue
		end
		GS = GS(1:etalimit,1:xilimit,:);


		%
		% Run for loop over time points
		%
		for j = 1:ngroups
			Gs = GS(:,:,j);

			Gs_fast = Gs(1,:)';
			Gs1 = Gs_fast(3:xilimit);
			xi = 0:(xilimit-1);
			eta = (0:(etalimit-1))';

			Xi = repmat(xi,etalimit,1);
			Eta = repmat(eta,1,xilimit);

			% ---------------------------------------------------------------------
			% First step: we fit the fast direction to the PSF (so D = 0)
			% ---------------------------------------------------------------------

			%
			% Set up the startpoint and boundaries
			%
			% 	A = max(Gs(1,1),-Gs(1,1)); % A is the amplitude.
			A0 = Gs1(1,1);
			AL = min(A0/3,A0*3); % it will be A*3 if A is negative
			AU = max(A0/3,A0*3); % it will be A/3 if A is negative

			if exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 1 && B0 == 0
				B = B0;
				BU =  1e-6;
				BL = -1e-6;
			elseif exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 3
				B = B0(1); % "B" is the background
				BU = B0(3);
				BL = B0(2);
			else
				B = 0; % "B" is the background
				BU = 1e6;
				BL = -1e6;
			end
			if exist('W0','var') && isnumeric(W0) && ~isempty(W0)
				w0 = W0;
				w0U = w0 + 1e-4; % hold fixed
				w0L = w0 - 1e-4;
			else
				w0 = 0.3;
				w0L = 0.1;
				w0U = 0.7;
			end

			Lower = [AL BL w0L];
			Upper = [AU BU w0U];
			StartPoint = [A0 B w0];

			%
			% Set up fittype and options.
			%
			tols = 1e-16;
			options = optimset('Display','off',...
				'MaxFunEvals',5e2,...
				'MaxIter',5e2,...
				'TolX',tols,...
				'TolFun',tols,...
				'TolCon',tols ,...
				'UseParallel',false);

			fhandle = @ftn_fit_Gaussian;
			[p,resnorm_PSF(j),residual_PSF,exitflag_PSF(j),output_PSF(j),~,J_PSF] = ...
				lsqcurvefit(fhandle,StartPoint,xi(3:xilimit)',Gs1,Lower,Upper,options,dr);

			A(j) = p(1);
			B_PSF(j) = p(2);
			w01(j) = p(3);

			%
			% Post fit analysis: R-square and errobars
			%
			if ~isempty(J_PSF)
				nu_PSF(j) = length(J_PSF)-length(p); % DOF

				%
				% Calc goodness of fit
				%
				R2_PSF(j) = 1 - resnorm_PSF(j)/norm(Gs1-mean(Gs1))^2;

				%
				% Make errorbars
				%
				rmse_PSF(j) = norm(residual_PSF(:))/sqrt(nu_PSF(j));
				C = (J_PSF'*J_PSF) \ eye(length(p));
				se_PSF(:,j) = sqrt(diag(C))*rmse_PSF(j);

				% 68% errorbar (one sigma)
				alpha = 1 - 0.68;
				delta68_PSF(:,j) = se_PSF(:,j)*tinv(1-alpha/2,nu_PSF(j));

				% 95% errorbar
				alpha = 0.05;
				delta95_PSF(:,j) = se_PSF(:,j)*tinv(1-alpha/2,nu_PSF(j));

			end


			% ---------------------------------------------------------------------
			% Second step, diffusion model: we fit the slow dir to the 1-cpt ACF
			% ---------------------------------------------------------------------

			%
			% Set up the startpoint and boundaries - D0 is reserved for the 2-cpt
			% model. I don't know if we ever use it here.
			%
			D0 = 1; % micron^2/s
			DU = 50;
			DL = 1e-4;

			AL = min(A(j)*(1-1e-4),A(j)*(1+1e-4)); % depending on if A is negative
			AU = max(A(j)*(1-1e-4),A(j)*(1+1e-4)); % depending on if A is negative

			if exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 1 && B0 == 0
				BU =  1e-6;
				BL = -1e-6;
			elseif exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 3
				BU = B0(3);
				BL = B0(2);
			else
				BL = -1e6;
				BU = 1e6;
			end


			Lower = [min(A(j)*[1-1e-4 1+1e-4]) BL DL w01(j)*(1-1e-4)];
			Upper = [max(A(j)*[1-1e-4 1+1e-4]) BU DU w01(j)*(1+1e-4)];
			StartPoint = [A(j) 0 D0 w01(j)];

			%
			% Define dist coord XI and target Gs
			%
			Gs1 = Gs; % These two lines are for the fitting eta-only case
			Gs1 = Gs1(:); Gs1(1) = [];
			XI = {Xi Eta};


			%
			% Set up fittype and options.
			%
			tols = 1e-16;
			options = optimset('Display','off',...
				'MaxFunEvals',5e2,...
				'MaxIter',5e2,...
				'TolX',tols,...
				'TolFun',tols,...
				'TolCon',tols ,...
				'UseParallel',false);

			%
			% perform the fitting
			%
			fhandle = @ftn_fit_1comp_diffusion;
			[p,resnorm1(j),residual1,exitflag1(j),output1(j),~,J1] = ...
				lsqcurvefit(fhandle,StartPoint,XI,Gs1,Lower,Upper,options,...
				dr,wz,taup,tauL,true);
			B1(j) = p(2);
			D(j) = p(3);

			%
			% Post fit analysis: R-square and errobars
			%
			if ~isempty(J1)
				nu1(j) = length(J1)-length(p); % DOF

				%
				% Calc goodness of fit
				%
				R2_1(j) = 1 - resnorm1(j)/norm(Gs1-mean(Gs1))^2;

				%
				% Make errorbars
				%
				rmse1(j) = norm(residual1(:))/sqrt(nu1(j));
				C = (J1'*J1) \ eye(length(p));
				se1(:,j) = sqrt(diag(C))*rmse1(j);

				% 68% errorbar (one sigma)
				alpha = 1 - 0.68;
				delta68_1(:,j) = se1(:,j)*tinv(1-alpha/2,nu1(j));

				% 95% errorbar
				alpha = 0.05;
				delta95_1(:,j) = se1(:,j)*tinv(1-alpha/2,nu1(j));

			end



			% ---------------------------------------------------------------------
			% Second step, 2-cpt model: we fit the slow dir to the 2-cpt ACF
			% ---------------------------------------------------------------------

			%
			% Set up the startpoint and boundaries
			%
			if exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 1 && B0 == 0
				BU =  1e-6;
				BL = -1e-6;
			elseif exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 3
				BU = B0(3);
				BL = B0(2);
			else
				BL = -1e6;
				BU = 1e6;
			end
			if exist('D20','var') && isnumeric(D20) && ~isempty(D20)
				D201 = D20;
				DU = D20*(1 + 1e-4);
				DL = D20*(1 - 1e-4);
			else
				D201 = 1; % micron^2/s
				DU = 25;
				DL = 1e-4;
			end
			if exist('W0','var') && isnumeric(W0) && ~isempty(W0)
				% 		w0 = w01;
				w0 = W0;
				w0U = w0 + 1e-4; % hold fixed
				w0L = w0 - 1e-4;
			else
				w0 = 0.3; % This is the mean, from data
				sigma_w0 = 0.064; % Variation, from data
				w0U = w0 + 2*sigma_w0; % allow for 2 sigma variation?
				w0L = w0 - 2*sigma_w0;
			end
			phi0 = 0.5; phiL = 0; phiU = 1;

			Lower = [AL BL DL phiL w0L];
			Upper = [AU BU DU phiU w0U];
			StartPoint = [A(j) 0 D201 phi0 w0];

			%
			% Set up fittype and options.
			%
			tols = 1e-16;
			options = optimset('Display','off',...
				'MaxFunEvals',5e2,...
				'MaxIter',5e2,...
				'TolX',tols,...
				'TolFun',tols,...
				'TolCon',tols ,...
				'UseParallel',false);

			%
			% perform the fitting
			%
			fhandle = @ftn_fit_diffusion_w_immobile_2step;
			[p,resnorm2(j),residual2,exitflag2(j),output2(j),~,J2] = ...
				lsqcurvefit(fhandle,StartPoint,XI,Gs1,Lower,Upper,options,...
				dr,wz,taup,tauL,true);
			B2(j) = p(2);
			D2(j) = p(3);
			phi(j) = p(4);
			w02(j) = p(5);

			%
			% Post fit analysis: R-square and errobars
			%
			if ~isempty(J2)
				nu2(j) = length(J2)-length(p); % DOF

				%
				% Calc goodness of fit
				%
				R2_2(j) = 1 - resnorm2(j)/norm(Gs1-mean(Gs1))^2;

				%
				% Make errorbars
				%
				rmse2(j) = norm(residual2(:))/sqrt(nu2(j));
				C = (J2'*J2) \ eye(length(p));
				se2(:,j) = sqrt(diag(C))*rmse2(j);

				% 68% errorbar (one sigma)
				alpha = 1 - 0.68;
				delta68_2(:,j) = se2(:,j)*tinv(1-alpha/2,nu2(j));

				% 95% errorbar
				alpha = 0.05;
				delta95_2(:,j) = se2(:,j)*tinv(1-alpha/2,nu2(j));
				
			end

		end
	end
	
	% ---------------------------------------------------------------------
	% Concluding remarks
	% ---------------------------------------------------------------------	
	data.(['model',append{i}]) = '2-step';
	data.(['A',append{i}]) = A;
	data.(['B_PSF',append{i}]) = B_PSF;
	data.(['B',append{i}]) = B1;
	data.(['D',append{i}]) = D;
	data.(['B2',append{i}]) = B2;
	data.(['D2',append{i}]) = D2;
	data.(['phi',append{i}]) = phi;
	data.(['w01',append{i}]) = w01;
	data.(['w02',append{i}]) = w02;
	data.(['R2_PSF',append{i}]) = R2_PSF;
	data.(['R2_1',append{i}]) = R2_1;
	data.(['R2_2',append{i}]) = R2_2;
	
	data.(['fit',append{i}]).errorbar68_PSF = delta68_PSF;
	data.(['fit',append{i}]).errorbar95_PSF = delta95_PSF;
	data.(['fit',append{i}]).J_PSF = J_PSF;
	data.(['fit',append{i}]).resnorm_PSF = resnorm_PSF;
	data.(['fit',append{i}]).nu_PSF = nu_PSF;
	data.(['fit',append{i}]).rmse_PSF = rmse_PSF;
	data.(['fit',append{i}]).se_PSF = se_PSF;
	data.(['fit',append{i}]).output_PSF = output_PSF;
	data.(['fit',append{i}]).exitflag_PSF = exitflag_PSF;
		
	data.(['fit',append{i}]).errorbar68_1 = delta68_1;
	data.(['fit',append{i}]).errorbar95_1 = delta95_1;
	data.(['fit',append{i}]).J1 = J1;
	data.(['fit',append{i}]).resnorm1 = resnorm1;
	data.(['fit',append{i}]).nu1 = nu1;
	data.(['fit',append{i}]).rmse1 = rmse1;
	data.(['fit',append{i}]).se1 = se1;
	data.(['fit',append{i}]).output1 = output1;
	data.(['fit',append{i}]).exitflag1 = exitflag1;
	
	data.(['fit',append{i}]).errorbar68_2 = delta68_2;
	data.(['fit',append{i}]).errorbar95_2 = delta95_2;
	data.(['fit',append{i}]).J2 = J2;
	data.(['fit',append{i}]).resnorm2 = resnorm2;
	data.(['fit',append{i}]).nu2 = nu2;
	data.(['fit',append{i}]).rmse2 = rmse2;
	data.(['fit',append{i}]).se2 = se2;
	data.(['fit',append{i}]).output2 = output2;
	data.(['fit',append{i}]).exitflag2 = exitflag2;


end