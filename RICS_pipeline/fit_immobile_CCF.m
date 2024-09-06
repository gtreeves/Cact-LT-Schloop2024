function data = fit_immobile_CCF(Gscc,data,B0,W0)
%Fits a model of an immobile CCF to the CCF data stored in Gscc.
%
%function data = fit_immobile_CCF(Gscc,data,B0,W0)
%
% This function takes the CCF, "Gscc", and the structure "data", which
% contains the metadata and the tertiary data, and fits a model of an
% immobile CCF to it.
%
% The output is the structure "data", with several more constants filled
% in.

[h,w,ngroups] = size(Gscc);

xilimit = 32;


%
% Unpack data
%
if ~isstruct(data)
	% then "data" is just "dr"
	dr = data;
	clear data
else
	dr = data.metadata.dr;
end

%
% Spatial coordinates
%
if ~exist('xilimit','var')
	xi = repmat(0:(w-1),h,1);
	eta = repmat((0:(h-1))',1,w);
	xilimit = h;
else
	xi = repmat(0:(xilimit-1),xilimit,1);
	eta = repmat((0:(xilimit-1))',1,xilimit);

end



A = NaN(1,ngroups);
B = NaN(1,ngroups);
w0 = NaN(1,ngroups);
dx = NaN(1,ngroups);
dy = NaN(1,ngroups);
R2 = NaN(1,ngroups);
delta1 = NaN(5,ngroups);
delta2 = NaN(5,ngroups);
J = NaN; % J is too big to keep saving. We never use it posthoc anyway
resnorm = NaN(1,ngroups);
nu = NaN(1,ngroups);
rmse = NaN(1,ngroups);
se = NaN(5,ngroups);
exitflag = NaN(1,ngroups);

for ii = 1:ngroups

	if any(isnan(Gscc(:)))
		break
	end

	%
	% To set up initial guess, use parabolic interpolation to estimate what
	% the first pixel "could be"
	%
	xi1 = 2:4;
	idx = xi1 + 1;
	A1 = [xi1'.^2 ones(length(idx),1)];
	b = Gscc(1,idx,ii)';
	x = A1 \ b;
	c = x(2);

	%
	% Set up the startpoint and boundaries
	%
	% A = Gscc(1,1);
	Ag = c;
	AL = min(Ag/3,Ag*3); % it will be A*3 if A is negative
	AU = max(Ag/3,Ag*3); % it will be A/3 if A is negative

	if exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 1 && B0 == 0
		Bg = B0;
		BU =  1e-6;
		BL = -1e-6;
	elseif exist('B0','var') && isnumeric(B0) && ~isempty(B0) && length(B0) == 3
		Bg = B0(1); % "B" is the background
		BU = B0(3);
		BL = B0(2);
	else
		Bg = 0; % "B" is the background
		BU = 1e6;
		BL = -1e6;
	end
	if exist('W0','var') && isnumeric(W0) && ~isempty(W0)
		w0g = W0;
		w0U = w0g + 1e-4; % hold fixed
		w0L = w0g - 1e-4;
	else
		w0g = 0.3;
		w0L = 0.1;
		w0U = 0.7;
	end

	dxg = 0.1; dxL = 0; dxU = 0.2;
	dyg = 0.1; dyL = 0; dyU = 0.2;


	StartPoint = [Ag Bg w0g dxg dyg];
	Lower = [AL BL w0L dxL dyL];
	Upper = [AU BU w0U dxU dyU];

	%
	% Define dist coord XI and target Gs
	%
	if exist('fitxionly','var') && fitxionly
		Gscc1 = Gscc(1,2:xilimit+1,ii)'; % These two lines are for the fitting xi-only case
		XI = {0:(xilimit-1) 0};
	else
		Gscc1 = Gscc(1:xilimit,1:xilimit,ii);
		Gscc1(1,1:2) = [NaN NaN];
		Gscc1 = Gscc1(:);
		Gscc1(isnan(Gscc1)) = [];
		XI = {xi eta};
	end



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
	fhandle = @ftn_fit_immobile_CCF;
	[p,resnorm(ii),residual,exitflag(ii),output(ii),~,J] = ...
		lsqcurvefit(fhandle,StartPoint,XI,Gscc1,Lower,Upper,options,dr,true);
	A(ii) = p(1);
	B(ii) = p(2);
	w0(ii) = p(3);
	dx(ii) = p(4);
	dy(ii) = p(5);


	%
	% Post fit analysis: R-square and errobars
	%
	if ~isempty(J)
		nu(ii) = length(J)-length(p); % DOF

		%
		% Calc goodness of fit
		%
		R2(ii) = 1 - resnorm(ii)/norm(Gscc1-mean(Gscc1))^2;

		%
		% Make errorbars
		%
		rmse(ii) = norm(residual(:))/sqrt(nu(ii));
		C = (J'*J) \ eye(length(p));
		se(:,ii) = sqrt(diag(C))*rmse(ii);

		% 68% errorbar (one sigma)
		alpha = 1 - 0.68;
		delta1(:,ii) = se(:,ii)*tinv(1-alpha/2,nu(ii));

		% 95% errorbar
		alpha = 0.05;
		delta2(:,ii) = se(:,ii)*tinv(1-alpha/2,nu(ii));
	end

end


% -----------------------------------------------------------------
% Concluding remarks
% -----------------------------------------------------------------
data.modelcc = 'Immobile CCF';
data.Acc = A;
data.Bcc = B;
data.w01cc = w0;
data.dxcc = dx;
data.dycc = dy;
data.R2_PSFcc = R2;

data.fitcc.errorbar68_PSF = delta1;
data.fitcc.errorbar95_PSF = delta2;
data.fitcc.J = J;
data.fitcc.resnorm = resnorm;
data.fitcc.nu = nu;
data.fitcc.rmse = rmse;
data.fitcc.se = se;
if exist('output','var')
	data.fitcc.output = output;
else
	data.fitcc.output = NaN;
end
data.fitcc.exitflag = exitflag;




