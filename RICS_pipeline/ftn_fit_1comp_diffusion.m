function F = ftn_fit_1comp_diffusion(P,XI,dr,wz,taup,tauL,long)
%Calculates a one componenet diffusion model for fitting to ACF data
%
%function F = ftn_fit_1comp_diffusion(P,XI,dr,wz,taup,tauL,long)
%
% This function calculates the theoretical ACF for a one-component pure
% diffusion model. It is used by the "lsqcurvefit" function to fit to the
% data ACF to give us the diffusivity D, A (the value of the ACF at xi =
% eta = 0), and background B. We also allow the PSF waist size, w0, to vary
% slightly for robustness. Although usually A and w0 are already determined
% by the first step in the two-step fitting, and B is tightly constrained
% around zero. So, effectively, the main adjustable parameter is D.
A = P(1);
B = P(2);
D = P(3);
w0 = P(4);

xi = XI{1}; eta = XI{2}; 

%
% Check wz input
%
if ~isnumeric(wz) || isnan(wz) || isempty(wz)
	wz = 3*w0;
end

%
% Check taup and tauL inputs to make sure they're in seconds instead of µs
% or ms, respectively.
%
if taup > 1e-4 % assume that oops I passed taup in µs instead of s
	taup = taup*1e-6;
end
if tauL > 1e-1 % assume that oops I passed tauL in ms instead of s
	tauL = tauL*1e-3;
end

%
% Check the size of xi and eta, to put them into 7
% 
[mxi,nxi] = size(xi);
[meta,neta] = size(eta);
if mxi > 1 && nxi == 1 && meta > 1 && neta == 1 % both col vecs, not necc same length
	xi = repmat(xi',meta,1);
	eta = repmat(eta,1,mxi);
elseif mxi == 1 && nxi > 1 && meta > 1 && neta == 1 % xi is row, eta col
	xi = repmat(xi,meta,1);
	eta = repmat(eta,1,nxi);
elseif mxi > 1 && nxi == 1 && meta == 1 && neta > 1 % xi col, eta row
	xi = repmat(xi',neta,1);
	eta = repmat(eta',1,mxi);
elseif mxi == 1 && nxi > 1 && meta == 1 && neta > 1 % both rows
	xi = repmat(xi,neta,1);
	eta = repmat(eta',1,nxi);
elseif meta == 1 && neta == 1 % eta scalar
	xi = xi(:);
elseif mxi == 1 && nxi == 1 % xi scalar
	eta = eta(:);
end

r = dr*sqrt(xi.^2 + eta.^2); % radius in microns 
tau = taup*xi + tauL*eta; % time tau in s

D1 = 1 + 4*D*tau/w0^2;
D2 = 1 + 4*D*tau/wz^2;

S = exp(-r.^2/w0^2./D1);
G = 1./D1./sqrt(D2);

F = (A-B)*S.*G + B;

if exist('long','var') && long
	% "long" format is for fitting, and we eliminate the first data point,
	% because it is spurious for fitting
	F = F(:);
	F(1) = [];
end

