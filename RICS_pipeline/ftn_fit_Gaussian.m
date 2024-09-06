function F = ftn_fit_Gaussian(P,xi,dr,long)
%Calculates Gaussian PSF model for fitting to fast-direction ACF data
%
%function F = ftn_fit_Gaussian(P,xi,dr,long)
%
% This function is used by the "lsqcurvefit" function to fit a Gaussian to
% the fast direction of the ACF to give us A (the value of the ACF at xi =
% eta = 0), background B (along the fast direction). We also allow the PSF
% size, w0, to vary slightly for robustness
A = P(1);
B = P(2);
w0 = P(3);

r = dr*xi; % radius in microns 

S = exp(-r.^2/w0^2);

F = (A-B)*S + B;

if exist('long','var') && long
	% "long" format is for fitting, and we eliminate the first data point,
	% because it is spurious for fitting
	F = F(:);
	F(1) = [];
end

