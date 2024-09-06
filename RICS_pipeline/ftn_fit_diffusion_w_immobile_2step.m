function F = ftn_fit_diffusion_w_immobile_2step(P,XI,dr,wz,taup,tauL,long)
%Two component model function for fitting purposes (called by lsqcurvefit)
%
%function F = ftn_fit_diffusion_w_immobile_2step(P,XI,dr,wz,taup,tauL,long)
%
% This function calculates the theoretical ACF for a two-component
% diffusion and binding model. It is used by the "lsqcurvefit" function to
% fit to the data ACF to give us the diffusivity D, the fraction immoble
% (phi), A (the value of the ACF at xi = eta = 0), and background B. We
% also allow the PSF waist size, w0, to vary slightly for robustness.
% Although usually A and w0 are already determined by the first step in the
% two-step fitting, and B is tightly constrained around zero. Furthermore,
% D is often assumed, based on other data. So, effectively, the main
% adjustable parameter is phi.




%
% Unpacking model parameters
% 
A = P(1);
B = P(2);
D = P(3);
phi = P(4);
w0 = P(5);


%
% Check the size of xi and eta, to put them into 7
% 
xi = XI{1}; eta = XI{2}; 

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

%
% Derive the ACF
%
if ~exist('long','var')
	long = false;
end
G1cpt = ftn_fit_1comp_diffusion([1 0 D w0],XI,dr,wz,taup,tauL,long);
PSF = ftn_fit_Gaussian([1 0 w0],sqrt(xi.^2 + eta.^2),dr,long);

F = (A-B)*((1-phi)*G1cpt + phi*PSF) + B;
