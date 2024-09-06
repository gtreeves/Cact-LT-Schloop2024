function F = ftn_fit_immobile_CCF(P,XI,dr,long)
%Calculates a zero diffusion CCF model for fitting to CCF data
%
%function F = ftn_fit_immobile_CCF(P,XI,dr,long)
%
% This function is used by the "lsqcurvefit" function to fit to the
% cross-correlation function to give us the amplitude, A and other
% parameters (but the amplitude is the most important)
A = P(1);
B = P(2);
w0 = P(3);
dx = P(4);
dy = P(5);

xi = XI{1}; eta = XI{2}; 

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

r = sqrt((dr*xi-dx).^2 + (dr*eta-dy).^2); % radius in microns 


F = (A-B)*exp(-r.^2/w0^2) + B;

if exist('long','var') && long
	% "long" format is for fitting, and we eliminate the first 2 data
	% points, because they are spurious for fitting
	F(1,1:2) = [NaN NaN];
	F = F(:);
	F(isnan(F)) = [];
end

