function [bg,sig0,lsminf2,A,B] = ftn_calcbg(filename)
%Calculates the background ("zero level") of image using a Gaussian model
%
%function [bg,sig0,lsminf2,A,B] = ftn_calcbg(filename)
%
% This function is designed to analyze the histogram of a frame of an image
% to get the background level and standard deviation, assuming the
% background (zero-level) pixels are the most prominent and are Gaussian
% distributed.
%
% Input can either be a filename or an image that you're calculating the
% background of.
%
% Outputs:
% "bg": the background level
% "sig0": the standard deviation of the background pixels
% "lsminf2": metadata, if the image is read-in (ie, if "filename" really is
%	a filename instead of the image itself)
% "A","B": amplitude and basal level of the Gaussian fit to the background
%	pixel count distribution

if ischar(filename)
	[IM,~,lsminf2] = openczi(filename);
else
	IM = filename;
	lsminf2 = NaN;
end

A = size(IM);
if length(A) >= 3
	numch = A(3);
else
	numch = 1;
end
bg = zeros(1,numch);
sig0 = zeros(1,numch);
A = zeros(1,numch);
B = zeros(1,numch);
for o = 1:numch

	I0 = double(IM(:,:,o,1)); % assume frame 1 is representative

	%
	% Extracting the slice and finding the "offset"; that is, the
	% background, or "dark current," which is characterized as the most
	% populated peak.
	%
	I0 = I0(:);
	m = max(I0);
	[n,x] = hist(I0(:),0:m);
	[A1,k] = max(n(2:end-1)); % cut off ends because of possibility of saturation
	k = k + 1;

	[~,ksig] = min(abs(n(k+1:end-1)-exp(-0.5)*A1)); ksig = ksig + k; % estimate 1 sigma
	xsig = x(ksig);
	sig = xsig-x(k);
	sig = max(sig,10); % 10 is the typical number; this is to make sure it isn't too small

	%
	% Now fit to a gaussian
	%
	[A1,B1,mu,sig] = fit_gaussian2(x(2:ksig)',n(2:ksig)',sig);
	bg(o) = mu;
	sig0(o) = sig;
	A(o) = A1;
	B(o) = B1;

end

if ischar(filename)
	save([filename(1:end-4),'_bg.mat'],'bg','sig0','lsminf2');
end








