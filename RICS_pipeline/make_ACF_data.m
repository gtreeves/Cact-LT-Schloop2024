function [Gs_out,Imean_out,Var_out,Gs_t] = make_ACF_data(IM,varargin)
%Calculates autocorrelation function from image (or time-series) IM
%
%function [Gs_out,Imean_out,Var_out,Gs_t] = make_ACF_data(IM,varargin)
%
% This function calculates the autocorrelation function for each frame in a
% time-series stack of H-by-W images found in IM.  There are several
% options to this function, which are passed in varargin.  
%
% Inputs:
%	"IM": the H-by-W-by-T image series.  The third dimension is time.
%
% Optional argument varargin can consist of these things, in this order:
%	* "mask": logical array of the same size as IM that will force the
%		autocorrelation function to be calculated only on the logical true
%		part of the mask.  This can be used to calculate the
%		autocorrelation function for any arbitrarily-shaped region, and
%		perhaps even disjoint regions.  If passed as a scalar logical
%		false, then the mask will not be used and the autocorrelation
%		function for the entire frame will be calculated.  Default, false.   
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "subtravg": Logical of whether to subtract off the average frame.  By
%		"average frame" we mean the pixel-by-pixel average of all frames in
%		the time series.  You average in the t-direction and you get an
%		H-by-W average image.  This is subtracted off of each frame in IM.
%		Then, the average intensity of that average frame is added back so
%		that the value of your function will not be very close to zero.  
%		This value of this procedure is it will get rid of the effect of
%		large-scale structures in your image, such as the presence of a
%		nucleus that has a different brightness than the cytoplasm.  This
%		is according to Digman et al., 2005, Biophys J, 88: L33-L36.
%		Default, true.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "sws": Number of frames on either side of frame "i" that you'd like to
%		include in a moving average.   If you specify this number to be
%		zero, then the average of the entire time series stack is
%		subtracted off, as explained above.  If a number greater than zero,
%		then an average frame in a sliding window centered at point "i"
%		will be subtracted from frame "i", in a similar manner to the
%		average of the whole stack being subtracted.  The total number of
%		frames in the sliding window average is 2*(sliding_window_size) +
%		1. Having a sliding window average to subtract off is crucial for a
%		proper calculation of the ACF. This argument is contingent on
%		subtravg being "true". Otherwise, this argument does nothing.
%		Default, five.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "Grpidx": Vector of group indices. Allows us to track which frames
%		belog to which mask group, and will also be used to average frames
%		together for time course purposes. Default, empty.
%
%		NOTE: In the future, these indices may also be used to determine
%		which frames are grouped together for sws procedure.
%
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% "Gs": The 3D autocorrelation matrix for each timepoint.
% "Imean_out": the average intensity of each frame (taking into account the
%	mask, if provided)
% "Var_out": the variance in intensity of each frame (taking into account
%	the mask, if provided), including shot noise. This is calculated by
%	multiplying the uncorrected Gs(1,1) by Imean for each frame.
%	



%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask = varargin{iArg}; else 
	mask = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	subtravg = varargin{iArg}; else 
	subtravg = true;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	sws = varargin{iArg}; else
	sws = 5;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Grpidx = varargin{iArg}; else
	Grpidx = [];
end%, iArg = iArg + 1;

%
% Consistency check for sliding window size
%
if ~isnumeric(sws) || sws < 0
	error('Sliding window size must be a non-negative number.')
end

[H,W,T,num_ch] = size(IM);

%
% Options
%
RICS_tiles = false;
if length(mask(:)) > 1
	usemask = true;
	
	%
	% Convert mask, if needed, into one that has the right number of
	% frames.
	% 
	if ~isempty(Grpidx) && length(Grpidx) == T
		RICS_tiles = true;
		mask1 = false(H,W,T);
		for i = 1:T
			k = Grpidx(i);
			mask1(:,:,i) = mask(:,:,k);
		end

	elseif size(mask,3) ~= T
		Grpidx = zeros(T,1);
		ngroups = size(mask,3);
		nslices = floor(T/ngroups);
		mask1 = false(H,W,T);
		for i = 1:T
			k = min(ceil(i/nslices),ngroups);
			mask1(:,:,i) = mask(:,:,k);
			Grpidx(i) = k;
		end
	else
		mask1 = mask;
	end
else
	usemask = false;
end

if subtravg
	I_avgframe = mean(IM,3); % need to wait to calc I_avgglobal.
else
	I_avgframe = 0;
end
if zerospad
	Z = zeros(H,W,num_ch); 
	Z1 = Z(:,:,1);
else
	Z = []; 
	Z1 = Z;
end

%
% Loop through each frame and find the autocorrelation function for each.
%
Gs = zeros(floor(H/2),floor(W/2),T);
Imean_out = zeros(T,1);

for i = 1:T
	if subtravg && sws > 0 && ~RICS_tiles
		i1 = max(1,i-sws);
		i2 = min(T,i+sws);
		I_avgframe = squeeze(mean(IM(:,:,i1:i2,:),3));
	elseif RICS_tiles
		idx = Grpidx == Grpidx(i);
		I_avgframe = squeeze(mean(IM(:,:,idx,:),3));

	end
	I = squeeze(IM(:,:,i,:)) - I_avgframe; % need to wait to add I_avgglobal.
	
	%
	% If asked for, use mask to set all non-nuclear to zero.  This also
	% allows us to calculate I_avgglobal and add it back in.
	%
	if usemask && subtravg && sws ~= 0
		mask = mask1(:,:,i);
		if num_ch > 1 % CCF
			I_avgframe1 = I_avgframe(:,:,1);
			I_avgframe2 = I_avgframe(:,:,2);
			I_avgglobal(1) = mean(I_avgframe1(mask)); % global average within nuclear mask
			I_avgglobal(2) = mean(I_avgframe2(mask));
			I1 = I(:,:,1);
			I2 = I(:,:,2);
			I1 = I1 + I_avgglobal(1); % add back on the global average.
			I2 = I2 + I_avgglobal(2);
			I1(~mask) = 0;
			I2(~mask) = 0;
			I(:,:,1) = I1;
			I(:,:,2) = I2;
			Imean1 = mean(I1(mask));
			Imean2 = mean(I2(mask));
			Imean = sqrt(Imean1*Imean2);
		else
			I_avgglobal = mean(I_avgframe(mask)); % global average within nuclear mask
			I = I + I_avgglobal; % add back on the global average.
			I(~mask) = 0; % setting to zero outside of the nucmask
			Imean = mean(I(mask)); % frame-specific average within the nucmask
		end
	elseif ~usemask && subtravg
		if num_ch > 1
			I_avgframe1 = I_avgframe(:,:,1);
			I_avgframe2 = I_avgframe(:,:,2);
			I_avgglobal(1) = mean(I_avgframe1(:)); % global average within nuclear mask
			I_avgglobal(2) = mean(I_avgframe2(:));
			I1 = I(:,:,1);
			I2 = I(:,:,2);
			I1 = I1 + I_avgglobal(1); % add back on the global average.
			I2 = I2 + I_avgglobal(2);
			I(:,:,1) = I1;
			I(:,:,2) = I2;
			Imean1 = mean(I1(:));
			Imean2 = mean(I2(:));
			Imean = sqrt(Imean1*Imean2);

		else
			I_avgglobal = mean(I_avgframe(:));
			I = I + I_avgglobal; % add back on the global average.
			Imean = mean(I(:)); % denominator
		end
	elseif usemask && (~subtravg || sws == 0) % have mask but no SWS?

		mask = mask1(:,:,i);
		if num_ch > 1
			I1 = I(:,:,1);
			I2 = I(:,:,2);
			I1(~mask) = 0;
			I2(~mask) = 0;
			I(:,:,1) = I1;
			I(:,:,2) = I2;
			Imean1 = mean(I1(mask));
			Imean2 = mean(I2(mask));
			Imean = sqrt(Imean1*Imean2);
		else
			I(~mask) = 0; % setting to zero outside of the nucmask
			Imean = mean(I(mask)); % frame-specific average within the nucmask
		end
	else
		Imean = mean(I(:)); % denominator
	end
	
	%
	% Perform the correlation
	%	
	I2 = [I Z; Z Z];
	if num_ch > 1
		J1 = fft2(I2(:,:,1));
		J2 = fft2(I2(:,:,2));
		G1 = ifft2(J2.*conj(J1));

	else
		J = fft2(I2);
		G1 = ifft2(J.*conj(J));
	end
	
	%
	% Normalizations. If there is a mask, the normalization is calculated
	% by transforming the mask.  If not, the normalization is just
	% performed by dividing by H*W.
	%
	if usemask
		%
		% Make the normalization from transforming the nucmask.  This gives
		% us the denominator...how many terms got summed into each pixel of
		% the autocorrel ftn G.
		%
		Im = [mask Z1; Z1 Z1];
		Jm = fft2(Im);
		D = ifft2(Jm.*conj(Jm));
		D = round(D);
		
		
	else
		D = H*W;
	end
	G11 = G1;
	G1 = G11./D; % normalizing
	G11 = G1(1:floor(H/2),1:floor(W/2))/Imean^2 - 1;

	G11(isinf(G11)) = 0; % fixing "Inf" elements.
	G11(isnan(G11)) = 0; % fixing "NaN" elements.
	
	Gs(:,:,i) = G11;
	Imean_out(i) = Imean;
end
Var_out = squeeze(Gs(1,1,:)).*(Imean_out.^2); % variance of each frame, including the shot noise

%
% Average the ACFs from each frame
%
Gs_out = mean(Gs,3);


if num_ch == 1

	%
	% Use parabolic method to fix pixels 1 and 2
	%
	xi1 = 2:4;
	idx = xi1 + 1;
	A = [xi1'.^2 ones(length(idx),1)];
	b = Gs_out(1,idx)';
	x = A \ b;
	a = x(1); c = x(2);

	xi1 = 0:1;
	idx = xi1 + 1;
	Gs_out(1,idx) = a*xi1.^2 + c;

end

%
% Note that Gs is an array with num_frames in the third dimension, and that
% we averaged all frames together to form Gs_out. Here we want to keep some
% of the timecourse data, but we still need to do some amount of averaging.
% So we will bin frames into groups, each with "nslices", and average those
% together.
%
if ~isempty(Grpidx)
	ngroups = max(Grpidx); 
else
	ngroups = 10;
	nslices = floor(T/ngroups); 
	for i = 1:ngroups

		if i < ngroups
			grpidx = nslices*(i-1)+1:nslices*i;
		elseif i == ngroups
			grpidx = nslices*(i-1)+1:T;
		end
		Grpidx(grpidx) = i;
	end
end

[h,w] = size(Gs_out);
Gs_t = zeros(h,w,ngroups);
for i = 1:ngroups
	grpidx = Grpidx == i;
	Gs_t(:,:,i) = squeeze(mean(Gs(:,:,grpidx),3));
end

M = sparse(Grpidx,(1:T),1,ngroups,T); sumM = sum(M,2);
Imean_out = (M*Imean_out)./sumM;
Var_out = (M*Var_out)./sumM;




