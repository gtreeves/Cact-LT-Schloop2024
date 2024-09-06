function [xbg,ybg,Vbg] = find_edge(IM,data,varargin)
%
%
%function [xbg,ybg,Vbg] = find_edge(IM,data,varargin)
%
% This function attempts to find the edge of the biological specimen in the
% case in which the boundary of the specimen is present in the image and in
% which the boundary cuts fully across the image (in other words, the
% specimen is not fully contained in the image).
%
% This function assumes the image has two channels, a data channel and a
% mask channel. The input, "data", contains metadata that tells us which
% channel is which. If there is only one channel, then as long as the
% metadata say that data_ch = 1 and mask_ch = 1, then that also will work.
%
% Optional argument varargin can consist of these things, in this order:
%	* "bg_lvl": The background level in the image. That is, what would the
%		brightness be if there were no embryo in the image? Sometimes, this
%		could be taken from the image itself, but not always. Furthermore,
%		it is not easy to make this an automated step in analyzing the
%		image. That is why this is an input. This background level will be
%		subtracted off of the image intensity to get more accurate absolute
%		readings. Default, zero.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "h": The cutoff for thresholding the border of a nucleus or embryo.
%		Default, 0.25.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs: 
% xbg,ybg: the x and y coordinates of the dividing line between
%	specimen and background (not-specimen)
% Vbg: a logical mask representing the background part of the image



%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	bg_lvl = varargin{iArg}; else 
	bg_lvl = 0;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	h = varargin{iArg}; else
	h = 0.20;
end%, iArg = iArg + 1;

data_ch = data.metadata.data_ch;
mask_ch = data.metadata.mask_ch;
H = data.metadata.H;
W = data.metadata.W;
scalings = data.metadata.scalings;

if bg_lvl > 0
	IM = imsubtract(IM,bg_lvl);
end
I0 = sum(double(IM(:,:,data_ch,:)),4); % sum of intensities
I00 = sum(double(IM(:,:,mask_ch,:)),4); % sum of intensities


% =========================================================================
% First make an averaged frame
% =========================================================================

%
% Pre-processing the image: erosion, dilation, gaussian filtering, then
% combining the the two channels
%
d_erode = 1; % erode by this many microns
d_dilate = 0.5; % dilate by this many microns
see = strel('disk',round(d_erode/scalings(1)));
sed = strel('disk',round(d_dilate/scalings(1)));

I1 = imerode(I0,see);
I1 = imdilate(I1,sed);
I11 = imgaussfilt(I1,round(0.5/scalings(1)));
I01 = imerode(I00,see);
I01 = imdilate(I01,sed);
I011 = imgaussfilt(I01,round(0.5/scalings(1)));

Itot = mDU(mDU(I11) + mDU(I011)); % combining the two channels..."mDU" is a
% function that zero-one normalizes

%
% We have to decide whether the embryo has an edge here, and where that
% edge is. We divide the summed image into several (16) segments
% horizontally and vertically. Then we will see if there is some hard
% boundary that clearly divides the image into bright and dark.
%
% Because we are doing this on a total (averaged) frame, so we are
% implicitly assuming the embryo does not move. Nuclei within the embryo
% might move, (meaning you need to take frame into account when doing
% nuclear vs cytoplasm), but here it's just once for the whole time series.
%
Vbg = false(H,W);
Ibg = mDU(double(imopen(Itot,strel('disk',round(1/scalings(1))))));

ngroups = 16; % how many groups do you want in the image?
nslicesx = floor(W/ngroups); % how many lines are avg'd together in each group?
nslicesy = floor(H/ngroups); % how many lines are avg'd together in each group?
x_out = zeros(ngroups,4);
y_out = zeros(ngroups,4);
for i = 1:ngroups
	
	%
	% Averaging a set of slices together.
	%
	if i < ngroups
		grpidxx = nslicesx*(i-1)+1:nslicesx*i;
		grpidxy = nslicesy*(i-1)+1:nslicesy*i;
	elseif i == ngroups
		grpidxx = nslicesx*(i-1)+1:W;
		grpidxy = nslicesy*(i-1)+1:H;
	end
	x_out(i,3:4) = round(mean(grpidxx));
	y_out(i,1:2) = round(mean(grpidxy));
	
	%
	% Take horizontal slices
	%
	Ih = squeeze(mean(Ibg(grpidxy,:),1));
	k1 = find(Ih(1:end-1) < h & Ih(2:end) >= h); % low then high
	if ~isempty(k1)
		k1 = k1(1);
		x_out(i,1) = interp1(Ih(k1:k1+1),k1:k1+1,h);
	else
		x_out(i,1) = NaN;
	end
	k1 = find(Ih(1:end-1) >= h & Ih(2:end) < h); % high then low
	if ~isempty(k1)
		k1 = k1(end);
		x_out(i,2) = interp1(Ih(k1:k1+1),k1:k1+1,h);
	else
		x_out(i,2) = NaN;
	end
	
	
	%
	% Take vertical slices
	%
	Iv = squeeze(mean(Ibg(:,grpidxx),2));
	k1 = find(Iv(1:end-1) < h & Iv(2:end) >= h);  % low then high
	if ~isempty(k1)
		k1 = k1(1);
		y_out(i,3) = interp1(Iv(k1:k1+1),k1:k1+1,h);
	else
		y_out(i,3) = NaN;
	end
	k1 = find(Iv(1:end-1) >= h & Iv(2:end) < h); % high then low
	if ~isempty(k1)
		k1 = k1(end);
		y_out(i,4) = interp1(Iv(k1:k1+1),k1:k1+1,h);
	else
		y_out(i,4) = NaN;
	end
	
end

%
% Now we will look at scores of how likely the embryo is to have a
% background. After we get a best-fit line and split the image in two, we
% will have an R2 value, mean of one half, mean of the other half, and std
% of each half. Also, m,b for each will be stored.
%
xbg = NaN; ybg = NaN;
for i = 1:4
	x = x_out(:,i); y = y_out(:,i);
	v = ~isnan(x) & ~isnan(y); x = x(v); y = y(v);
	if length(x) > 1
		
		%
		% Demarcate the corner of the image that we think might be
		% background.
		%
		[V,xp,yp] = connectdots_line(x,y,H,W);
		if isnan(V)
			continue
		end
		
		a{1} = Itot(V); b(1) = mean(a{1});
		a{2} = Itot(~V); b(2) = mean(a{2});
		[~,imin] = min(b);
		
		%
		% Double-check to see if the region that we think is bg really is
		% background (below threshold).
		%
		Y = prctile(a{imin},97);
		if Y < h
			xbg = xp;
			ybg = yp;
			if imin == 1
				Vbg = V;
			else
				Vbg = ~V;
			end
			break
		end
	end
end

