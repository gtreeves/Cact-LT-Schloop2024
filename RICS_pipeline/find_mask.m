function [Nucmask,Cytmask,Nucstats,Cytstats,Grpidx] = find_mask(IM,data,varargin)
%Segments nuclei in a zoomed-in image of Drosophila blastoderm embryo
%
%function [Nucmask,Cytmask,Nucstats,Cytstats,Grpidx] = find_mask(IM,data,varargin)
%
% This function attempts to segment the nuclei just on a zoomed-in time
% series that is compatible with RICS.
%
% 
% Optional argument varargin can consist of these things, in this order:
%	* "Vbg": Logical array (mask) that tells which pixels are background
%		(pixels are true for background and false for embryo). This
%		argument only makes sense if the embryo edge is in view (ie, if
%		"hasbackground" is true). Default, scalar NaN.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "bg_lvl": The background level in the image. That is, what would the
%		brightness be if there were no embryo in the image? Sometimes, this
%		could be taken from the image itself, but not always. Furthermore,
%		it is not easy to make this an automated step in analyzing the
%		image. That is why this is an input. This background level will be
%		subtracted off of the image intensity to get more accurate absolute
%		readings. Default, zero.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "nslices": the number of frames you'd like to average together to better
%		calculate the mask. This can also be a vector of group indices,
%		which then will be directly passed on to the output "Grpidx". Default,
%		10.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "data_ch": 
%		Default, non-existent.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "mask_ch": 
%		Default, non-existent.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs:
% Nucmask,Cytmask: H-by-W-by-ngroups arrays (where H and W are the height
%	and width of the original image, and ngroups is the number of groups
%	that you decided on, and where a group is a certain number of frames
%	that you chose to average together to make segmentation more robust).
%	Each frame of these arrays is a nuclear or cytoplasmic mask,
%	respectively, of a group of frames from the original image, as
%	specified by the "nslices" varargin input
% Nucstats,Cytstats: the output of regionprops for each mask frame in
%	Nucmask,Cytmask, containing the following fields: 
%		'Meanintensity' (nucmask only),'PixelIdxList'
%	Additionally, these structures also contain 'PerimXYList' not found in
%	regionprops (see function m-file "traceobject" for details)



%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Vbg = varargin{iArg}; else 
	Vbg = false(size(IM,1),size(IM,2));
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	bg_lvl = varargin{iArg}; else 
	bg_lvl = 0;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nslices = varargin{iArg}; else 
	nslices = 10;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	data_ch = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else
	%
end%, iArg = iArg + 1;

%
% Try to make adaptive threshold based on possible bleaching. We will use
% these as starting values and decrease the thresholds if we see bleaching
% happen.
%
h_cyt0 = 0.15;
h_nuc0 = 0.45;

if exist('data_ch','var') && exist('mask_ch','var')
	if isfield(data,'metadata')
		scalings = data.metadata.scalings;
	else % assume it's "lsminf2"
		scalings = data.cziinfo.Scalings*1e6;
	end
elseif isfield(data,'metadata')
	data_ch = data.metadata.data_ch;
	mask_ch = data.metadata.mask_ch;
	scalings = data.metadata.scalings;
else % assume it's "lsminf2" and that the data_ch is green and mask_ch is red
	
	cziinf = data.cziinfo.cziinf.Channel;
	fnames = fieldnames(cziinf);
	v1 = contains(fnames,'Wavelength1');
	v2 = contains(fnames,'Wavelength2');
	if any(v1) && any(v2)
		lambda = [cziinf.Wavelength1 cziinf.Wavelength2];
		data_ch = find(lambda > 450 & lambda < 500);
		mask_ch = find(lambda > 500 & lambda < 600);
	else
		data_ch = 1;
		mask_ch = 2;
	end
	scalings = data.cziinfo.Scalings*1e6;
end

[H,W,num_ch,num_frames] = size(IM);
if any(bg_lvl > 0)
	for o = 1:num_ch
		IM(:,:,o,:) = imsubtract(IM(:,:,o,:),bg_lvl(o));
	end
end



% =========================================================================
% Loop through each frame to get a mask
% =========================================================================

I_data = squeeze(IM(:,:,data_ch,1:num_frames)); % sum of intensities
I_mask = squeeze(IM(:,:,mask_ch,1:num_frames)); % sum of intensities


d_erode = 1; % erode by this many microns
see = strel('disk',round(d_erode/scalings(1)));
d_dilate = 0.5; % dilate by this many microns
sed = strel('disk',round(d_dilate/scalings(1)));

d_nucerode = 0.5; % erode nuclei by this many microns
s_nuc = strel('disk',round(d_nucerode/scalings(1)));
d_cyterode = 0.5; % erode cytoplasm by this many microns
s_cyt = strel('disk',round(d_cyterode/scalings(1)));

if length(nslices) == num_frames % then nslices is really a Group index vector
	Grpidx = nslices;
	Grp_u = unique(Grpidx);
	ngroups = length(Grp_u);
else
	if num_frames > nslices
		ngroups = floor(num_frames/nslices); % how many groups do you want in the image?
	else
		ngroups = 1;
		nslices = num_frames;
	end
	Grpidx = zeros(num_frames,1);
end
Nucmask = false(H,W,ngroups);
Cytmask = false(H,W,ngroups);
Nucstats = cell(ngroups,1);
Cytstats = cell(ngroups,1);
for i = 1:ngroups
	if length(nslices) == num_frames
		
		grpidx = find(Grpidx == Grp_u(i));
	else
		if i < ngroups
			grpidx = nslices*(i-1)+1:nslices*i;
		elseif i == ngroups
			grpidx = nslices*(i-1)+1:num_frames;
		end
		Grpidx(grpidx) = i; % index vector that later will allow us to track
		% which frame belongs to which group
	end
	
	if numel(I_mask) > 1
		I01 = squeeze(sum(double(I_mask(:,:,grpidx)),3));
	else
		I01 = squeeze(sum(double(I_data(:,:,grpidx)),3));
	end
	
	
	%
	% Gaussian blur, erode, dilate, then saturating some pixels.
	%
	Igf = mDU(imgaussfilt(I01,round(0.5/scalings(1))));
	Igf2 = Igf(~Vbg);
	Y = prctile(Igf2,[2 98]); % saturate 2% at each end
	Igf = (Igf - Y(1))/(Y(2) - Y(1));
	Igf(Igf > 1) = 1; Igf(Igf < 0) = 0;

	%
	% Erode, then dilate. The erosion is with a bigger strel than the
	% dilation.
	%
	Ie = imerode(Igf,see);
	Id = imdilate(Ie,sed);
	I1 = Id;

	%
	% Hard threshold to determine the nuclei.
	%
	if i <= 2
		h_cyt = h_cyt0;
		h_nuc = h_nuc0;
	else
		MI2 = mean([Nucstats{i-2}.MeanIntensity]);
		MI1 = mean([Nucstats{i-1}.MeanIntensity]);
		chi = (MI2 - MI1)/MI2; % frac decrease in intensity
		
		h_cyt = h_cyt*(1 - chi);
		h_nuc = h_nuc*(1 - chi);		
	end
	mask_cyt = I1 < h_cyt & ~Vbg; % cytoplasm
	mask_cyt = imerode(mask_cyt,s_cyt); % erode the thresholded cyt mask
	
	mask_nuc = I1 >= h_nuc & ~Vbg; % nuclei
	mask_nuc = imfill(mask_nuc,'holes'); 
	mask_nuc = imdilate(mask_nuc,s_nuc); % dilate the thresholded nuc mask
	
	%
	% Getting outlines of the nuclei and creating final "nucstats"
	% structure
	%
	nucstats = regionprops(mask_nuc,I1,'Meanintensity','PixelIdxList');
	nucstats = traceobject(nucstats,H,W);
	Nucstats{i} = nucstats;
	
	%
	% Getting stats about the cytoplasm
	%
	cytstats = regionprops(mask_cyt,'PixelIdxList');
	cytstats = traceobject(cytstats,H,W);
	Cytstats{i} = cytstats;
		
	Nucmask(:,:,i) = mask_nuc;
	Cytmask(:,:,i) = mask_cyt;
end


