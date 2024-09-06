function [IM,lsminf1,lsminf2,bg,bgsig,istart,iend,Grpidx,tstart] = ...
	ftn_imread(filename,subtr_bg,stabilize,istart,iend)
%Read in microscopy image file, perform background subtr & image stabilization if asked for
%
% function [IM,lsminf1,lsminf2,bg,bgsig,istart,iend,Grpidx,tstart] = ...
%	ftn_imread(filename,subtr_bg,stabilize,istart,iend)
%
% This function takes an lsm or czi file (the full path to the file, and the
% filename should be passed in the character variable "filename"), reads it
% in and reads in the metadata. Then, if asked for, the image data is background
% subtracted with the mode of the first frame/z-slice. Alternatively, one
% can provide a background subtraction value. This function also can
% stabilize the image if it is a time series, if asked for. Finally, the
% function can choose to load-in only certain frames of the time series,
% starting with istart and ending with iend (see more below).
%
% Inputs:
% "filename": character variable with full path to the file, including the
%	filename and extension
% "subtr_bg": can be a logical or numeric variable. If it is false, then no
%	background subtraction will take place. If it is logical true, then
%	background subtraction will take place with the mode of each channel's
%	first frame/z-slice. If passed as a numeric, non-zero, non-unity
%	variable, then this value will be used to subtract the background. If
%	the numeric value is non-scalar, then the number of elements should
%	match the number of channels in the image data, and each value will be
%	used to subtract the background of the corresponding channel.
% "stabilize": can be a logical or numeric variable. If it is false, then no
%	image stabilization will take place. If it is logical true, then
%	stabilization will take place assuming the image stack is a time series
%	and that the objects in the image move together as a group. If passed as
%	a numeric, non-zero, non-unity variable, then this value will be used to
%	determine the minimum frame size after stabilization. 
% "istart","iend": Sometimes the user wants only certain frames to be used.
%	There are two ways to do this:
%	(1) The user could pass "istart" and "iend" to this file as numerical
%	indices. If so, that overrides everything.
%	(2) Altenatively, if the user did not provide "istart" or "iend" to
%	this function, or if the user passed either empty arrays or NaN's to
%	both of these, then we will check for a text file that delimits the the
%	frames that should be used. In that case, the user must save a text
%	file under the filename (including the full path) with a ".txt"
%	extension instead of ".lsm" or ".czi". In this text file, the user can
%	specify "i_start=<M>;" and/or "i_end=<N>;" where "M" is the first and
%	"N" is the last frame the user wants included in the analysis. This
%	text file is only used if BOTH passed "istart" and "iend" are
%	nonsensical.
%
% There is another scenario for "istart" and "iend" that also may come into
% play, which is when the czi file (or other microscopy file) has tiles (or
% "scenes", which are separate image stacks all taken in the single
% experiment and stored together inside of a single czi file). In that
% case, "istart" and "iend" need to both be two-element vectors, with the
% first element being the tile number and the second being the frame number
% within that tile. 
%
% There are three subsets of this scenario:
%	Subset 1: The tiles are each single snapshots of different locations in
%	the specimen. In this case, istart and iend mean nothing and should not
%	be used. We have not tested what would happen in this scenario if
%	istart and iend were attempted.
%
%	Subset 2: The tiles are each timecourses, and the microscope is
%	programmed to take a single snap at a tile before moving on to the next
%	tile. After the final tile is imaged for a single frame, the process
%	repeats for "num_frames" cycles. In this Subset, istart and iend only
%	make sense if they isolate a single tile (so the first element of
%	istart has to be numerically equal to the first element of iend).
%
%	Subset 3:  The tiles are each timecourses, as in Subset 2. However, in
%	this Subset, the microscope is programmed to take a short timecourse of
%	N frames at a tile before moving on to the next tile. After the final
%	tile is imaged for N frames, the process repeats for "num_frames/N"
%	cycles. In this Subset, istart and iend can refer to different tiles,
%	and each frame in the experiment is eventually ordered by the time
%	stamp in which it was actually taken. This Subset is called
%	"RICS_tiles" and is useful for a RICS timecourse experiment.
%
%[IM,lsminf1,lsminf2,bg,bgsig,istart,iend,Grpidx,tstart]
% Outputs:
% IM: numerical array or cell array. If there are no tiles, this is a
%	numerical array. The order of dimensions is XYCZT, so by default this
%	is a 5D array. However, if the final dimension(s) is(are) singleton,
%	then Matlab cuts them off. So if the image stack is a single time point
%	z-stack, this will be a XYCD 4D array. However, if it is a single
%	z-slice time course, it will be a XYCZT 5D array, with a singleton 4th
%	dimension. If there are tiles, this will be an ntiles-by-1 cell array,
%	where "ntiles" is the number of tiles. Each element of this cell array
%	obeys the dimension rules laid out above. EXCEPTION: if the experiment
%	was set up to be a "RICS_tiles" experiment, then all frames are
%	concatenated into on 5D numerical array in the order in which they were
%	physically taken by the microscope (ie, by their time stamp)
% "lsminf1": deprecated, do not use
% "lsminf2": image metadata
% "bg","bgsig": if background was subtracted, these vectors contain the
%	background level or standard deviation for each channel
% "istart,iend": the actual istart and iend used
% "Grpidx": only makes sense in a RICS_tiles scenario. The indices of each
%	grouping of frames from the original tiles/cycles set up (useful
%	because in the RICS_tiles scenario, all frames are concatenated into
%	one large 5D numerical array, so this vector "remembers" the groups
%	whence they came). If not RICS_tiles, this is scalar false.
% "tstart": absolute time stamp of the first frame in IM in the format
%	'yyyy-MM-dd''T''HH:mm:ss.SSS'


if ~exist('istart','var') || ~isnumeric(istart) || any(istart <= 0)
	istart = [];
end
if ~exist('iend','var') || ~isnumeric(iend) || any(iend <= 0)
	iend = [];
end
if isempty(istart) && isempty(iend) && exist([filename(1:end-4),'.txt'],'file')
	s = fileread([filename(1:end-4),'.txt']);
	eval(s)
end



%
% Load image and metadata
%
if strcmp(filename(end-2:end),'lsm')
	lsm =  lsmRead2(filename);
	lsminf1 = lsm.lsm;
	IM = lsm.data;
	clear lsm % to save memory.
	lsminf2 = lsminfo(filename);
	
	%
	% Reshape is the image stack is 5D with order 'XYCZT'
	%
	H = lsminf2.DimensionY;
	W = lsminf2.DimensionX;
	D = lsminf2.DimensionZ;
	T = lsminf2.DimensionTime;
	numch = lsminf2.DimensionChannels;

	IM = reshape(IM,[H,W,numch,D,T]);
elseif strcmp(filename(end-2:end),'czi')
	istart0 = istart; iend0 = iend;
	[IM,lsminf1,lsminf2,istart,iend,RICS_tiles,tstart] = openczi(filename,istart0,iend0);
else
	error('Only lsm or czi files supported')
end
if ~iscell(IM)
	IM = {IM};
end

%
% Subtract background, if asked for.
%
if isfield(lsminf2,'DimensionT')
	T = lsminf2.DimensionT;
elseif isfield(lsminf2,'DimensionTime')
	T = lsminf2.DimensionTime;
else
	T = 1;
end
num_ch = lsminf2.NUMBER_OF_CHANNELS;
if exist('subtr_bg','var') && ~isempty(subtr_bg)
	if isnumeric(subtr_bg)

		if length(subtr_bg) < num_ch
			bg = subtr_bg(1)*ones(1,num_ch);

			for i = 1:length(IM)
				IM1 = IM{i};
				IM1 = imsubtract(IM1,bg(1));
				IM{i} = IM1;

			end
		else
			bg = subtr_bg(1:num_ch);

			for i = 1:length(IM)
				IM1 = IM{i};
				for o = 1:num_ch
					IM1(:,:,o,:,:) = imsubtract(IM1(:,:,o,:,:),bg(o));
				end
				IM{i} = IM1;
			end
		end
		bgsig = NaN;

	elseif islogical(subtr_bg) && subtr_bg

		%
		% We assume the background intensity level (i.e., the true black
		% level) for each channel is equal to the mode of intensities seen
		% in a slice.
		%
		for i = 1:length(IM)
			IM1 = IM{i};
			[bg,bgsig] = ftn_calcbg(IM1(:,:,:,1,1),false);
			for o = 1:num_ch
				IM1(:,:,o,:) = imsubtract(IM1(:,:,o,:,:),bg(o));
			end
			IM{i} = IM1;
		end
	else
		bg = NaN;
		bgsig = NaN;
	end
else
	bg = NaN;
	bgsig = NaN;
end



%
% Stabilize movement in the image time series
%
if exist('stabilize','var') && ((isnumeric(stabilize) && stabilize > 0) || ...
		(islogical(stabilize) && stabilize)) && T > 1

	for i = 1:length(IM)
		IM1 = IM{i};
		IM2 = squeeze(sum(IM1,3));
		a = size(IM2);
		if length(a) > 2 && a(3) > 1

			[dx,dy] = image_movement(IM2); % this is probably fine because we won't saturate
			if isnumeric(stabilize)
				minframesize = stabilize;
			else
				minframesize = [];
			end

			IM1 = image_stabilize(IM1,dx,dy,minframesize);
			IM{i} = IM1;
		end
	end
end

%
% Convert cell array
%
Grpidx = false;
if length(IM) == 1
	IM = IM{1};
elseif RICS_tiles
	numdims = cellfun(@ndims,IM);
	numdims = max(numdims);
	H = cellfun(@(x)size(x,1),IM);
	W = cellfun(@(x)size(x,2),IM);
	T = cellfun(@(x)size(x,numdims),IM);

	h = min(H);
	w = min(W);

	%
	% Resize all tiles to make them the same (cutting off a few pixels on
	% the borders...this happened because of image_stabilize). Also forming
	% the "Grpidx" vector
	%
	Grpidx = zeros(sum(T),1);
	count = 1;
	for i = 1:length(IM)
		IM1 = IM{i};

		switch numdims
			case 3
				IM1 = IM1(1:h,1:w,:);
			case 4
				IM1 = IM1(1:h,1:w,:,:);
			case 5
				IM1 = IM1(1:h,1:w,:,:,:);
			otherwise
				error('check dimensions of your image array')
		end
		IM{i} = IM1;

		Grpidx(count:count+T(i)-1) = i;
		count = count + T(i);
	end
	IM = cat(numdims,IM{:});
end











