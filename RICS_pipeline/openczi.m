function [IM,lsminf1,lsminf2,istart,iend,RICS_tiles,tstart] = openczi(filename,istart,iend)
%Uses bio-formats to open czi files and store their metadata
%
%function [IM,lsminf1,lsminf2,istart,iend,RICS_tiles,tstart] = openczi(filename,istart,iend)
%
% This function opens czi files using the bio-formats package from ImageJ:
% https://www.openmicroscopy.org/site/support/bio-formats5.5/users/matlab/index.html
%
% Inputs:
% "filename": a string containing the filename, including path, of the czi
%	file you want to open.
% "istart","iend": see "ftn_imread" for more details
%
% Outputs are designed to interface with analyze_xs:
%	"IM": the 5D image data: H-by-W-by-num_ch-by-D-by-T, where H and W are
%		the height and width of your frames, num_ch is the number of channels 
%		in your data, D is the number of z-slices, and T is the number of
%		frames. In other words, XYCZT order of dimensions
%	"lsminf1": deprecated structure that used to contain metadata when read
%		from lsm files. Now is empty, but continues to be a placeholder.
%	"lsminf2": structure containing the metadata. It may appear to be set
%		up in an illogical way, but it is organized to match how the
%		function "lsminfo" spits out metadata for lsm files, to ensure
%		back-compatibility with lsm files.
% "istart,iend": the actual istart and iend used; see "ftn_imread" for more
%	details 
% "RICS_tiles": logical that tells us whether we are in the RICS_tiles
%	scenario. See "ftn_imread" for more details
% "tstart": absolute time stamp of the first frame in IM in the format
%	'yyyy-MM-dd''T''HH:mm:ss.SSS'


%
% Check to see if we have a czi file:
%
if ~strcmp(filename(end-2:end),'czi')
	error('You must only pass czi files to "openczi"')
end

%
% Pull the simple and global czi metadata. This will allow us to access a
% subset of the image rather than having to load-in the whole thing
%
[H,W,T,num_ch,D,dr,tauf,timestamp0,omeXML,omeMeta,reader,...
	TimeStamps,PositionXYZ,numSeries,FramesPerIteration,numIterations] = simple_czimetadata(filename);
cziinf = global_czimetadata(filename);

%
% Parse metadata
%
[~,lsminf1,lsminf2] = organize_cziinf(cziinf,H,W,T,num_ch,D,PositionXYZ);
lsminf2.TimeStamps.AvgStep = tauf;

%
% Class and dimension order
%
fnames = fieldnames(cziinf);
v = find(contains(fnames,'ComponentBitCount'));
clas = ['uint',num2str(cziinf.(fnames{v(1)}))];
DimensionOrder = char(omeMeta.getPixelsDimensionOrder(0).getValue());

% =========================================================================
% Extract image data
% =========================================================================

%
% NOTE: there could be a series of "tiles" (or "scenes") in the same czi
% file. There are two scenarios: 
% (1) If there is no series, then create "IM" as you normally would.
% (2) If there is a series, then IM is a cell variable (exception: scenario
% 2c below), and the dimensions of that cell variable, and what each
% element means, is dependent on whether the scenes have more than one time
% point each, or whether the scenes have multiple iterations (cycles
% through).
% 
% (2a) If each scene/tile in the series has only one time point, OR 
% (2b) if there is only one iteration through the series, then:
% IM has dimensions numSeries-by-1, and each element is a scene/tile.
%
% (2c) HOWEVER, if each scene or tile has multiple time points AND the
% imaging looped through the series more than once, then IM will be one 5D
% numerical array, which is a concatenation of all frames, sorted by time
% stamp. This is done for RICS timecourse imaging ("RICS_tiles").
%
if ~exist('istart','var')
	istart = 1;
elseif ~isnumeric(istart) || isempty(istart) || any(isnan(istart))
	if size(istart,2) <= 1 || numSeries == 1
		istart = 1;
	elseif size(istart,2) == 2
		istart = [1,1];
	end
end
if ~exist('iend','var')
	iend = T;
elseif ~isnumeric(iend) || isempty(iend) || any(isnan(iend))
	if size(iend,2) <= 1 || numSeries == 1
		iend = T;
	elseif size(iend,2) == 2
		iend = [numSeries,T];
	end
end



%
% Preallocate
%
if numSeries == 1 || FramesPerIteration == 1 || T == 1
	IM = cell(numSeries,1);
	RICS_tiles = false;
	tstart = TimeStamps{1}(istart);
else
	IM = cell(numSeries,numIterations);
	RICS_tiles = true;
		
	% Start and end time for given stage in this acquisition
	tstart = TimeStamps{istart(1)}(istart(2));
	tend = TimeStamps{iend(1)}(iend(2));
	if isnan(tend)
		tend = max([TimeStamps{:}]);
	end
end

for ii = 1:numSeries

    reader.setSeries(ii-1);

	if RICS_tiles % then istart and iend should be 1x2 row vecs
		
		%
		% Compare timestamps of current scene to start and end times
		%
		timestamps = TimeStamps{ii};
		if tstart > timestamps(end) || tend < timestamps(1)
			continue
		end
		istart1 = find(timestamps >= tstart,1);
		iend1 = find(timestamps <= tend,1,'last');

		%
		% Group the iterations into different bins and eliminate any iterations
		% that have less than 5 members. This includes deleting the
		% corresponding timestamps
		%
		iteration_bins = repmat((1:numIterations),FramesPerIteration,1);
		iteration_bins(1:istart1-1) = NaN;
		iteration_bins(iend1+1:end) = NaN;
		counts = sum(~isnan(iteration_bins));
		iteration_bins(:,counts <= 5) = NaN;
		timestamps(isnan(iteration_bins)) = NaN;
		TimeStamps{ii} = timestamps;
		iteration_bins(isnan(iteration_bins)) = [];

		%
		% Redefine istart and iend to match the frames that were not thrown
		% out
		%
		istart1 = find(~isnan(timestamps),1);
		iend1 = find(~isnan(timestamps),1,'last');
	else
		istart1 = istart; % not yet addressing what to do if we have 
		% multiple scenes but not "RICS_tiles". So this "else" is only for
		% "regular" right now. It might work for series that aren't
		% RICS_tiles, but not sure
		iend1 = iend;
	end

	%
	% Preallocate
	%
	if strcmp(DimensionOrder(4:5),'ZT')
		IM1 = zeros(H,W,num_ch,D,(iend1-istart1+1),clas);
	elseif strcmp(DimensionOrder(4:5),'TZ')
		IM1 = zeros(H,W,num_ch,(iend1-istart1+1),D,clas);
		IM1 = permute(IM1,[1 2 3 5 4]);
		warning('Trying new dimension order')
	else
		error(['dimension order is ',DimensionOrder])
	end

	%
	% Read-in the image frame-by-frame. See also "script_loadpart" in
	% "Matlab_archive/TAMU/Ablations"
	%
	for iT = istart1:iend1
		for iC = 1:num_ch
			for iZ = 1:D
				iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
				IM1(:,:,iC,iZ,iT-istart1+1) = bfGetPlane(reader, iPlane);
			end
		end
		if iT == istart1+1
		    fprintf('Reading series #%d\n',ii);
		end
		if iT/1000 == round(iT/1000)
			disp(['iT = ',num2str(iT),' out of ',num2str(T)])
		end
	end

	if RICS_tiles %size(IM,2) > 1
		IM11 = cell(1,numIterations); % RICS timecourse tiles
		for i = 1:numIterations
			IM11{i} = IM1(:,:,:,:,iteration_bins == i);
		end
		IM(ii,:) = IM11;
	else
		IM{ii} = IM1; % normal
	end
end

if numSeries == 1 % No need for a cell variable if there's only one element
	IM = IM{1};
elseif RICS_tiles
	IM = IM(:);
	v = cellfun(@isempty,IM);
	IM(v) = [];

	timestamps = sort([TimeStamps{:}])';
	istart = find(timestamps >= tstart,1);
	iend = find(timestamps <= tend,1,'last');
end









