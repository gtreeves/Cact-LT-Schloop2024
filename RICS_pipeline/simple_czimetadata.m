function [H,W,T,num_ch,D,dr,tauf,timestamp0,omeXML,omeMeta,reader,...
	TimeStamps,PositionXYZ,numSeries,FramesPerIteration,numIterations] = simple_czimetadata(filename)
%Extacts simple OME metadata from czi file using bioformats
%
%function [H,W,T,num_ch,D,dr,tauf,timestamp0,omeXML,omeMeta,reader,...
%	TimeStamps,PositionXYZ,numSeries,FramesPerIteration,numIterations] = simple_czimetadata(filename)
%
% This function reads simple metadata from a czi file without having to
% read the whole image. This is not a complete set of metadata, but it is
% enough for some purposes.
%
% Metadata extracted (outputs):
% H,W,T,num_ch,D: the dimensions of the image stack (height, width, number
%	of frames, number of channels, number of z-slices)
% dr: pixel size (meters per pixel)
% tauf: frame time
% timestamp0: absolute time stamp of first frame, in the format of
%	'yyyy-MM-dd''T''HH:mm:ss.SSS' (maybe no ".SSS")
% omeXML: metadata stored in text-based XML format
% omeMeta: metadata in java hashtable format (I think)
% reader: java object that metadata were pulled from. Can be used to pull
%	more information, such as a more complete set of metadata (see
%	"global_czimetadata") or image data
% Timestamps: a list of absolute time stamps of each frame in the image
%	stack(s). Plural if there are multiple tiles/scenes
% PositionXYZ: the relative coordinates of the stage for each frame. Only
%	useful if there are multiple tiles/scenes. They are "relative"
%	coordinates because the stage could have been "zeroed" at any location
% numSeries: the number of tiles/scenes in a series
% FramesPerIteration: if there is more than one tile/scene in the file, how
%	many frames are taken in one tile before moving to the next tile?
% numIterations: if there is more than one tile/scene in the file, how
%	many times do we loop through the tiles?
%
% Other metadata that I am not extracting right now, but you could, is your
% dimension order, the laser excitation wavelength, which PMT detector
% you're using, pinhole size, channel index (which might have some
% information about the track?)

reader = bfGetReader(filename);
omeMeta = reader.getMetadataStore();
omeXML = char(omeMeta.dumpXML());

% =========================================================================
% From the XML metadata, we can extract name/value pairs
% =========================================================================

%
% Start with acquisition date
%
key1 = '<AcquisitionDate>';
keylength1 = length(key1);
k1 = strfind(omeXML,key1);
key2 = '</AcquisitionDate>';
k2 = strfind(omeXML,key2);
timestamp1 =  omeXML(k1+keylength1:k2(1)-1);
try
	timestamp0 = datetime(timestamp1,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
catch
	timestamp0 = datetime(timestamp1,'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
end

%
% Count the number of "series", or "scenes", or "tiles" by counting the
% number of instances of "<Image ID="
%
kImageStart = strfind(omeXML,'<Image ID=');
kImageEnd = strfind(omeXML,'</Image>');
numSeries = length(kImageStart);

TimeStamps = cell(numSeries,1);
PositionXYZ = zeros(numSeries,3);
for ii = 1:numSeries

	%
	% Extract chunk of the XML text corresponding to series ii
	%
	omeXML1 = omeXML(kImageStart(ii):kImageEnd(ii));

	%
	% Extract some parameters directly from the omeMeta variable...for now,
	% we will ignore that these could change with each series/scene/tile.
	% We don't want to have to deail with that right now, and this code is
	% the same as it was before (when the five lines below were outside of
	% a for loop like this, since we were ignoring series at the time), but
	% at least here inside the for loop, we have taken the next step. If we
	% do need separate values for each series, we can put these into cell
	% variables.
	%
	W = omeMeta.getPixelsSizeX(ii-1).getValue(); % image width, pixels
	H = omeMeta.getPixelsSizeY(ii-1).getValue(); % image height, pixels
	D = omeMeta.getPixelsSizeZ(ii-1).getValue(); % number of Z slices
	T = omeMeta.getPixelsSizeT(ii-1).getValue(); % number of time points
	num_ch = omeMeta.getPixelsSizeC(ii-1).getValue(); % number of channels
	
	dr = omeMeta.getPixelsPhysicalSizeX(0).value;
	dr = double(dr)*1e-6; % meters per pixel
	
	%
	% Define locations of "<Plane", then extract all "plane" data. In
	% particular, timestamps.
	%
	kplanestart = strfind(omeXML1,'<Plane ') + length('<Plane ');
	kplaneend = strfind(omeXML1(kplanestart(1):end),'/>') + kplanestart(1) - 2;

	timestamps = zeros(D,T,num_ch);
	for i = 1:D*T*num_ch
		planeattributes = omeXML1(kplanestart(i):kplaneend(i));
		C = str2cell(planeattributes,' ');
		for j = 1:length(C)
			C1 = str2cell(C{j},'"');
			try
				eval([C1{1},C1{2},';'])
			catch
				eval([C{j},';'])
			end
		end
		if length(C) > 3
			timestamps(TheZ+1,TheT+1,TheC+1) = DeltaT;
		else % no data
			timestamps(TheZ+1,TheT+1,TheC+1) = NaN;
		end
	end
	PositionXYZ(ii,:) = [PositionX PositionY PositionZ];

	%
	% Convert TimeStamps into an array that has "num_tracks" elements in
	% the third dimension. And remove third dimension elements (ie, entire
	% "stacks") that are duplicated because two or more color channels were
	% on the same track.
	%
	a = diff(squeeze(timestamps(1,1,:)));
	if ~isempty(a)
		v = a == 0;
		timestamps(:,:,v) = [];
	end
	TimeStamps{ii} = timestamps;

	%
	% Now calculate tauf. We do it this roundabout way because if there is
	% more than one tile/scene, and if there are multiple loops of the series,
	% then you can't trust the "TimeIncrement" attribute of "Pixels".
	% Except for apparently sometimes it glitches and dt = 0 (identically),
	% so we'll throw those out
	%
	dt = diff(timestamps);
	m = median(dt(dt > 0),'omitmissing');
	v = abs(dt - m) > 0.1*m; % you have to be within 10% of median
	tauf = mean(dt(~v),'omitmissing');

	%
	% Calculate the number of frames per iteration of a series. 
	%
	k = find(v);
	if ~isempty(k)
		FramesPerIteration = k(1);
	else
		FramesPerIteration = T;
	end
	numIterations = ceil(size(timestamps,2)/FramesPerIteration);

end













