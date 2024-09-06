function [cziinfo,lsminf1,lsminf2] = organize_cziinf(cziinf,H,W,T,num_ch,D,PositionXYZ)
%Organize metadata from cziinf and make it back-compatible with lsm files.
%
%function [cziinfo,lsminf1,lsminf2] = organize_cziinf(cziinf,H,W,T,num_ch,D,PositionXYZ)
%
% This function examines the cziinf metadata and extracts certain metadata
% values to place into lsminf2 for easier access. This is also to ensure
% backwards compatability so that earlier versions of the code, that
% reference these metadata values in lsminf2, could find the right
% metadata.


%
% The metadata are stored in the structure "cziinf". Let's mine the data
% for what we need. 
% 
% NOTE: for now, we are missing number of time points, and time stamps.
% Also, note that some of these fields might not exist if we are looking at
% a z-stack or a bleach box or something else.
%

% Imagesize
cziinfo.DimensionXYZ = [H W D]';
cziinfo.num_ch = num_ch;
cziinfo.DimensionT = T;

%
% Scalings. Number of entries can depend on whether you have a zstack or
% not. Presumably, although this isn't tested, if you have a line-scan,
% maybe you'll only have one entry?
%
fieldnamesscaling = fieldnames(cziinf.Scaling.Distance);
k = strfind(fieldnamesscaling,'Value');
nscalings = sum(~cellfun(@isempty,k));
scalings = zeros(1,nscalings);
for i = 1:nscalings
	if nscalings == 1
		scalings(i) = cziinf.Scaling.Distance.Value;
	else
		scalings(i) = cziinf.Scaling.Distance.(['Value',num2str(i)]);
	end
	% meters per pixel
end
cziinfo.Scalings = scalings;

%
% Laserpower, Detector Gain, Amplifier Offset, Pixeltime, Linetime,
% Frametime (all depend on num_ch).
%
% Laserpower is called "transmission", and it's confusing because the
% number of transmission values given does not always match the number of
% color channels there actually are in the image. So for now we will assume
% if there are num_ch channels, then the first num_ch transmission values
% are actually the ones we need.
%
sLP = 'Transmission';
if isfield(cziinf.Channel,'Gain1')
	sDG = 'Gain';
else
	sDG = 'Voltage';
end
LP = zeros(1,cziinfo.num_ch);
DG = zeros(1,cziinfo.num_ch);
AO = zeros(1,cziinfo.num_ch);
PT = zeros(1,cziinfo.num_ch);
LT = zeros(1,cziinfo.num_ch);
FT = zeros(1,cziinfo.num_ch);
if num_ch > 1
	for i = 1:cziinfo.num_ch
		LP(i) = cziinf.Channel.([sLP,num2str(i)]);
		DG(i) = cziinf.Channel.([sDG,num2str(i)]);
		AO(i) = cziinf.Channel.(['Offset',num2str(i)]);
		PT(i) = cziinf.Channel.LaserScanInfo.(['PixelTime',num2str(i)]);
		
		if isfield(cziinf.Channel.LaserScanInfo,['FrameTime',num2str(i)])
			FT1 = cziinf.Channel.LaserScanInfo.(['FrameTime',num2str(i)]);
			LT(i) = FT1/H*1000;
			FT(i) = FT1;
		else % in the case of a line scan
			LT(i) = cziinf.Channel.LaserScanInfo.(['LineTime',num2str(i)]);
			FT(i) = NaN;
		end
	end
elseif isfield(cziinf,'SizeX1')
	LP = cziinf.Channel.([sLP,'1']);
	DG = cziinf.Channel.([sDG,'1']);
	AO = cziinf.Channel.Offset1;
	PT = cziinf.Channel.LaserScanInfo.PixelTime1;
	if isfield(cziinf.Channel.LaserScanInfo,'FrameTime1')
		FT1 = cziinf.Channel.LaserScanInfo.FrameTime1;
		LT = FT1/H*1000;
		FT = FT1;
	else % in the case of a line scan
		LT = cziinf.Channel.LaserScanInfo.LineTime1;
		FT = NaN;
	end
	
else
	try
		LP = cziinf.Channel.(sLP);
	catch
		LP = cziinf.Channel.([sLP,'1']);
	end
	DG = cziinf.Channel.(sDG);
	AO = cziinf.Channel.Offset;
	PT = cziinf.Channel.LaserScanInfo.PixelTime;
	if isfield(cziinf.Channel.LaserScanInfo,'FrameTime')
		FT1 = cziinf.Channel.LaserScanInfo.FrameTime;
		LT = FT1/H*1000;
		FT = FT1;
	else % in the case of a line scan
		LT = cziinf.Channel.LaserScanInfo.LineTime;
		FT = NaN;
	end
	
end
cziinfo.Laserpower = LP;
cziinfo.Detectorgain = DG;
cziinfo.Amplifieroffset = AO;
cziinfo.Pixeltime = PT*1e6;
cziinfo.Linetime = LT;
cziinfo.Frametime = FT;

cziinfo.PositionXYZ = PositionXYZ;

% ROI stuff
if isfield(cziinf,'Experiment') && isfield(cziinf.Experiment.AcquisitionBlock.AcquisitionModeSetup,'UseRois1')
	cziinfo.UseRois = ...
		cziinf.Experiment.AcquisitionBlock.AcquisitionModeSetup.UseRois1;
	cziinfo.FitFramesizeToRoi = ...
		cziinf.Experiment.AcquisitionBlock.AcquisitionModeSetup.FitFramesizeToRoi1;
else
	cziinfo.UseRois = false; % GTR: changed to "false" from "[]" on 6/17/2022
	cziinfo.FitFramesizeToRoi = true;
end

cziinfo.cziinf = cziinf; % all the metadata in case needed for future reference

%
% Now re-work the structure to make it back-compatible with the original
% formulation of analyze_xs
%
lsminf2.NUMBER_OF_CHANNELS = cziinfo.num_ch;
lsminf2.VOXELSIZES = cziinfo.Scalings;
lsminf2.ScanInfo.USE_ROIS = cziinfo.UseRois;
lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS = ~cziinfo.FitFramesizeToRoi;
lsminf2.DimensionX = cziinfo.DimensionXYZ(2); % double check to see if this is true
lsminf2.DimensionY = cziinfo.DimensionXYZ(1);
lsminf2.DimensionZ = cziinfo.DimensionXYZ(3);
lsminf2.DimensionT = T;
lsminf2.ScanInfo.PIXEL_TIME = num2cell(cziinfo.Pixeltime);
lsminf2.TimeInterval = cziinfo.Frametime(1);
lsminf2.ScanInfo.POWER = num2cell(cziinfo.Laserpower);
lsminf2.ScanInfo.DETECTOR_GAIN = num2cell(cziinfo.Detectorgain);
lsminf2.cziinfo = cziinfo;
if isfield(cziinf.T,'Interval') && isfield(cziinf.T.Interval,'Increment')
	lsminf2.TimeStamps.AvgStep = cziinf.T.Interval.Increment;
elseif isfield(cziinf.T,'Interval') && isfield(cziinf.T.Interval,'Increment1')
	lsminf2.TimeStamps.AvgStep = cziinf.T.Interval.Increment1;
else
	lsminf2.TimeStamps.AvgStep = NaN;	
end

lsminf1 = []; % deprecated placeholder
