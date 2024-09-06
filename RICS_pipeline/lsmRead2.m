function s = lsmRead2(varargin)
%Ftn to read in image z-stacks from lsm fmt, from tiffread ver 2.5
%
%function s = lsmRead2(varargin)
%
% This function differs from "lsmReadDU" because it is capable of
% constructing an output with an arbitrary number of color channels.
%
% NOTE: commented out the "|| width == 128" condition to end the data
% collection. This was originally to prevent thumbnails from trying to be
% combined with the original image. However, some of my images now do have
% a width of 128, but it looks like the 710 doesn't put thumbnails. I also
% allowed for a "info.lsm.DimensionTime" factor in the preallocation of the
% "info.data" field.
%
% Reads 8,16,32 bits uncompressed grayscale and color tiff stacks out of
% the lsm format from Zeiss.  The entire infoF standard is not supported,
% but you may extend it by editing this file. If you do so, it would be
% nice to return your modification to F. Nedelec, so that it can be
% included in future releases.
%
% The function can be called with a file name in the current directory, or
% without argument, in which case it pops up a file openning dialog to
% allow manual selection of the file. If the stack contains multiple
% images, loading can be restricted by specifying the first and last images
% to read, or just one image to read.
%
% At return, "s" is a structure containing the different images with some
% additional information. The images themselves are stored in the field
% .data, a cell array.  Each element of the cell array is a z-slice image,
% and is stored as a 3D array.  Note that the identity of RGB may be
% switched around. The pixel values are returned in the native (integer)
% format, and must be converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.lsm');
% imshow(im.data{1});
%
% Optional argument varargin can consist of the following things:
%	(1) "filename": string containing the filename (including relative
%		path) to the lsm file.  If not specified, a dialog box will appear
%		for user to browse to the file.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% "s": structure containing the information about the images, as well as
% the images themselves in the field ".data": 
%	"data": field of "s", which is a cell array.  Each element of the cell
%		array is a z-slice of the image in the lsm file.
%
% Greg Reeves, Caltech, Feb 2008.  Adapted to read in the whole z-stack of
% lsm files, from:
%
% Francois Nedelec, EMBL, Copyright 1999-2007.
% rewriten July 7th, 2004 at Woods Hole during the physiology course.
% last modified December 21, 2007.
% With contributions from:
%   Kendra Burbank for the waitbar
%   Hidenao Iwai for the code to read floating point images,
%   Stephen Lang to be more compliant with PlanarConfiguration
%   Jan-Ulrich Kreft for Zeiss LSM support
%   Elias Beauchanp and David Kolin for additional Metamorph support
%
% Please, help us improve this software: send us feedback/bugs/suggestions
% This software is provided at no cost by a public research institution.
%
% Francois Nedelec
% nedelec (at) embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg;
%	Germany
% http://www.embl.org
% http://www.cytosim.org

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	filename = varargin{iArg}; else
	% if this isn't specified, the user browses to the file:
	[filename,pathname] = ...
		uigetfile('*.tif;*.stk;*.lsm','select image file');
	filename = [pathname,filename];
end
% end, iArg = iArg + 1;

%
% Preamble.  The structure info is returned to the user, while info is not.
% so tags usefull to the user should be stored as fields in info, while
% those used only internally can be stored in info.
%
% info.data = [];
info.SampleFormat = 1; % setting default values
info.SamplesPerPixel = 1;
info.BOS = 'ieee-le'; % byte order string

info.file = fopen(filename,'r','l');
if info.file == -1
	error(['file <',filename,'> not found.']);
end

%
% Reading in header, including byte order (II = little endian, MM = big
% endian), whether we have tif images (id must be 42), and then the byte
% offset for the first image file directory (IFD)
%
byte_order = fread(info.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
	info.BOS = 'ieee-le'; % normal PC format
elseif ( strcmp(byte_order','MM') )
	info.BOS = 'ieee-be';
else
	error('This is not a tiff file (no MM or II).');
end

tiff_id = fread(info.file,1,'uint16', info.BOS);
if (tiff_id ~= 42) % if it's 42, then it's infoF fmt.
	% Update by GTR on 7/28/2020: changed from error to warning
    warning('This is not a tiff file (missing 42).');
end

info.img_pos = fread(info.file,1,'uint32',info.BOS); % byte offset, 1st IFD
info.filename = fullfile(pwd,filename);
imageCount = 1;

while info.img_pos ~= 0

	%
	% Reading in all entries in this IFD, then the next IFD address.
	%
	info = readIFDentries(info);
    info.img_pos = fread(info.file, 1,'uint32',info.BOS);
	
	%
	% Preallocate the 4D matrix "info.data", which will contain our image
	% stack.
	%
	if ~isfield(info,'data')
		if info.bits == 8
			info.data = zeros(info.height,info.width,...
				info.SamplesPerPixel,info.lsm.DimensionZ,'uint8');
		elseif info.bits == 16
			% EDIT by GTR on 2020-07-14: removed "info.lsm.DimensionTime"
% 			info.data = zeros(info.height,info.width,...
% 				info.SamplesPerPixel,info.lsm.DimensionZ*info.lsm.DimensionTime,'uint16');
			info.data = zeros(info.height,info.width,...
				info.SamplesPerPixel,info.lsm.DimensionZ,'uint16');
		else
			disp('Format other than uint8 or uint16?')
		end
		[imageH,imageW,n_colors,imageD] = size(info.data);
	end

	info.StripCnt = 1;
	Image = zeros(imageH,imageW,n_colors,class(info.data));
	for i = 1:n_colors
		[info,Image1] = read_plane(info,i);
		if isempty(Image1) || size(Image1,1) ~= imageH
			thumbnail = true;
			break
		else
			thumbnail = false;
		end
		Image(:,:,i) = Image1;
% 		planeCnt = 1; [info,red] = read_plane(info,planeCnt);
% 		planeCnt = 2; [info,green] = read_plane(info,planeCnt);
% 		planeCnt = 3; [info,blue] = read_plane(info,planeCnt);
% 		Image = cat(3,red,green,blue);
	end
	if ~isempty(Image) && ~thumbnail
		info.data(:,:,:,imageCount) = Image;
		% 		info.data = [info.data;{Image}];
		imageCount = imageCount + 1;
	end


end

%
% Cleaning up
%
fclose(info.file);
if exist('waitbar_handle', 'var')
	delete( waitbar_handle );
	clear waitbar_handle;
end
s.filename = filename;
s.data = info.data;
s.lsm = info.lsm;


% --------- subfunction to read the image plane (z-slice) --------
function [info,plane] = read_plane(info,planeCnt)

offset = 0;
width = info.width;
height = info.height;

%
% Returning an empty array if the sample format has zero bits
%
if info.BitsPerSample(planeCnt) == 0% || width == 128
	plane = [];
	return;
end

% fprintf(1,'reading plane %i size %i %i\n',planeCnt,width,height);

%
% Determining the type needed to store the pixel values:
%
switch info.SampleFormat
	case 1
		classname = sprintf('uint%i', info.BitsPerSample(planeCnt));
		% This would be uint8, I think.
	case 2
		classname = sprintf('int%i', info.BitsPerSample(planeCnt));
	case 3
		if info.BitsPerSample(planeCnt) == 32
			classname = 'single';
		else
			classname = 'double';
		end
	otherwise
		error('unsuported infoF sample format %i', info.SampleFormat);
	%
end

%
% Preallocate a matrix to hold the sample data:
%
plane = zeros(width,height,classname);

%
% Reading the strips and concatenating them:
%
line = 1;
while info.StripCnt <= info.StripNumber

	[info,strip] = read_strip(info,offset,width,planeCnt,classname);
	info.StripCnt = info.StripCnt + 1;

	%
	% Copying the strip onto the data
	%
	plane(:,line:(line+size(strip,2)-1)) = strip;

	line = line + size(strip,2);
	if line > height
		break
	end

end

%
% Extract valid part of data if needed
%
if ~all(size(plane) == [width height])
	plane = plane(1:width,1:height);
	disp('Cropping data: found more bytes than needed...') % used to be a
	% "warning", but changed to "disp".
end

%
% transposing the image (otherwise display is rotated in matlab)
%
plane = plane';


% --------- subfunction to read strip --------
function [info,strip] = ...
	read_strip(info,offset,width,planeCnt,classname)

% fprintf(1,'reading strip at position %i\n',...
%	info.StripOffsets(stripCnt)+offset)
stripCnt = info.StripCnt;
StripLength =info.StripByteCounts(stripCnt)./info.BytesPerSample(planeCnt);

% fprintf(1,'reading strip %i\n',stripCnt);
fseek(info.file,info.StripOffsets(stripCnt)+offset,-1);
bytes = fread(info.file,StripLength,classname,info.BOS );

if any(length(bytes) ~= StripLength)
	error('End of file reached unexpectedly.');
end

strip = reshape(bytes,width,StripLength/width);




% --------- subfunction to store the type of data --------
function [nbBytes, matlabType] = convertType(tiffType)

switch (tiffType)
	case 1
		nbBytes=1;
		matlabType='uint8';
	case 2
		nbBytes=1;
		matlabType='uchar';
	case 3
		nbBytes=2;
		matlabType='uint16';
	case 4
		nbBytes=4;
		matlabType='uint32';
	case 5
		nbBytes=8;
		matlabType='uint32';
	case 11
		nbBytes=4;
		matlabType='float32';
	case 12
		nbBytes=8;
		matlabType='float64';
	otherwise
		% Update by GTR on 7/28/2020: changed from error to warning. Also
		% added the fake nbBytes and matlabType variables.
		warning('tiff type %i not supported', tiffType)
		
		nbBytes=NaN;
		matlabType='';
	%
end




% --------- subfunction to read IFD entry --------
function  info = readIFDentries(info)

%
% Move in the file to the beginning of the IFD
%
fseek(info.file, info.img_pos,-1);
% disp(strcat('reading img at pos :',num2str(info.img_pos)));

%
% Reading in the number of IFD entries
%
num_entries = fread(info.file,1,'uint16', info.BOS);
% disp(strcat('num_entries =', num2str(num_entries)));

%
% Read and process each IFD entry
%
for i = 1:num_entries
	
	file_pos  = ftell(info.file); % saving current position
	info.entry_tag = ...
		fread(info.file,1,'uint16',info.BOS); % read entry tag

	entry.tiffType = fread(info.file,1,'uint16',info.BOS);
	entry.cnt = fread(info.file,1,'uint32',info.BOS);
	% disp(['tiffType =',num2str(entry.tiffType),...
	%	',cnt = ',num2str(entry.cnt)])

	[entry.nbBytes,entry.matlabType] = convertType(entry.tiffType);
	
	if entry.nbBytes*entry.cnt > 4
		% next field contains an offset:
		offset = fread(info.file,1,'uint32',info.BOS);
		% disp(strcat('offset = ',num2str(offset)));
		fseek(info.file,offset,-1);
	end

	if info.entry_tag == 34412
		entry.val = readLSMinfo(info);
	else
		if entry.tiffType == 5
			entry.val=fread(info.file,2*entry.cnt,entry.matlabType,info.BOS);
		else
			entry.val = fread(info.file,entry.cnt,entry.matlabType,info.BOS);
		end
	end

	if entry.tiffType == 2;
		entry.val = char(entry.val');
	end

	% disp(strcat('reading entry <',num2str(info.entry_tag),'>'));

	%
	% Depending on the value of the entry tag (stored in info.entry_tag),
	% the structure "entry", which was read in from the lsm file using the
	% subfunction "readIFDentry", will mean different things.  Here, we
	% unpack those different meanings.
	%
	switch info.entry_tag
		case 254
			info.NewSubfiletype = entry.val;
		case 256 % image width - number of columns
			info.width = entry.val;
		case 257 % image height - number of rows
			info.height = entry.val;
			info.ImageLength = entry.val;
		case 258 % BitsPerSample per sample
			info.BitsPerSample = entry.val;
			info.BytesPerSample = info.BitsPerSample / 8;
			info.bits = info.BitsPerSample(1);
			% fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
		case 259 % compression
			if entry.val ~= 1
				error('Compression format not supported.');
			end
		case 262 % photometric interpretation
			info.PhotometricInterpretation = entry.val;
			if ( info.PhotometricInterpretation == 3 )
% 				disp('Ignoring infoF look-up table'); % used to be a 
				% "warning", but changed to a "disp".
			end
		case 269
			info.document_name = entry.val;
		case 270 % comments:
			info.info = entry.val;
		case 271
			info.make = entry.val;
		case 273 % strip offset
			info.StripOffsets = entry.val;
			info.StripNumber = entry.cnt;
			% fprintf(1,'StripNumber = %i,...
			%	size(StripOffsets) = %i %i\n',...
			%	info.StripNumber,size(info.StripOffsets));
		case 277 % sample_per pixel
			info.SamplesPerPixel = entry.val;
			% fprintf(1,'Color image: sample_per_pixel=%i\n',...
			%	info.SamplesPerPixel);
		case 278 % rows per strip
			info.RowsPerStrip = entry.val;
		case 279 % strip byte counts - number of bytes in each strip
			% after any compression
			info.StripByteCounts = entry.val;
		case 282 % X resolution
			info.x_resolution = entry.val;
		case 283 % Y resolution
			info.y_resolution = entry.val;
		case 284 % planar configuration describe the order of RGB
			info.PlanarConfiguration = entry.val;
		case 296 % resolution unit
			info.resolution_unit = entry.val;
		case 305 % software
			info.software = entry.val;
		case 306 % datetime
			info.datetime = entry.val;
		case 315
			info.artist = entry.val;
		case 317 % predictor for compression
			if entry.val ~= 1
				error('unsuported predictor value');
			end
		case 320 % color map
			info.cmap = entry.val;
			info.colors = entry.cnt/3;
		case 339
			info.SampleFormat = entry.val;
		case 34412 % Zeiss LSM data
			info.lsm = entry.val;
		otherwise
			fprintf(1,'Ignored infoF entry with tag %i (cnt %i)\n',...
				info.entry_tag,entry.cnt);
			%
	end


	fseek(info.file, file_pos+12, -1); % move to next IFD entry
end


% --------- subfunction to parse LSM info --------
function R = readLSMinfo(info)

%
% Only the first table is read! This provides only very partial info, as
% the offset indicates that additional data is stored in the file.
%

R.MagicNumber = sprintf('0x%x',fread(info.file, 1, 'uint32', info.BOS));
R.StructureSize = fread(info.file, 1, 'int32', info.BOS);
R.DimensionX = fread(info.file, 1, 'int32', info.BOS);
R.DimensionY = fread(info.file, 1, 'int32', info.BOS);
R.DimensionZ = fread(info.file, 1, 'int32', info.BOS);
R.DimensionChannels = fread(info.file, 1, 'int32', info.BOS);
R.DimensionTime = fread(info.file, 1, 'int32', info.BOS);
R.DataType = fread(info.file, 1, 'int32', info.BOS);
R.ThumbnailX = fread(info.file, 1, 'int32', info.BOS);
R.ThumbnailY = fread(info.file, 1, 'int32', info.BOS);
R.VoxelSizeX = fread(info.file, 1, 'float64', info.BOS);
R.VoxelSizeY = fread(info.file, 1, 'float64', info.BOS);
R.VoxelSizeZ = fread(info.file, 1, 'float64', info.BOS);
R.ScanType = fread(info.file, 1, 'uint32', info.BOS);
R.DataType = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetVectorOverlay = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetInputLut = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetOutputLut = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetChannelColors = fread(info.file, 1, 'uint32', info.BOS);
R.TimeInterval = fread(info.file, 1, 'float64', info.BOS);
R.OffsetChannelDataTypes = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetScanInformation = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetKsData = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetTimeStamps = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetEventList = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetRoi = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetBleachRoi = fread(info.file, 1, 'uint32', info.BOS);
R.OffsetNextRecording = fread(info.file, 1, 'uint32', info.BOS);
R.Reserved = fread(info.file, 1, 'uint32', info.BOS);

% 
% % --------- subfunction to read IFD entry --------
% function  [info,info] = readIFDentry(info,info)
% 
% entry.tiffType = fread(info.file, 1, 'uint16', info.BOS);
% entry.cnt      = fread(info.file, 1, 'uint32', info.BOS);
% % disp(['tiffType =',num2str(entry.tiffType),',cnt = ',num2str(entry.cnt)])
% 
% [entry.nbBytes,entry.matlabType] = convertType(entry.tiffType);
% 
% if entry.nbBytes * entry.cnt > 4
% 	% next field contains an offset:
% 	offset = fread(info.file,1,'uint32',info.BOS);
% 	% disp(strcat('offset = ',num2str(offset)));
% 	fseek(info.file,offset,-1);
% end
% 
% if info.entry_tag == 34412
% 	entry.val = readLSMinfo(info);
% else
% 	if entry.tiffType == 5
% 		entry.val = fread(info.file,2*entry.cnt,entry.matlabType,info.BOS);
% 	else
% 		entry.val = fread(info.file,entry.cnt,entry.matlabType,info.BOS);
% 	end
% end
% 
% if entry.tiffType == 2;
% 	entry.val = char(entry.val');
% end
% 
% % disp(strcat('reading entry <',num2str(info.entry_tag),'>'));
% 
% %
% % Depending on the value of the entry tag (stored in info.entry_tag),
% % the structure "entry", which was read in from the lsm file using the
% % subfunction "readIFDentry", will mean different things.  Here, we
% % unpack those different meanings.
% %
% switch info.entry_tag
% 	case 254
% 		info.NewSubfiletype = entry.val;
% 	case 256 % image width - number of columns
% 		info.width = entry.val;
% 	case 257 % image height - number of rows
% 		info.height = entry.val;
% 		info.ImageLength = entry.val;
% 	case 258 % BitsPerSample per sample
% 		info.BitsPerSample = entry.val;
% 		info.BytesPerSample = info.BitsPerSample / 8;
% 		info.bits = info.BitsPerSample(1);
% 		% fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
% 	case 259 % compression
% 		if entry.val ~= 1
% 			error('Compression format not supported.');
% 		end
% 	case 262 % photometric interpretation
% 		info.PhotometricInterpretation = entry.val;
% 		if ( info.PhotometricInterpretation == 3 )
% 			warning('Ignoring infoF look-up table');
% 		end
% 	case 269
% 		info.document_name = entry.val;
% 	case 270 % comments:
% 		info.info = entry.val;
% 	case 271
% 		info.make = entry.val;
% 	case 273 % strip offset
% 		info.StripOffsets = entry.val;
% 		info.StripNumber = entry.cnt;
% 		% fprintf(1,'StripNumber = %i,...
% 		%	size(StripOffsets) = %i %i\n',...
% 		%	info.StripNumber,size(info.StripOffsets));
% 	case 277 % sample_per pixel
% 		info.SamplesPerPixel = entry.val;
% 		% fprintf(1,'Color image: sample_per_pixel=%i\n',...
% 		%	info.SamplesPerPixel);
% 	case 278 % rows per strip
% 		info.RowsPerStrip = entry.val;
% 	case 279 % strip byte counts - number of bytes in each strip
% 		% after any compression
% 		info.StripByteCounts = entry.val;
% 	case 282 % X resolution
% 		info.x_resolution = entry.val;
% 	case 283 % Y resolution
% 		info.y_resolution = entry.val;
% 	case 284 % planar configuration describe the order of RGB
% 		info.PlanarConfiguration = entry.val;
% 	case 296 % resolution unit
% 		info.resolution_unit = entry.val;
% 	case 305 % software
% 		info.software = entry.val;
% 	case 306 % datetime
% 		info.datetime = entry.val;
% 	case 315
% 		info.artist = entry.val;
% 	case 317 % predictor for compression
% 		if entry.val ~= 1
% 			error('unsuported predictor value');
% 		end
% 	case 320 % color map
% 		info.cmap = entry.val;
% 		info.colors = entry.cnt/3;
% 	case 339
% 		info.SampleFormat = entry.val;
% 	case 34412 % Zeiss LSM data
% 		info.lsm = entry.val;
% 	otherwise
% 		fprintf(1,'Ignored infoF entry with tag %i (cnt %i)\n',...
% 			info.entry_tag,entry.cnt);
% 		%
% end
% 
% 
