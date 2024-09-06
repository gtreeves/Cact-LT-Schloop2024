function [cziinf,metadata] = global_czimetadata(id, useallfields)
%load metadata quickly from a czi file without loading the images
%
%function [cziinf,metadata] = global_czimetadata(id, useallfields)
%
% This function has the first few lines of "bfopen" ver 5.5.1. But for some
% reason, in that version of the bioformats package (and even the most
% updated one as of Jul 2022) there is no function that allows you to
% extract all of the metadata without also loading the image. You can
% quickly load *some* of the metadata using "bfGetReader", but not all. So
% here we are extracting all of the metadata and arranging it into the
% matlab structure cziinf.
%
% The first input variable, "id", can be either a filename, or a "metadata
% cell array". If it's a filename, then this function will use the
% bioformats pacakge to read *just* the metadata from the file (so you
% don't have to read-in the entire set of image data). Once the metadata
% are read-in from the file, you get an N-by-2 cell array (the "metadata
% cell array") with the first column being the metadata key and the second
% being the metadata value.
%
% The second part of this function takes the metadata keys and creates a
% recursive structure that contains the metadata in the way that we're used
% to seeing it.
%
% If the input variable "id" is the "metadata cell array" already, then the
% first part of the function is skipped and the only thing this function
% does is create the recursive structure from the cell array.

%
% Prompt for a file if not input and if the input is not a metadata cell
% already
%
if (~exist('id','var')  || ((ischar(id) || isstring(id)) && ~exist(id, 'file'))) || (iscell(id) && size(id,2) ~= 2)
	[file, path] = uigetfile(bfGetFileExtensions, 'Choose a file to open');
	id = [path file];
	if isequal(path, 0) || isequal(file, 0), return; end
end

%
% If the input is a filename and not a metadata cell already, then load the
% metadata
%
if ~iscell(id) || size(id,2) ~= 2
	autoloadBioFormats = 1;
	stitchFiles = 0;
	
	% load the Bio-Formats library into the MATLAB environment
	status = bfCheckJavaPath(autoloadBioFormats);
	assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
		'to the static Java path or add it to the Matlab path.']);
	
	
	% Initialize logging
	bfInitLogging();
	
	% Get the channel filler
	if isjava(id)
		% then assume it's the reader
		r = id;
	else
		% then assume it's a filename
		r = bfGetReader(id, stitchFiles);
	end
		
	globalMetadata = r.getGlobalMetadata();
	r.close();
	
	jentries = globalMetadata.entrySet.toArray;
	nentries = length(jentries);
	
	metadata = cell(nentries,2);
	for i = 1:nentries
		metadata{i,1} = jentries(i).getKey;
		metadata{i,2} = jentries(i).getValue;
	end
	
else
	metadata = id;
	nentries = size(metadata,1);
end



% =========================================================================
% Now that we have all of the metadata extracted from the java hashtable,
% we need to organize it so it makes sense to the human eye. Sometimes the
% file really begins here if the "id" variable is already the "metadata"
% variable.
% =========================================================================
cziinf = struct;
for i = 1:nentries
	a = metadata{i,1}; % fieldname
	if strcmp(a(1:7),'Global ')
		a(1:7) = [];
	end
	
	%
	% It turns out that the only information worth saving is in
	% "Information|Image", so we'll gate on that
	%
	if length(a) >= 18 && strcmp(a(1:18),'Information|Image|') 
		usefield = true;
		a(1:18) = [];
	elseif length(a) >= 28 && strcmp(a(1:28),'Experiment|AcquisitionBlock|') 
		usefield = true;
		a(1:11) = [];
	elseif strcmp(a(1:8),'Scaling|') 
		usefield = true;
	elseif exist('useallfields','var') && useallfields
		usefield = true;
	else
		usefield = false;
	end
		
	if usefield
		val = metadata{i,2}; % value...but is contained in a string.
		try
			val_ = str2num(val);
			if ~isempty(val_)
				val = val_;
			end
		catch
		end
		
		%
		% Convert the java entries separated by "|" into ones separated
		% by dots to be consistent with Matlab structure/field notation
		%
		v_divider = a == '|';
		s = a; s(v_divider) = '.';
		v_space = s == ' ' | s == '-' | s == '#' | s == '_'; % removes spaces, dashes, and
		% pound signs and, most recently, underscores
		s(v_space) = [];
		
		
		%
		% The character vbl "s" also has many instances of repeated
		% underscores. The next set of codes removes those
		%
		[Vcheck,scheck] = repeatcheck(s);
		Vremove = false(1,length(s));
		for j = 1:length(scheck)
			if strcmp(scheck(j),'_')
				Vremove(Vcheck{j}(2:end)) = true;
			end
		end
		s(Vremove) = [];
		
		%
		% Now store the value in the field. Sometimes it doesn't work, hence
		% the "try-catch" block. We will proceed assuming you don't need the
		% metadata that doesn't work.
		%
		try
			eval(['cziinf.',s,' = val;']); % It has to be this instead of
			% "cziinf.(s) = val;" because you can't reference nested
			% structures in using dynamic field names. Also, the char vbl
			% "s" has too many characters to be called as a dynamic field
			% name.
			
		catch ME % For some reason, a field value can be a structure and a number? This catch block tries to fix that
			
		end
	end
	1;
	
	
end


%
% Now cziinf should be populated with fields, and eight of those will be
% structures themselves. We will sort those names alphabetically. Also an
% important substructure of Channel is LaserScanInfo.
%
cziinf = orderfields(cziinf);
structnames = {'Track','Session','Scaling','Channel','ObjectiveRef','Y','T','X','MicroscopeRef'};
for i = 1:length(structnames)
	cziinf.(structnames{i}) = orderfields(cziinf.(structnames{i}));
end
cziinf.Channel.LaserScanInfo = orderfields(cziinf.Channel.LaserScanInfo);




