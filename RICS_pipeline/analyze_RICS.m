function Gs = analyze_RICS(filename,data_ch,genotype,w0,wz,varargin)
%Calculates the RICS experimental ACF for total, nuc, and cyt.
%
%function data = analyze_RICS(filename,channels,channelnames,genotype,varargin)
%
% This function first reads in an image file of a time series acquisition
% suitable for RICS analysis.  Next, after image ROI is cropped and the
% image is background subtracted, the ACF is calculated from the data.  (A
% mask may be calculated first, if asked for).
%
% Inputs:
% "filename": the name of the file, including the path to it
% "data_ch": the channel of your image stack where the data (on which RICS
%	should be performed) reside.
% "genotype": a descriptive genotype of the biological sample
% "w0": the e^(-2) xy waist of the Gaussian PSF
% "wz": the axial waist 
%
% Optional argument varargin can consist of these things, in this order:
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
%	* "sliding_window_size": Number of frames on either side of frame "i"
%		that you'd like to include in a moving average.   If you specify
%		this number to be zero, then the average of the entire time series
%		stack is subtracted off, as explained above.  If a number greater
%		than zero, then an average frame in a sliding window centered at
%		point "i" will be subtracted from frame "i", in a similar manner to
%		the average of the whole stack being subtracted.  The total number
%		of frames in the sliding window average is 2*(sliding_window_size)
%		+ 1.  This argument is contingent on subtravg being "true".
%		Otherwise, this argument does nothing.  Default, five.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "mask_ch": numeric input that tells the function which channel has
%		the image data from which to derive the mask.  This can be used to
%		calculate the autocorrelation function for any arbitrarily-shaped
%		region, and perhaps even disjoint regions.  If passed as a scalar
%		numerical zero, or logical false, then the mask will not be used
%		and the autocorrelation function for the entire frame will be
%		calculated. Default, 2.
%		If this is not specified, but you still want to
%		specify other arguments, put empty brackets -- [] -- in place of
%		this argument.
%	* "usemask": input that tells the function whether to use a mask as an
%		ROI to calculate the data-ACF.  This input can take on several
%		different values and classes. If this input is a string that starts
%		with 'c' (insensitive to case), then the mask will be calculated as
%		the cytoplasmic portion of the image, and that will be used to
%		build the ACF. If it is a string that starts with anything else,
%		then the mask will be the nuclei in the image. If it is anything
%		else that is non-empty and non-logical false, then both a
%		cytoplasmic and nuclear mask will be calculated and used to get the
%		nuclear and cytoplasmic intensities. Default, false.
%		If this is not specified, but you still want to
%		specify other arguments, put empty brackets -- [] -- in place of
%		this argument.
%	* "hasbackground": Logical variable that tells whether the edge of the
%		embryo is present in the image, and if so, that we need to
%		calculate the border between embryo and background (outside the
%		embryo). Default, false.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "stabilize": Logical of whether you'd like the time series to
%		stabilize small movements upon loading. Default, true.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%
%
% The function returns the structure 'Gs' with following the fields:
%
% pth: the full path to the file analyzed
% filenameshort: the filename without the path
% filename: the filename including the full path
% genotype: user-supplied description of the genotype of specimen
% side: placeholder for what side of the embryo was imaged
% metadata: a structure containing the metadata that the user
%	doesn't often need to see or interact with (see below for explanation
%	of what goes in there).
% Gstot,Gsnuc,Gscyt,Gscc,Gsnuc_ch: average ACFs for different aspects of
%	the image (tot = total image; nuc = nuclear mask; cyt = cytoplasmic
%	mask; cc = cross correlation [nuclei only]; nuc_ch = nuclear channel
%	[red channel; nuclei only])
% nucmask,cytmask: h-by-w-by-ngroups masks of the nuclei and cytoplasm,
%	where h and w are the height and width of the image stack after
%	stabilization, ngroups is the number of groups in the timecourse,
%	and each group is composed of "nslices" frames that are averaged
%	together to make a robust mask
% totsignal_all,nucsignal_all,cytsignal_all,nuc_chsignal_all: the
%	background-subtracted intensities of each part of the image, averaged
%	over the whole time course
% totsignal,nucsignal,cytsignal,nuc_chsignal: same as above, but these are
%	vectors (one element for each time point)
% totvar_all,nucvar_all,cytvar_all,nuc_chvar_all: same as two above, but
%	variance
% totvar,nucvar,cytvar,nuc_chvar: same as two above, but variance
% Intensity: average image (data channel only)
% Variance: variance frame across the entire time course
% I_nb,N_nb,B_nb,n_nb,ep_nb,V_nb: placeholders for future N&B analysis (not
%	currently being used)
%
%
% The metadata field, which is a structure, has the following fields:
% data_ch, mask_ch: numerical indices about which channel contains the data
%	and which the mask (usually nuclei)
% lsminf1: empty; deprecated
% lsminf2: structure with more detailed metadata
% scalings
% bg
% std_bg
% S_factor: placeholder for future N&B analysis
% subtravg,sliding_window_size,usemask,hasbackground: recording the
%	varargin status
% xbg,ybg: x and y coordinates of curve cutting across the image that
%	separates the biological specimen from background (only makes sense if
%	"hasbackground" is true)
% H: the number of y-pixels (after stabilization)
% W: the number of x-pixels (after stabilization)
% z_depth: the number of z-slices
% num_frames: number of frames
% num_ch: number of channels
% dr: pixel size in microns per pixel
% w0: xy waist of PSF in microns
% wz: z waist of PSF in microns
% taup: pixel dwell time in microseconds
% tauL: line time in milliseconds
% tauf: frame time in seconds


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	subtravg = varargin{iArg}; else 
	subtravg = true;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	sliding_window_size = varargin{iArg}; else
	sliding_window_size = 5;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else 
	mask_ch = 2;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	usemask = varargin{iArg}; else 
	usemask = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	hasbackground = varargin{iArg}; else 
	hasbackground = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	stabilize = varargin{iArg}; else
 	stabilize = 1;
end%, iArg = iArg + 1;

%
% Check type of argument for usemask
%
if ~islogical(usemask)
	usemask = logical(usemask);
end

%
% filename and load in image and metadata
%
if ~ischar(filename)
	error('wrong argument type')
end
v = find(filename == '\' | filename == '/');
v = v(end);
filenameshort = filename(v+1:end);
pth = filename(1:v);

bg1 = true; % subtract background within ftn_imread (calls ftn_calcbg)
[IM,lsminf1,lsminf2,bg_lvl,sig0] = ftn_imread(filename,bg1,stabilize);


%
% Size and scalings
%
[H,W,num_ch,z_depth,num_frames] = size(IM);
H_orig = lsminf2.DimensionY;
if z_depth == 1
	IM = reshape(IM,[H,W,num_ch,num_frames]); % instead of squeeze because maybe only one channel
end
scalings = lsminf2.VOXELSIZES*1e6; % microns/pixel
dr = scalings(1);

%
% Mask ch and data ch override. This works only for the czi files, because
% I don't know if I can trust the metadata of the lsm files. Luckily, all
% of our lsm files have data_ch = 1.
%
if isfield(lsminf2,'cziinfo') && (isempty(data_ch) || isnan(data_ch))
	cziinf = lsminf2.cziinfo.cziinf.Channel;
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
end


%
% Check for the "channels" input
%
if data_ch > num_ch
	fprintf('Number of channels in the z-stack: %d\n',num_ch)
	error(['Make sure your "data_ch" input is included between 1 and ',num2str(num_ch)])
end



%
% Pixel, line, and frame time
%
taup = lsminf2.ScanInfo.PIXEL_TIME{data_ch}*1e-6; % pixel dwell time in s
tauL = lsminf2.TimeInterval/H_orig; % line time in s
tauf = lsminf2.TimeInterval; % total frame time in s


%
% Middle remarks:
%
Gs.pth = pth;
Gs.filenameshort = filenameshort;
Gs.filename = filename;
Gs.genotype = genotype;
Gs.side = ''; % placeholder
Gs.metadata.data_ch = data_ch;
Gs.metadata.mask_ch = mask_ch;
Gs.metadata.lsminf1 = lsminf1;
Gs.metadata.lsminf2 = lsminf2;
Gs.metadata.scalings = scalings;
Gs.metadata.bg = bg_lvl;
Gs.metadata.std_bg = sig0; % std in intensity due to shot noise
Gs.metadata.S_factor = stabilize; % relationship between digital level and photons
Gs.metadata.subtravg = subtravg;
Gs.metadata.sliding_window_size = sliding_window_size;
Gs.metadata.usemask = usemask;
Gs.metadata.hasbackground = hasbackground;
Gs.metadata.xbg = NaN; % placeholder
Gs.metadata.ybg = NaN; % placeholder
Gs.metadata.H = H; 
Gs.metadata.W = W; 
Gs.metadata.z_depth = z_depth;
Gs.metadata.num_frames = num_frames;
Gs.metadata.num_ch = num_ch;
Gs.metadata.dr = dr;
Gs.metadata.w0 = w0;
Gs.metadata.wz = wz;
Gs.metadata.taup = taup*1e6; % pixel time in microseconds
Gs.metadata.tauL = tauL*1000; % line time in ms
Gs.metadata.tauf = tauf; % frame time in s




% =========================================================================
% Now that the pre-processing is done, it is time to do the analysis here
% in the second half of this function.
% =========================================================================


% 
% Find the mask, if asked for.  Note that this part is not well optimized.
% 
if usemask && isnumeric(mask_ch) && mask_ch > 0 && mask_ch <= num_ch
	if hasbackground
		[xbg,ybg,Vbg] = find_edge(IM,Gs,0);
		Gs.metadata.xbg = xbg;
		Gs.metadata.ybg = ybg;
	else		
		Vbg = [];
	end
	[nucmask,cytmask] = find_mask(IM,Gs,Vbg,0);
else
	nucmask = NaN;
	cytmask = NaN;
end

%
% Calculate the autocorrelation function from our data
%
I = double(squeeze(IM(:,:,data_ch,:)));
[Gstot,Imean_tot,Var_tot] = make_ACF_data(I,false,subtravg,sliding_window_size);	
totsignal = mean(Imean_tot); totvar = mean(Var_tot);
if usemask
	[Gsnuc,Imean_nuc,Var_nuc] = make_ACF_data(I,nucmask,subtravg,sliding_window_size);
	[Gscyt,Imean_cyt,Var_cyt] = make_ACF_data(I,cytmask,subtravg,sliding_window_size);
	nucsignal = mean(Imean_nuc); nucvar = mean(Var_nuc);
	cytsignal = mean(Imean_cyt); cytvar = mean(Var_cyt);
	if num_ch > 1 && mask_ch ~= data_ch
		
		%
		% Cross correlation
		%
		I2 = double(permute(IM,[1 2 4 3]));
		Gscc = make_ACF_data(I2,nucmask,subtravg,sliding_window_size);
		
		%
		% Nuclear channel only
		%
		I_nuc = double(squeeze(IM(:,:,mask_ch,:)));
		[Gsnuc_ch,Imean_nuc_ch,Var_nuc_ch] = make_ACF_data(I_nuc,nucmask,subtravg,sliding_window_size);
		nuc_chsignal = mean(Imean_nuc_ch); nuc_chvar = mean(Var_nuc_ch);
	else
		Gscc = NaN; Gsnuc_ch = NaN;
		Imean_nuc_ch = NaN; Var_nuc_ch = NaN;
		nuc_chsignal = NaN; nuc_chvar = NaN;
	end
else
	Gscyt = NaN; Gsnuc = NaN; Gscc = NaN; Gsnuc_ch = NaN;
	Imean_nuc = NaN; Imean_cyt = NaN; Imean_nuc_ch = NaN;
	nucsignal = NaN; cytsignal = NaN; nuc_chsignal = NaN;
	Var_nuc = NaN; Var_cyt = NaN; Var_nuc_ch = NaN;
	nucvar = NaN; cytvar = NaN; nuc_chvar = NaN;
end

%
% Average fluorescence intensity and variance of each pixel in the time
% course (for N+B)
%
I0 = double(squeeze(IM(:,:,data_ch,:)));
I = mean(I0,3);
V = var(I0,[],3);


%
% Closing remarks
%
Gs.Gstot = Gstot;
Gs.Gsnuc = Gsnuc;
Gs.Gscyt = Gscyt;
Gs.Gscc = Gscc;
Gs.Gsnuc_ch = Gsnuc_ch;
Gs.nucmask = nucmask;
Gs.cytmask = cytmask;
Gs.totsignal_all = Imean_tot;
Gs.nucsignal_all = Imean_nuc;
Gs.cytsignal_all = Imean_cyt;
Gs.nuc_chsignal_all = Imean_nuc_ch;
Gs.totsignal = totsignal;
Gs.nucsignal = nucsignal;
Gs.cytsignal = cytsignal;
Gs.nuc_chsignal = nuc_chsignal;
Gs.totvar_all = Var_tot;
Gs.nucvar_all = Var_nuc;
Gs.cytvar_all = Var_cyt;
Gs.nuc_chvar_all = Var_nuc_ch;
Gs.totvar = totvar;
Gs.nucvar = nucvar;
Gs.cytvar = cytvar;
Gs.nuc_chvar = nuc_chvar;
Gs.Intensity = I;
Gs.Variance = V;


save([Gs.filename(1:end-4),'_Gs.mat'],'Gs');











