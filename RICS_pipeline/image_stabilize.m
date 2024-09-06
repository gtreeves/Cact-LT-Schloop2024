function IM_out = image_stabilize(IM,dx,dy,minframesize)
%Stabilizes a time series of images given calculated movement
%
%function IM_out = image_stabilize(IM,dx,dy)
%
% This function takes vectors dx and dy, which represent how much the image
% globally moves in the x and y dir (respectively) with repect to the first
% frame in the time series, and translates the subsequent images to remove
% the movement and align with the first image. As a result, the output
% time series of images will have smaller frame sizes.
%
% Note that this function makes no sense in fixed images, and should not be
% used in that case.
%
% Note also that this function assumed the dimension order is XYCZT.
% However, there is a chance that the dimension order is XYT or XYCT
% (leaving out either num_ch or D). 

if ~exist('minframesize','var') || isempty(minframesize)
	minframesize = 0;
end
n = ndims(IM);
if n == 3
	[H,W,num_frames] = size(IM); % XYT
	num_ch = 1;
	D = 1;
	IM = permute(IM,[1 2 4 5 3]);
elseif n == 4
	[H,W,num_ch,num_frames] = size(IM); % XYCT
	D = 1;
	IM = permute(IM,[1 2 3 5 4]);
elseif n == 5
	[H,W,num_ch,D,num_frames] = size(IM); % XYCZT
elseif n == 2
	s = sprintf(['You need to have at least three dimensions for me to take ',...
		'seriously that you have a time course.\n',...
		'(Unless, of course, this is a pCF, and if so, you need to investigate '...
		'"image_stabilize.m" to account for pCFs)\n']);
	error(s)
else
	error('You have too many dimensions (greater than 5). Expected: XYCZT')
end
% Regardless of the original composition of the image, we now have a 5D
% array of size [H W num_ch D num_frames]


%
% make the stabilization 
%
IM_out = IM;
for j = 2:num_frames
	IM_out(:,:,:,:,j) = imtranslate(IM(:,:,:,:,j),-[dx(j-1),dy(j-1)],'FillValues',0);
end

%
% Reduce frame size to take into account the movement. If there is a
% minimum frame size, then we will cut the number of frames.
% 
i1 = cummax(max(1,-round(dy))); i2 = cummin(min(H,H-round(dy)));
j1 = cummax(max(1,-round(dx))); j2 = cummin(min(W,W-round(dx)));
h = i2 - i1 + 1; w = j2 - j1 + 1;
v = find(h < minframesize | w < minframesize);
if ~any(v)
	IM_out = IM_out(i1(end):i2(end),j1(end):j2(end),:,:,:);
else
	lastframe = max(1,v(1)-1); % this leaves the possibility there are very
	% few frames (maybe only one frame). I don't want to have to plan for
	% those contingencies, so we'll leave it the way it is.
	IM_out = IM_out(i1(lastframe):i2(lastframe),j1(lastframe):j2(lastframe),:,:,1:lastframe);
end



%
% Re-permute the IM back to its original dimensions
%
if n ~= 5
	IM_out = squeeze(IM_out);
end



