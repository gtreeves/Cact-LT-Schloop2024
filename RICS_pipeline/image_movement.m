function [dx,dy] = image_movement(I)
%Calculates the global movement of the image over a time series
% 
%function [dx,dy] = image_movement(I)
%
% This function take a time series of images and calculates the global
% movement of objects in the image over time. The output is a pair of
% vectors, dx and dy, that represent the total movement of each frame (in
% pixels) for the x and y direction, respectively, with respect to the
% first frame in the time series. 
%
% Note that this function does NOT fix the images. It just calculates the
% movement.


[h,w,num_frames] = size(I);

%
% Loop through each frame
%
dx = zeros(num_frames-1,1);
dy = dx;
for j = 1:num_frames-1
	I1 = double(I(:,:,j));
	I2 = double(I(:,:,j+1));
	
	%
	% open the image
	%
	I1 = imopen(I1,strel('disk',5));
	I2 = imopen(I2,strel('disk',5));
	
	%
	% Gauss blur
	%
	I1 = imgaussfilt(I1,10);
	I2 = imgaussfilt(I2,10);
	
	%
	% Create gradients
	%
	Ix = diff(I1(1:h-1,:),1,2);
	Iy = diff(I1(:,1:w-1));
	It = I2(1:h-1,1:w-1) - I1(1:h-1,1:w-1);
	
	%
	% Calculate best-fit image translation
	%
	A = [Ix(:) Iy(:) -ones((h-1)*(w-1),1)];
	B = -It(:);
	X = A \ B;
	dx(j) = X(1);
	dy(j) = X(2);
end
dx = cumsum(round(dx));
dy = cumsum(round(dy));

