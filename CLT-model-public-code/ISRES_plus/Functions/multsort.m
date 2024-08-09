function [A, I] = multsort(options,params)
% This function takes one argument - Options and one optional argument
% params. There are 2 modes in which this function operates. 
%       Mode 1: options is cell array with each element being a
%               structure that has in its fields all params. In this mode,
%               params must be specified. The fields are sorted in the
%               order in which params is specified.                
%
%       Mode 2: options is a m-by-n matrix. In this mode, each
%               column is taken as a parameter. The matrix is sorted in the
%               order of columns.
% Author: Prasad Bandodkar
% Date: 24 Sept 2020 


if exist('params','var')
    n       = length(options);
    numpar  = length(params);
    A       = zeros(n,numpar);
    for i=1:n
        for j=1:numpar        
            A(i,j) = options{i}.(params(j));
        end
    end
else
    [n,numpar] = size(options);
    A          = options;
end


[~,i0]  = sort(A(:,1));
A       = A(i0,:);
I       = i0;

% Now begin sorting from second column onwards
for i=2:numpar
   b  = diff(A(:,i-1));
   i0 = find(b);
   i0 = [1;i0+1;n+1];
   
   a  = A(:,i);
   for j=1:length(i0)-1
       rindex       = i0(j):i0(j+1)-1;
       [~,i1]       = sort(a(rindex));
       i1           = i1 + i0(j) - 1;
       I(rindex)    = I(i1);
       
       % update all columns with this sort order
       for k=i:numpar
            A(rindex,k)  = A(i1,k);
       end
   end
end


end

