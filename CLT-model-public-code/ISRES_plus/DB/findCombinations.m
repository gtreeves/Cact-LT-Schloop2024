function T = findCombinations(names,vals)
% This functions gives all possible combinations of names and vals. For
% example, if names = ['a','b'], and vals = {[1,2,3],[50,500]}, then the
% result would be a table as follows: 
%     a     b 
%     _    ___
%     1     50
%     1    500
%     2     50
%     2    500
%     3     50
%     3    500

    
    % Remove those variables that are of zero length! 
    ncond     = length(names);
    len       = cellfun(@length,vals);
    [len,i0]  = sort(len,'descend');
    if ~isempty(find(len==0, 1))
       ilen0      = find(len==0, 1);
       len(ilen0) = [];
       i0(ilen0)  = [];
       ncond      = length(len);
    end
    
    % Figure out how many times to repeat an element in the table and how
    % to repeat it.
    vals      = vals(i0);
    names     = names(i0);
    irep      = cumprod(len,'reverse');
    irep(ncond+1) = 1;

    
    % make table
    for i=1:ncond
        v0  = repelem(vals{i},irep(i+1));
        v0  = reshape(v0,[],1);
        T.(names{i}) = repmat(v0,irep(1)/length(v0),1);
    end
    
    T = struct2table(T);
end
