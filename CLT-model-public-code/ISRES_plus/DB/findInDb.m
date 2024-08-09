function v = findInDb(dat,opts)

% load('sample.mat')
fnames = fieldnames(opts);
v0     = true;

for i=1:length(fnames)
   val = opts.(fnames{i});
   if ischar(val) 
      i0  = {dat.(fnames{i})}';
   else
      i0  = [dat.(fnames{i})]';
   end
   if isnumeric(val) || islogical(val)
       cond = i0 == val;
   else
       fun  = @(x,temp) contains(x,temp); 
       temp = cellfun(@(x) fun(x,val),i0,'UniformOutput',false);
       cond = cell2mat(temp);
   end
   v0   = v0 & cond;      
end

v = find(v0==1); 
 
end

