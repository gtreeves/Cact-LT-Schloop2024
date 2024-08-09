function output = findField(s,field,ind)

if ~exist('ind','var')
   ind = {[]}; 
end

narg = length(ind);
flag = isfield(s,field);
if flag
    output = s.(field);
    opsize = size(output);
    nsize  = length(opsize);
    if narg > nsize
       error('Excess arguments supplied') 
    end
    index    = cell(1,nsize);
    index(:) = {':'};    
    for i=1:narg    
       if ~isempty(ind{i})
           index{i} = ind{i};
       end
    end
 
    output = output(index{:});
    return
else
    names = fieldnames(s);
    for i=1:length(names)
       s0 =  s.(names{i});
       if isstruct(s0)
           output = findField(s0,field,ind);
           if ~isempty(output)
                return;
           end
       end
    end
    output = [];
end
end