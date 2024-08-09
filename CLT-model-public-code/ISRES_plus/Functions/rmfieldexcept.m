function opts = rmfieldexcept(opts,exceptfields)
% This functions removes all fields of a structure except the fields
% specified in exceptfields.

    names     = fieldnames(opts);
    diffnames = setdiff(names,exceptfields);
    for j=1:length(diffnames)      
        opts = rmfield(opts,diffnames{j});
    end
  
end