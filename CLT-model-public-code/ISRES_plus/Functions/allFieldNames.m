function output = allFieldNames(s)

    names = fieldnames(s);
    n  = length(names);
    output = {};
    for i=1:n
       name = cell2mat(names(i)); 
       if isstruct(s.(name))
          op = allFieldNames(s.(name));
          output = [output; op];
       else
          output = [output; name]; 
       end
    end
    
    


end