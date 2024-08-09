function d = getDate(filename,str)
% This function takes as input a character array that is usually the name
% of the file or the location of the file including its name. It takes an
% optional argument called str which is a boolean. If str is true, a date
% string is returned and if str is false, a number is returned that can
% later be transformed into a date string. Default value of str is true.

if ~exist('str','var')
    str = true;
end

d = strsplit(filename,{'/','\','.','_'});
d = d(end-6:end-1);
d = join(d,'-');
d = datenum(d,'yyyy-mm-dd-HH-MM-SS');

if str
   d = datestr(d); 
end

end
