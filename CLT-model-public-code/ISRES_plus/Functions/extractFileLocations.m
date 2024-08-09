function FileLocations = extractFileLocations(directory,ext,checkSubFolder)
% This functions extracts all ".ext" files from a given parent directory 
%
% Input :
%       "directory": character/string vector contains name of parent 
%                     directory to search for ".ext" files
%       "ext": Searches for files with the extension ".ext"
%       "checksubFolder" (optional): accepts boolean values- true or false
%                       If yes, the function includes subfolders in the
%                       directory. Default is true.
%
% Output :
%       "FileLocations": String column vector that contains file locations

if ~exist('ext','var')
   ext = 'lsm';
end
if ~exist('checkSubFolder','var')
   checkSubFolder = true;
end



% extract extension if user entered ".ext" instead of ext.
c   = strsplit(ext,".");
ext = char(c(end));


% check subfolders if asked for
if checkSubFolder
    a = dir([directory,filesep,'**',filesep,'*.',ext]);
else
    a = dir([directory,filesep,'*.',ext]);
end


a             = struct2table(a);
FileLocations = string(strcat(a.folder,filesep, a.name));

end



