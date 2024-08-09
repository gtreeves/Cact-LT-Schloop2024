function [output,fieldvalues] = getStats(dat,opts,stats,index)
% [output,fieldvalues] = getStat(dat,opts,stats,index)
% [cell, struct] = getStat(struct,struct,char/string,double/int,double/int)
%
% This function will first extract all files from the dat database that
% satisfy the options is opts. In those files, it would extract the stat 
% and store it in the output cell array. It also takes in two additional
% arguments row and col as typically this variable might be a matrix. If
% these aren't specified then all rows and columns are extracted. 
% In addition to the output variable the function also sends as output
% fieldvalues. This refers to the name-value pairs of the dat structure
% that satisfies the opts options. 


%
% Check inputs
%
if ~exist('index','var')
    nstats = length(stats);
    index  = cell(nstats,1);
end

% narg    = length(varargin);
% indices = length(narg);
% for i=1:narg
%    if isempty(varargin{i}) 
%       indices(i) = nan;
%    else
%       indices(i) = varargin{i}; 
%    end
% end


%
% Check if vanilla ISRES is being asked for
%
% if isfield(opts,'yeslin') && isfield(opts,'yesnew')
%    if ~opts.yeslin && ~opts.yesnew
%       exceptfields  = {'yeslin','yesnew','lambda','G','nIslands','fileLocation'};
%       opts         = rmfieldexcept(opts,exceptfields); 
%    end
% end



%
% Find files that satisfy opts and extract stat
%
v       = findInDb(dat,opts);
nfiles  = length(v);
nstats  = length(stats);
output  = cell(nfiles,nstats);
for i=1:nfiles
    if mod(i,100) == 0
        fprintf('File: %i/%i \n',i,nfiles)
    end
   a         = load(dat(v(i)).fileLocation);
   for j=1:nstats
      output{i,j} = findField(a,stats{j},index(j));
   end
end



%
% Get Field values
%
if ~isempty(output)
    dat     = dat(v);
    dat     = rmfield(dat,'Gm');
    dat     = rmfield(dat,'fileLocation');
    dat     = rmfield(dat,'BestMin');
    if isfield(dat,'filename')
       dat  = rmfield(dat,'filename'); 
    end
    fnames  = fieldnames(dat);
    for i=1:length(fnames)      
       vals = [dat.(fnames{i})];
       if length(unique(vals))~=1 && ~ischar(vals)
          disp('Files dont have uniq algorithm params') 
       end
    end   
    fieldvalues = dat(1);  
else
    fieldvalues = [];
end
%}

end


