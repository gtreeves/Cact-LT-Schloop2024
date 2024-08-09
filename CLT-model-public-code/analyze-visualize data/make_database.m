
clear
clc
close all
% addpath /scratch/user/razeen/Manuplus/ISRES_plus/Functions/

% mat_dir   = '/Users/razeenshaikh/Library/CloudStorage/GoogleDrive-razeen@tamu.edu/My Drive/Research_TAMU/Projects/Dorsal/Dl_Cact_Raz/dsat/11-04-2023/dsat/Results';
mat_dir =     '/Users/razeenshaikh/Library/CloudStorage/GoogleDrive-razeen@tamu.edu/My Drive/Research_TAMU/Projects/Cact_LT/code_clt/two component Gtot uniform/Results';
dbname    = 'two_component_model_simple_normalized_GFPsd_rng2.mat';

%% Run from here

try
   load(['Mats/',dbname],'dat')
catch
   fprintf('%s does not exist yet\n',dbname) 
end

Files = extractFileLocations(mat_dir,'mat',false);


%% Find list of files not in the db

if exist('dat','var')
    for i=1:length(dat)
       i0 = find(strcmp(Files,dat(i).fileLocation)); 
       Files(i0) = [];
    end
    fprintf('# of new files = %i\n',length(Files))
end

if isempty(Files)
   return 
end

%% Go through all files one by one and extract relevant info and add to the 
%  database

tic
clear datnew

nfiles  = length(Files);
b = load(Files(1),'opts');
b.options = rmfield(b.opts,'model');
fnames = allFieldNames(b.options);
%fnames = [fieldnames(evo);fieldnames(plus)];
%fnames  = fieldnames(optArray);
nfnames = length(fnames);
%datnew  = optArray;


c = 1;
i = 1;
datnew = [];
while i<=nfiles
    if mod(i,10) == 0
       fprintf('Running: %i/%i\n',i,nfiles) 
    end
    
    try
        b = load(Files(i),'opts','BestMin','xb','Gtot','Y0','KnucG');
%         if b.BestMin < 0.4
        datnew(c).('fileLocation') = Files{i};
        datnew(c).('BestMin') = b.BestMin;
        datnew(c).('xb') = b.xb;
        datnew(c).('Gtot') = b.Gtot;
        datnew(c).('Y0') = b.Y0;
        datnew(c).('KnucG') = b.KnucG;
        for j=1:nfnames 
           val = findField(b,(fnames{j}));
           if  ~isempty(val)
               datnew(c).(fnames{j}) = val;
           else
              datnew(c).(fnames{j}) = -1;              
           end        
        end
        c = c+1;
%         end
    catch ME     
        disp(ME)
%         c = c-1;
    end
    i = i+1;
    
end


toc


%% save
if exist('dat','var')
    dat = [dat, datnew];
else
    dat = datnew;
end


%% remove duplicates

BestMin = [dat.BestMin];
[~,i0]  = unique(BestMin);
dat     = dat(i0);


%% Save
save(dbname)

