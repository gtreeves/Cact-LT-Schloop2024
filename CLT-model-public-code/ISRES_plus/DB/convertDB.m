clear
clc
close all


mainloc  = '/Users/prasadbandodkar/Desktop/Dl-Cact local/Results/Dl_Cact3-2021_10_14/';
Files   = extractFileLocations(mainloc,'mat');


load('optstruct.mat');
options     = rmfield(options,'model');
fnames      = allFieldNames(options);
nfnames     = length(fnames);

c = 1;
i = 1;
nfiles = length(Files);
while i<=nfiles
    if mod(i,100) == 0
       fprintf('Running: %i/%i\n',i,nfiles) 
    end

    try
        b = load(Files(i),'options','BestMin','Gm');
        filename = strsplit(Files{i},{'/'});
        filename = char(filename(end));
        dati.('fileLocation')  = Files{i};
        dati.('BestMin')       = b.BestMin;
        dati.('Gm')            = b.Gm;

        for j=1:nfnames 
           val = findField(b,(fnames{j}));
           if  ~isempty(val)
              dati.(fnames{j}) = val; 
           else
              dati.(fnames{j}) = -1;  
           end        
        end 
        dati    = convert3to3(dati,b.options);
        dat(c)  = dati;
        
    catch ME     
        disp(ME)
        c = c-1;
    end
    i = i+1;
    c = c+1;
end


% remove duplicates
BestMin = [dat.BestMin];
[~,i0]  = unique(BestMin);
dat     = dat(i0);


save('db_Dl_Cact3.mat','dat')


function dat = convert1to3(dat,options)

struct2vars(options);
dat.mu          = round((pEvo*lambda)/oneByParents);
dat.liveUpdates = false;

dat.sortPrevParamsByError = false;
dat.useFullNewtonStep     = false;


if ~yeslin && ~yesnew
   nlin       = 0;
   nnewt      = 0;
   nlinPar    = -1;
   nnewtPar   = -1;
   startnew   = -1;
   endnew     = -1;
   startlin   = -1;
   endlin     = -1;
elseif yeslin && yesnew
   nlin       = round((lambda-round(pEvo*lambda))/2);
   nnewt      = nlin;
   nlinPar    = ceil(pPar*lambda);
   nnewtPar   = nlinPar;
   startnew   = 1;
   endnew     = G;
   startlin   = 1;
   endlin     = G; 
elseif yeslin && ~yesnew
   nlin       = round((lambda-round(pEvo*lambda))/2);
   nnewt      = 0;
   nlinPar    = ceil(pPar*lambda);
   nnewtPar   = -1;
   startnew   = -1;
   endnew     = -1;
   startlin   = 1;
   endlin     = G; 
elseif ~yeslin && yesnew
   nlin       = 0;
   nnewt      = round((lambda-round(pEvo*lambda))/2);
   nlinPar    = 0;
   nnewtPar   = ceil(pPar*lambda);
   startnew   = 1;
   endnew     = G;
   startlin   = -1;
   endlin     = -1; 
end
dat.nlin       = nlin;
dat.nnewt      = nnewt;
dat.nlinPar    = nlinPar;
dat.nnewtPar   = nnewtPar;
dat.startnew   = startnew;
dat.endnew     = endnew;
dat.startlin   = startlin;
dat.endlin     = endlin; 


end

function dat = convert2to3(dat,options)

struct2vars(options.evo);
struct2vars(options.plus);

dat.mu          = round((pEvo*lambda)/oneByParents);
dat.liveUpdates = false;

dat.sortPrevParamsByError = false;
dat.useFullNewtonStep     = false;


if ~yeslin && ~yesnew
   nlin       = 0;
   nnewt      = 0;
   nlinPar    = -1;
   nnewtPar   = -1;
elseif yeslin && yesnew
   nlin       = round((lambda-round(pEvo*lambda))/2);
   nnewt      = nlin;
   nlinPar    = ceil(pPar*lambda);
   nnewtPar   = nlinPar; 
elseif yeslin && ~yesnew
   nlin       = round((lambda-round(pEvo*lambda))/2);
   nnewt      = 0;
   nlinPar    = ceil(pPar*lambda);
   nnewtPar   = -1;
elseif ~yeslin && yesnew
   nlin       = 0;
   nnewt      = round((lambda-round(pEvo*lambda))/2);
   nlinPar    = 0;
   nnewtPar   = ceil(pPar*lambda);
end
dat.nlin       = nlin;
dat.nnewt      = nnewt;
dat.nlinPar    = nlinPar;
dat.nnewtPar   = nnewtPar;
dat.startnew   = startnew;
dat.endnew     = endnew;
dat.startlin   = startlin;
dat.endlin     = endlin; 


end

function dat = convert3to3(dat,options)

struct2vars(options.evo);
struct2vars(options.plus);
dat.sortPrevParamsByError = sortByError;
dat.useFullNewtonStep     = true;


end
