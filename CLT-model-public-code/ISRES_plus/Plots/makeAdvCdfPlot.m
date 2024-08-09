function [h,txt,leg] = makeAdvCdfPlot(dat,opts,stat2plot,index)

%
% Check inputs
%
if ~exist('index','var')
    nstats = length(stat2plot);
    index  = cell(nstats,1);
elseif isempty(index)
    nstats = length(stat2plot);
    index  = cell(nstats,1);
end
rng('shuffle')
leg  = [];


%
% Find if the opts array is a unique set of options. If so, use the cdf plot
% for a single variable
%
nU        = findUnique(dat,opts,false);
if isempty(nU)
   [h,txt]  = makeCdfPlot(dat,opts,stat2plot,index); 
   return
end


%
% If there are varibles in the database that vary outside of the ones
% specified by the opts structure, find those and create combinations.
%
pnames    = fieldnames(nU);
for i=1:length(pnames)
   pvals{i} =  nU.(pnames{i})(1,:);
end
T               = findCombinations(pnames,pvals);
[nvals,nvars]   = size(T);



%
% make cdfplots
%
txt  = [];
leg  = [];  
for j=1:nvals 
    l = [];
    for i=1:nvars
        pname        = pnames{i};
        pval         = T.(pname)(j);
        opts.(pname) = pval;
        l = [l, [pname,' = ',num2str(pval),', ']];         
    end
    l        = l(1:end-2);   
    [h,txt1] = makeCdfPlot(dat,opts,stat2plot,index);
    if ~isempty(txt1)
        leg   = [leg,{l}];
        txt   = [txt,txt1];
    end 
    hold on
end







% %
% % plot vanilla isres
% %
% % load('Mats/db_Dl_Cact-isres.mat','dat')
% % stat2plot           = {'Min'};
% opts.lambda         = 250;
% opts.G              = 75;
% opts.nIslands       = 4;
% % index               = {25};
% 
% % load('Mats/db_Dl_Cact-restart0.1.mat','dat')
% 
% opts.yeslin         = false;
% opts.yesnew         = false; 
% exceptfields        = {'yeslin','yesnew','lambda','G','nIslands'};
% opts                = rmfieldexcept(opts,exceptfields);
% [BestMin,vals]      = getStats(dat,opts,stat2plot,index);
% BestMin             = min(cell2mat(BestMin),[],2);   
% nFiles              = length(BestMin);
% h                   = cdfplot(BestMin);
% h.LineWidth         = 2;
% leg                 = [leg,"ISRES"];
% set(h,'color','black')
% hold off
% names = fieldnames(vals);
% txt1  = cell(length(names),1);
% for i=1:length(names)   
%    txt1{i} = [names{i},' = ',num2str(vals.(names{i}))];
% end
% txt1{end+1}  = ['nFiles = ',num2str(nFiles)];
% txt{end} = txt1;



end

