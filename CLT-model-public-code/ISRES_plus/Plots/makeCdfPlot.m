function [h,txt] = makeCdfPlot(dat,opts,stat2plot,index)  

h = [];
txt = [];

[~,nUall]  = findUnique(dat,opts);
if ~isempty(nUall)
    [BestMin,vals]      = getStats(dat,opts,stat2plot,index);
    BestMin             = min(cell2mat(BestMin),[],2);   
    nFiles              = length(BestMin);
    h                   = cdfplot(BestMin);
    h.LineWidth         = 2;
    names               = fieldnames(vals);
    txt                 = cell(length(names),1);
    %set(h,'color','black')
    hold on
    for i=1:length(names)   
       txt{i} = [names{i},' = ',num2str(vals.(names{i}))];
    end
    txt{end+1}  = ['nFiles = ',num2str(nFiles)];
end













   
















end