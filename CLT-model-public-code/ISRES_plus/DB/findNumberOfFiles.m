function T = findNumberOfFiles(dat,opts)
% This function create a table of number of files based on the conditions
% laid out in the structure variable "opts"

    nU   = findUnique(dat,opts,false);
    if ~isempty(nU)
        names  = fieldnames(nU);
        vals   = struct2cell(nU);
        vals   = cellfun(@(x) x(1,1:end),vals,'UniformOutput',false);
        T      = findCombinations(names,vals);
    end

    S       = table2struct(T);
    nComb   = size(T,1);
    nFiles  = [];
    BestMin = [];
    keepi   = true(nComb,1);
    f = waitbar(0, 'Starting');
    for i=1:nComb
        v0   = findInDb(dat,S(i));
        if ~isempty(v0)
            nFiles   = [nFiles;length(v0)];
            BestMin  = [BestMin; mean([dat(v0).BestMin])];
        else
            keepi(i)    = false;
        end
        waitbar(i/nComb, f, sprintf('Progress: %d %%', floor(i/nComb*100)));
    end
    T = T(keepi,:);
    T.('#files') = nFiles;
    T.('BestMinMedian') = BestMin;

    

%     names     = fieldnames(opts);
%     vals      = struct2cell(opts);
% %     i0        = cellfun(@isempty,vals);
% %     i0        = find(i0);
% %     vals0     = findAllValues(dat,names,i0);
% %     vals(i0)  = vals0;
% 
% 
%     T = findCombinations(names,vals);
%     i = 1;
%     while i <= size(T,1)     
%        a    = table2struct(T(i,:)); 
%        v0   = findInDb(dat,a);
%        BestMin(i) = mean([dat(v0).BestMin]);
%        v(i) = length(v0);
%        nU   = findUnique(dat,a,false);
%        if ~isempty(nU)
%             tnames  = fieldnames(nU);
%             tvals   = struct2cell(nU);
%             tvals   = cellfun(@(x) x(1,1:end),tvals,'UniformOutput',false);
% %             for j=1:length(tnames)                
% %                 opts.(tnames{j}) = tvals{j};
% %             end
%             names = [names;tnames];
%             vals  = [vals; tvals];
%             T = findCombinations(names,vals);
%             i = 0;
%        end
%        i = i+1;
%     end
%     
%     v        = reshape(v,[],1);
%     BestMin  = reshape(BestMin,[],1);
%     i0       = (v > 0)';
%     T(~i0,:) = [];
%     
%     T.('#files') = v(i0);
%     T.('BestMinMedian') = BestMin(i0);
%     
%     T = removevars(T,fieldnames(opts));
   
end


