function [nU,nUall] = findUnique(dat,opts,showdisp)
% This function use the database 'dat' and the options 'opts' used in the 
% analysis to find if other parameters the in the options array have
% different values. The output is a structure that captures in the first
% row the value of the variables and in the second row the number of files
% with that value.
    
    if ~exist('showdisp','var')
       showdisp = true; 
    end

    nU          = [];
    nUall       = [];
    dbnames     = fieldnames(dat);
     i0          = cell2mat(cellfun(@(x) contains(x,'fileLocation') || contains(x,'Gm') || ...
         contains(x,'BestMin'), dbnames,'UniformOutput',false)); 
    dbnames(i0) = [];        
    v           = findInDb(dat,opts);
    nFiles      = length(v);
    dat         = dat(v);
    
    
    if ~isempty(dat)
        for i=1:length(dbnames)
           mat  = {dat(1:end).(dbnames{i})};
           mat  = cellfun(@(x) (x(~isnan(x))),mat,'UniformOutput',false);
           if ischar(mat{1})
               un   = unique(mat,'stable');
           elseif isnumeric(mat{1}) || islogical(mat{1})
               mat  = cell2mat(mat);
               un   = unique(mat,'stable');
           end
           
           
%            mat  = string(mat);
%            mat(isnan(mat)) = -1;
%            un   = unique(mat,'stable');
           in   = zeros(1,length(un));
           for j=1:length(un)
               temp = un(j);
               if isnumeric(temp)                 
                   if temp==0
                       cond = temp == mat;
                   else
                       cond = abs(temp - mat)./max(temp,mat) < 0.01;
                   end
               elseif islogical(temp)
                   cond = temp == mat;
               else
                   temp = cell2mat(temp);
                   fun  = @(x,temp) strcmp(x,temp); 
                   temp = cellfun(@(x) fun(x,temp),mat,'UniformOutput',false);
                   cond = cell2mat(temp);
                   %cond = temp == mat;
               end
              in(j) = length(find(cond)); 
           end
%           if ischar(un)
%               un = false;
%               in = nFiles;
%            end
           if isnumeric(temp) || islogical(temp)
               nUall{i} = [un;in]; 
           else
               nUall{i} = {un;in};
           end
        end    
        nUall = cell2struct(reshape(nUall,[],1),dbnames,1);


        for i=1:length(dbnames)
            n = nUall.(dbnames{i})(2);
            if iscell(n)
               n = cell2mat(n); 
            end
           if  n~= nFiles
              nU.(dbnames{i}) =  nUall.(dbnames{i});
           end
        end

        if showdisp
            if ~isempty(nU)
                nUnames    = fieldnames(nU);
                fprintf('#Variables with multiple values = %i \n',length(nUnames))
                for i=1:length(nUnames)
                   T = table(nU.(nUnames{i}),'VariableNames',nUnames(i),'RowNames',{'Value','#files'});
                   disp(T)
                end
            else
                fprintf('All variables are unique \n') 
            end
        end
    end
    
end