function [X,F,Phi] = getDatafromMatfile(file)

load(file);
struct2vars(Stats)
X    = permute(X,[1,3,2]);
X    = reshape(X,[],nParams); 
F    = reshape(F,[],1);
Phi  = permute(Phi,[1,3,2]);
Phi  = reshape(Phi,[],2); 
    
end