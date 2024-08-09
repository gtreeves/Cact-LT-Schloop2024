%clear
clc
close all


% first load a param set that is a result of an isres run
path        = '../Results/Mats/Rng/10/';
files       = extractFileLocations(path,'mat');


% 
% logging time taken to run dsat2_gradcalc
%
nfiles = length(files);
Grads2 = [];
for i=1:length(files)
   disp(['Running ',num2str(i),'/',num2str(nfiles)])
    
   filename    = files(i);
   load(filename,'xb','addParams','modelname')
   params = [xb,addParams];
   params = 10.^params;
   
   
   % first, run and get soln
   [F,~,solnwt] = errorCalc(log10(params)); 
    
   
   % finite differences
    tic
    fdgrad = grad_FiniteDifferences(params,addParams,F);
    T1(i) = toc;
    FD(:,i) = fdgrad;
   
   % adjoint method
    tic
    adgrad = grad_AdjointMethod2a(solnwt);
    T2(i) = toc;
    Adgrad(:,i) = adgrad;
   

   % dsat2_gradcalc
   tic
   fdot = dsat2_gradcalc(solnwt)';
   T3(i) = toc;
   Fdot(:,i) = fdot;
   
   grads = table(fdgrad,adgrad,fdot,fdot./adgrad,'variablenames',{'FD','Adjoint','fdot','fdot./adgrad'});
   Grads{i} = grads;
   Grads2 = [Grads2,fdgrad,adgrad,fdot,fdot./adgrad];
end

t1 = mean(T1);
t2 = mean(T2);
t3 = mean(T3);

T = table(T1',T2',T3','variablenames',{'T_FD','T_Adjoint','T_fdot'});


% 
% Comparing solutions from different techniques
%
 for i=1:length(files)
    
    
    


end














