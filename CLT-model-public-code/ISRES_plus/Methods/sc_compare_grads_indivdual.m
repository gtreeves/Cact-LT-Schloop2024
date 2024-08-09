clear
clc
close all


% first load a param set that is a result of an isres run
path        = '../Results/Mats/Rng/10/';
files       = extractFileLocations(path,'mat');
filename    = files(3);

load(filename,'xb','addParams')
params = [xb,addParams];
params = 10.^params;
clear files

[F,~,solnwt] = errorCalc(log10(params));

%check_ydot_minus_f(solnwt)

%%
% gradient from finite differences
%
fprintf('Finite Differences\n')
tic
fdgrad = grad_FiniteDifferences(params,addParams,F);
T(1) = toc;
%fdgrad = fdgrad/norm(fdgrad);

fprintf('Time taken: %.1f s\n', T(1))



%% activity
% gradient from forward sensitivity method
%
fprintf('\ndsat2_gradcalc\n')
tic
fdot = dsat2_gradcalc(solnwt)';
T(2) = toc;
% fdot = fdot./norm(fdot);



differ  = norm(fdot - fdgrad);
angl    = acosd(dot(fdot,fdgrad)/(norm(fdot)*norm(fdgrad)));

fprintf('Time taken: \t\t%.1f s\n', T(2))
fprintf('norm(fdgrad - fdot): \t%.3f \n',differ)
fprintf('Angle (fdgrad, fdot): \t%.3f\n', angl)


Grads = table(fdgrad,fdot,'VariableNames',{'Finite Diff', 'fdot'});

%% gradient from adjoint method

clear AllGrads

%
% 2. inteprolating lambda
%
fprintf('\nAdjoint Method\n')
tic
[adgrad2,GradSum,Q,Lambda,lamsoln] = grad_AdjointMethod2(solnwt);
%adgrad2 = grad_AdjointMethod2(solnwt);
T(3) = toc;
%adgrad1 = adgrad1/norm(adgrad1);

check_lambdadot_minus_f(lamsoln,solnwt)




differ1  = norm(fdgrad - adgrad2);
angl1    = acosd(dot(fdgrad,adgrad2)/(norm(fdgrad)*norm(adgrad2)));
differ2  = norm(fdot - adgrad2);
angl2    = acosd(dot(fdot,adgrad2)/(norm(fdot)*norm(adgrad2)));



fprintf('Time taken: \t\t%.1f s\n', T(3))
fprintf('norm(fdgrad - adgrad2): %.3f\n',differ1)
%fprintf('norm(fdot - adgrad2): \t%.3f\n',differ2)fprintf('Angle(fdgrad, adgrad2): %.3f\n', angl1)
fprintf('Angle(fdot, adgrad2): \t%.3f\n', angl2)



adgrad_table = table(adgrad2, fdot./adgrad2,'VariableNames',{'Adjoint2 (interp lambda)'; ...
                    'fdot./adgrad'});
AllGrads = [Grads, adgrad_table];



%
% 1. interpolating concentrations
%{
fprintf('\nAdjoint1: Interpolating concentrations\n')
tic
adgrad1 = grad_AdjointMethod1(solnwt);
T(4) = toc;
%adgrad1 = adgrad1/norm(adgrad1);

differ1  = norm(fdgrad - adgrad1);
angl1    = acosd(dot(fdgrad,adgrad1)/(norm(fdgrad)*norm(adgrad1)));
differ2  = norm(fdot - adgrad1);
angl2    = acosd(dot(fdot,adgrad1)/(norm(fdot)*norm(adgrad1)));

fprintf('Time taken: \t\t%.1f s\n', T(4))
fprintf('norm(fdgrad - adgrad1): %.3f \n',differ1)
%fprintf('norm(fdot - adgrad1): \t%.3f \n',differ2)
fprintf('Angle(fdgrad, adgrad1): %.3f\n', angl1)
%fprintf('Angle(fdot, adgrad1): \t%.3f\n', angl2)

adgrad_table = table(adgrad1,'VariableNames',{'Adjoint1 (interp conc)'});
AllGrads = [AllGrads, adgrad_table];

%}










