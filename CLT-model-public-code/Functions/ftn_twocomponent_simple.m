function F = ftn_twocomponent_simple(t,Y,p,Vn,Vc,Gtot,KnucG,is2xLT)
%
%
%
%
% Function used by fsolve to find the amount of Cact-LT bound to GFP vs
% free.

%
% Unpacking state variables
%
% Gc Gn Cc Cn CLGc CLGn

Cn = Y(1);
CLGn = Y(2);
Gc = Y(3);

%
% Parameters
%
kin = p(1);
kout = p(2);
kon = p(3);
KD = p(4);
Ctot = p(5);

if exist("is2xLT","var") && is2xLT == 1
    Ctot = 2*Ctot; %2xCactus
else
    is2xLT = 0;
end

%
% Mass balance relationships
%
Gn = KnucG*Gc;

fun = @(x)solveCyt(x,p,Cn,CLGn,Gc,Gtot,KnucG,Vc,Vn,is2xLT);
x0 = [KD*CLGn/(kon*kin*Cn),kout*Cn/kin];
X = fsolve(fun,x0);

CLGc = X(1);
Cc = X(2);

%
% Diffy Q's
%

dCndt = kin*Cc - kout*Cn - kon*(Cn*Gn + KD*CLGn);
dCLGndt = kin*CLGc - kout*CLGn + kon*(Cn*Gn - KD*CLGn);
dGcdt = -(kon/(1+KnucG))*((Cn - KD*CLGn) + (Cc*Gc - KD*CLGc));
F = [dCndt; dCLGndt; dGcdt];

end

