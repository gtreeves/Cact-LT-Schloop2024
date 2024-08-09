function [phi,NCRI,Y,Fval,exitflag,output] = calc_phi_NCRI(p,Y0,Gtot,KnucG,is2xLT)

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

stage = 'interphase'; m = 51;
[~,~,~,Vn,Vc] = nuclearSize(1,'static',m,stage);


% Gc Gn Cc Cn CLGc CLGn


fhandle = @(Y,p) ftn_twocomponent_simple(0,Y,p,Vn,Vc,Gtot,KnucG,is2xLT);
options = optimoptions(@fsolve,'display','off','FunctionTolerance',1e-20,'Algorithm','trust-region','OptimalityTolerance',1e-20);

[Y,Fval,exitflag,output] = fsolve(fhandle,Y0,options,p);

% Gc Gn Cc Cn CLGc CLGn

Cn = Y(1);
CLGn = Y(2);
Gc = Y(3);
Gn = KnucG*Gc;

fun = @(x)solveCyt(x,p,Cn,CLGn,Gc,Gtot,KnucG,Vc,Vn,is2xLT);
x0 = [KD*CLGn/(kon*kin*Cn),kout*Cn/kin];
F = fsolve(fun,x0);

CLGc = F(1);
Cc = F(2);

% pack these solutions in Y
Y(4) = CLGc;
Y(5) = Cc;

phi = CLGn/(Gn + CLGn);
NCRI = (CLGn + Gn)/(CLGc + Gc);
end