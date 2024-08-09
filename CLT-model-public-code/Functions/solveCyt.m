function F = solveCyt(x,p,Cn,CLGn,Gc,Gtot,KnucG,Vc,Vn,is2xLT)

Vtot = Vc + Vn;
Ctot = p(5);

if exist("is2xLT","var") && is2xLT == 1
    Ctot = 2*Ctot; %2xCactus
end

Gn = KnucG*Gc;

CLGc = x(1);
Cc = x(2);

F(1) = (Vc/Vtot)*(CLGc + Cc) + (Vn/Vtot)*(CLGn + Cn) - Ctot;
F(2) = (Vc/Vtot)*(CLGc + Gc) + (Vn/Vtot)*(CLGn + Gn) - Gtot;

end