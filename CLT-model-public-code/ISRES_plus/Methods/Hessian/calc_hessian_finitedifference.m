function [Hess,Grad,f] = calc_hessian_finitedifference(x,fhandle)
% Hessian is calculated as - 
% hij = (f(x1,x2,..xi+h,xj+k..xn) - f(x1,x2,..xi+h,xj..xn) -
%        f(x1,x2,..xi,xj+k..xn) + f(x1,x2,..xi,xj..xn))/hk
% where hij is the (i,j)th element of the Hessian matrix
% 
% Note that the Hessian is in linear space. x is the list of params in
% logspace.


np      = length(x);   
pert    = 1e-3;


%
% Method 1: Vectorizing Hessian calculation
%
tperti  = pert*rand*x;
Xi      = repmat(x,np,1) + diag(tperti); 
Xij     = repmat(x,(np*np-np)/2 + np,1);
istart  = 0;
row     = [];
col     = [];
for i = 1:np
    iend            = istart+ np - i + 1;
    i0              = istart+1:iend;
    col             = [col;(i:np)'];
    row             = [row;repmat(i,np - i + 1,1)];
    Xij(i0,i)       = Xij(i0,i) + repmat(tperti(i),length(i0),1);
    Xij(i0,i:end)   = Xij(i0,i:end) + diag(tperti(i:end));
    istart          = iend;
end



%
% Get error values for all parameter sets
%
X       = [Xi;Xij;x];
F       = fhandle(X); 


%
% Calculate Hessian
%
Fi      = F(1:np);
Fij     = F(np+1:end-1);
f       = F(end);
Fr      = [];
Fc      = [];
tpertr  = [];
tpertc  = [];
for i=1:np
    rep    = np-i+1;
    Fr     = [Fr; repmat(Fi(i),rep,1)];
    tpertr = [tpertr; repmat(tperti(i),rep,1)];
    Fc     = [Fc; Fi(i:end)];
    tpertc = [tpertc; tperti(i:end)'];
end
Hessl      = (Fij - Fr - Fc + f)./(tpertr.*tpertc);
Hess       = zeros(np,np);
ind        = sub2ind([np,np],row,col);
Hess(ind)  = Hessl;
Hess       = (Hess+Hess') - eye(np).*diag(Hess);
Grad       = (Fi - f)./tperti';
Grad       = Grad';  



%
% Method 2: Calculating Hessian from finite differences
%{
Hess    = zeros(np,np);
Grad    = zeros(np,1);
istart       = 0;
tic
for i=1:np
   xi       = x;
   tperti   = pert*rand*xi(i);
   xi(i)    = xi(i) + tperti; 
   Fi       = calc_error(log10(xi));
   istart        = istart+1;

   Grad(i)  = (Fi - f)/tperti;
   for j=i:np
       xj      = x;
       tpertj  = pert*rand*xj(j);
       xj(j)   = xj(j) + tpertj; 
       Fj      = calc_error(log10(xj));
       istart       = istart+1;
       
       xij     = xi;
       xij(j)  = xij(j) + tpertj; 
       Fij     = calc_error(log10(xij));
       istart       = istart+1;

       Hess(i,j)  = (Fij - Fi - Fj + f)/(tperti*tpertj);
   end  
end
toc
Hess  = (Hess+Hess') - eye(np).*diag(Hess);
%}   
   
end  



%{

istart = 0;
np= 16;
for i=fliplr(1:np)
    i0    = istart+1:istart+i
    Terpi = repelem(tperti,i,1);
end


Xij     = repelem(Xi,np,1);
tpert   = repmat(diag(tpertj),np,1);
Xij     = Xij + tpert;

%}












