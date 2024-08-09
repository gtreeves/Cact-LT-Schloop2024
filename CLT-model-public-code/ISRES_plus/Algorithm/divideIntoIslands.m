function x = divideIntoIslands(x,nIslands)

if nIslands == 1
   return 
end


lambda         = size(x,1)/nIslands;
[idx,C,sumd,d] = kmeans(x,nIslands);

for i=1:nIslands
   csize(i) = length(find(idx==i)); 
end
[~,imax] = max(csize);

flag  = true;
while flag
   i0    = idx == imax; 
   i1    = find(i0);
   d0    = d(i0,i);
   npop  = size(d0,1);
   if  npop > lambda
     [~,id]         = sort(d0);
     i1             = i1(id);
     nadd           = length(id) - lambda;
     d(:,imax)      = inf;
     dextra         = d(i1(lambda+1:lambda+nadd),:);
     [dextra,imin]  = min(dextra,[],2);
     
     idx(i1(lambda+1:lambda+nadd))  = imin;   
     for i=1:nIslands
        csize(i) = length(find(idx==i)); 
     end
     [cmax,imax] = max(csize);
     if cmax == lambda
        flag = false; 
     end
   end
end

[~,i0]  = sort(idx);
x       = x(i0,:);


%{
lb          = lu(1,:);
ub          = lu(2,:);
for i=1:nIslands
   i0    = idx == i; 
   d0    = d(i0,i);
   xtemp = x(i0,:);
   npop  = size(xtemp,1);
   if  npop > lambda
       [~,id]   = sort(d0);
       xtemp    = xtemp(id,:);
       xtemp    = xtemp(1:lambda,:);
   elseif npop < lambda
       nnew     = lambda - npop;
       dir      = rand(nnew,size(xtemp,2));
       dir      = dir./(vecnorm(dir')');   
       
       % Find the max value that can be multiplied to the recodir so the resulting
       % individuals hit the boundary.
       direct   = d0(randi(npop,nnew,1)).*dir;
       llimit   = (repmat(lb,nnew,1)- C(i,:))./(direct);
       ulimit   = (repmat(ub,nnew,1)- C(i,:))./(direct);
       betamax  = zeros(nnew,1);
       for j=1:nnew
           lim         = llimit(j,:);
           lbeta       = min(lim(lim > 0));
           lim         = ulimit(j,:);
           ubeta       = min(lim(lim > 0));
           if isempty(lbeta)
              lbeta = inf; 
           end
           if isempty(ubeta)
               ubeta = inf; 
           end
           betamax(j) = min(lbeta,ubeta);
       end 
       
       % The individuals that have betamax<1 will be unsuccessful in taking a full
       % step
       v           = betamax < 1; 
       betarand    = rand(nnew,1);
       x_          = C(i,:)  + direct;
       x_(v,:)     = C(i,:)  + betarand(v).*betamax(v).*direct(v,:);
       xtemp       = [xtemp; x_];
   end
   X = [X; xtemp];
end

[idx,C,~,d] = kmeans(X,nIslands);
%}
end