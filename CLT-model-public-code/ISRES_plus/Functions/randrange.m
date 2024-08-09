function N = randrange(r0,r1,n)

if r0 > r1
   t  = r1;
   r1 = r0;
   r0 = t;
end

if ~exist('n','var')
   n = 100; 
end

N = r0 + (r1-r0).*rand(n,1);

end

