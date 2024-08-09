function I = stoch(f,phi,pf)

if ~isequal(length(f),length(phi))
	error('f and phi must be of equal length')
end

N = length(f);
I = 1:N;
%rng shuffle
for i = 1:N
	yesbreak = true;
	for j = 1:N-1
		u = rand;
		if isequal(phi(I(j)),phi(I(j+1)),0) || u < pf
			if f(I(j)) > f(I(j+1))
				a = I(j); b = I(j+1);
				I(j) = b; I(j+1) = a;
				yesbreak = false;
			end
		else
			if phi(I(j)) > phi(I(j+1))
				a = I(j); b = I(j+1);
				I(j) = b; I(j+1) = a;
				yesbreak = false;
			end
		end
	end
	if yesbreak
		break
	end
	
end
