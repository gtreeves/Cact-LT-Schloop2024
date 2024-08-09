function makeTimeCourseMovie(data)
% This function makes time course movie of the soln


% make model data ready
N11D = data.NC11;
N12D = data.NC12;
N13D = data.NC13;
N14D = data.NC14;

M = size(N14D,1);

xspan11  = linspace(0,1,size(N11D,1));
xspan12  = linspace(0,1,size(N12D,1));
xspan13  = linspace(0,1,size(N13D,1));
xspan14  = linspace(0,1,size(N14D,1));

NC11     = zeros(M,size(N11D,2));
NC12     = zeros(M,size(N11D,2));
NC13     = zeros(M,size(N11D,2));
NC14     = zeros(M,size(N11D,2));


for i=1:size(N14D,1)
    x = xspan14(i);
    for j = 1:size(N11D,2)
        NC11(i,j) = interp1(xspan11, N11D(:,j), x);
    end
    for j = 1:size(N12D,2)
        NC12(i,j) = interp1(xspan12, N12D(:,j), x);
    end 
    for j = 1:size(N13D,2)
        NC13(i,j) = interp1(xspan13, N13D(:,j), x);
    end
    for j = 1:size(N14D,2)
        NC14(i,j) = interp1(xspan14, N14D(:,j), x);
    end
end

Dlmodel = [NC11, NC12, NC13, NC14];
Dlmodel = Dlmodel'/max(max(Dlmodel));

figure
for i=1:size(Dlmodel,1)
   plot(xspan14,Dlmodel(i,:))
   ylim([0,1])
   hold on
   pause(0.005)
   %drawnow limitrate
end

end