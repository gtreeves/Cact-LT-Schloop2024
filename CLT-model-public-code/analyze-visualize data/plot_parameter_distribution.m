clc
clear
close all

load("two_component_model_simple_normalized_GFPsd_rng2.mat")

% figure
for i=1:length(datnew)
    BM(i) = datnew(i).BestMin;
    XB(i,:) = datnew(i).xb;
    G(i) = datnew(i).Gtot;
%     hold on
%     plot(dat(i).Stats.Min)
end

figure
lb = [0.016  0.1676 1e-3  1e-2   1];
ub = [0.0680  0.3101  1e2   1e2   1e3];

plot_sets = XB;
boxchart(plot_sets)

x = 1:size(plot_sets,2);
% X = {'1','2','3','4','5','6','7','8'};
% X = categorical(X);
% T = [x' plot_sets'];
% T = array2table(T);
x1 = repmat(x,size(plot_sets,1),1);
% BMrep = repmat(BM,1,length(BM));
for i = 1:length(x)
    hold on
    swarmchart(i*ones(size(plot_sets,1),1),plot_sets(:,i),[],BM,'filled')
end
% p = [kin kout kon koff Ctot]
ax = gca;
ax.YScale = 'log';
xticklabels({'k_{in} [min^{-1}]','k_{out} [min^{-1}]','k_{on} [nM^{-1} min^{-1}]','K_{D}','Cact_{tot} [nM]'})
hold on
xint = [1 x(2:end-1) x(end) x(end) x(end-1:-1:2) 1];
% lb = [1e-2*ones(1,14) 1e-2*ones(1,2) 1e2*ones(1,4)];
% ub = [1e2*ones(1,14) 1e1*ones(1,2) 1e4*ones(1,4)];
yint = [lb(1:end) ub(end:-1:1)];
p1 = fill(xint,yint,[0.6 0.6 0.6]);
p1.FaceColor = [0.6 0.6 0.6];      
p1.EdgeColor = 'none'; 
p1.FaceAlpha = 0.2;
cb = colorbar;
cb.Label.String = 'Error';
set(gca,'FontSize',14,'FontName','Arial')