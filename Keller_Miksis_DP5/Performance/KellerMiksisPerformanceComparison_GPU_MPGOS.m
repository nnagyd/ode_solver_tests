function KellerMiksisPerformanceComparison_GPU_MPGOS

MPGOS_REG.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
% MPGOS_REG.y = [32.1 32.0 31.9 31.7 31.8 31.8 32.1 34.0 49.2 73.8 98.6 122.5 146.6 195.4 294.1 493.2 1221.5 6557.4];
MPGOS_REG.z = [34.1 34.0 33.9 33.7 33.8 33.8 34.1 36.2 52.4 78.6 105.0 130.3 156.0 208.1 313.1 524.7 1299.4 6977.5];

fig=figure(1); hold on;
% p1=plot(MPGOS_REG.x,MPGOS_REG.y);
p2=plot(MPGOS_REG.x,MPGOS_REG.z);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e1 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
% set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2,'DisplayName','GPU, MPGOS, Titan Black (1707 GFLOPS)');
set(p2,'Color',[1 0 0], 'LineStyle','-','LineWidth',2,'DisplayName','GPU, MPGOS, Titan Black (1707 GFLOPS)');

xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('runtime (s)','FontSize',20,'FontWeight','bold','FontName','Times');