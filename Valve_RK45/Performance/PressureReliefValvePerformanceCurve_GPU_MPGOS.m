function PressureReliefValvePerformanceCurve_GPU_MPGOS

MPGOS_REG.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
% MPGOS_REG.z = [3.2 3.1 3.2 3.2 3.2 3.2 3.3 4.1 6.6 9.2 11.9 14.6 17.4 23.2 35.1 59.0 147.2 796.8];
MPGOS_REG.y = [3.3 3.3 3.4 3.4 3.5 3.6 3.8 4.9 8.1 11.8 15.0 18.8 22.0 29.3 44.1 74.0 185.0 1000.7];

fig=figure(1); hold on;
p1=plot(MPGOS_REG.x,MPGOS_REG.y);

set(gca,'XScale','log','YScale','log','XLim',[1e0 1e7],'YLim',[1e-0 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2,'DisplayName','GPU, MPGOS, Titan Black (1707 GFLOPS)');

xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('1000\timesRK4 step (s)','FontSize',20,'FontWeight','bold','FontName','Times');