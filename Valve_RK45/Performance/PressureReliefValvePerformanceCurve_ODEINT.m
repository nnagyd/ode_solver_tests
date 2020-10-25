function PressureReliefValvePerformanceCurve_ODEINT

MPGOS_REG.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
MPGOS_REG.z = [3.3 3.3 3.4 3.4 3.5 3.6 3.8 4.9 8.1 11.8 15.0 18.8 22.0 29.3 44.1 74.0 185.0 1000.7];

ODEINT_CPU.x = [256 768 1536 3072 3840 5120 7680 15360 30720]; %46080 61440]; % 76800 92160 122880 184320 307200 768000 4147200];
ODEINT_CPU.z = [9618 28925 57665 114692 143240 191204 286491 574039 1201000]/1000;

ODEINT_VCL.x = [256 768 1536 3072 3840 5120 7680 15360 30720]; %46080 61440]; % 76800 92160 122880 184320 307200 768000 4147200];
ODEINT_VCL.z = [11029 27909 52764 100145 125182 161897 236735 447463 865461]/1000;


fig=figure(1); hold on;
p1=plot(MPGOS_REG.x,MPGOS_REG.z);
p2=plot(ODEINT_CPU.x,ODEINT_CPU.z);
p4=plot(ODEINT_VCL.x,ODEINT_VCL.z);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e-0 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2);
set(p2,'Color',[0 1 0], 'LineStyle','-','LineWidth',2);
set(p4,'Color',[1 1 0], 'LineStyle','-','LineWidth',2);

xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('runtime of PRV (s)','FontSize',20,'FontWeight','bold','FontName','Times');
hold off;