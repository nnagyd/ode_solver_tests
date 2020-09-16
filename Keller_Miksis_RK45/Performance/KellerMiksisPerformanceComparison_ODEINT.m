function KellerMiksisPerformanceComparison_ODEINT

MPGOS_REG.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
MPGOS_REG.y = [32.1 32.0 31.9 31.7 31.8 31.8 32.1 34.0 49.2 73.8 98.6 122.5 146.6 195.4 294.1 493.2 1221.5 6557.4];

julia_cpu.x = [4        8       16      32      64      128     256 	1024        2048        4096     10000];
julia_cpu.z = [3.898	6.911	13.139	25.709	50.941	102.874	201.296	771.332     1562.726	3161.66  NaN]; 

ODEINT_CPU.x = [256 768 1536 3072 3840 5120 7680]; % 15360 30720 46080 61440 76800 92160 122880 184320 307200];
ODEINT_CPU.z = [184605 545272 1082017 2163337 2704680 3606777 5411924]/1000; 
% ODEINT_CPU.y = [4.114	7.288	13.803	26.862	53.259	105.801	222.343	853.743 1709.705	3384.723	8231.672863];

ODEINT_VCL.x = [256 768 1536 3072 3840 5120 7680 15360 30720];%46080 61440 76800 92160 122880 184320 307200];
ODEINT_VCL.z = [113602 299860 551228 1018598 1246025 1615194 2345539 4576752 8719558]/1000; 

ODEINT_THRUST.x = [2 8 16 32];
ODEINT_THRUST.z = [345 578 902 1438]; 


fig=figure(1); hold on;
p1=plot(MPGOS_REG.x,MPGOS_REG.y);

p3=plot(julia_cpu.x,julia_cpu.z);

p6=plot(ODEINT_CPU.x,ODEINT_CPU.z);
p8=plot(ODEINT_VCL.x,ODEINT_VCL.z);
p10=plot(ODEINT_THRUST.x,ODEINT_THRUST.z);

set(gca,'XScale','log','YScale','log','XLim',[1 1e7],'YLim',[1e1 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2);

set(p3,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);

set(p6, 'Color', 'yellow','LineStyle','-','LineWidth',2);
set(p8, 'Color', [1 102/255 0],'LineStyle','-','LineWidth',2);

set(p10, 'Color', 'black','Marker','o','MarkerSize', 2, 'LineStyle', '-');

xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('Keller-Miksis Runtime (s)','FontSize',20,'FontWeight','bold','FontName','Times');
hold off