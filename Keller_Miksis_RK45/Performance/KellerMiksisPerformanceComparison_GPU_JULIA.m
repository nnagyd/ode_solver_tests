function KellerMiksisRuntimeComparison

%y: transient only, z: full run

%mpgos...
mpgos.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
%mpgos.y = [32.1 32.0 31.9 31.7 31.8 31.8 32.1 34.0 49.2 73.8 98.6 122.5 146.6 195.4 294.1 493.2 1221.5 6557.4];
mpgos.z = [34.1 34.0 33.9 33.7 33.8 33.8 34.1 36.2 52.4 78.6 105.0 130.3 156.0 208.1 313.1 524.7 1299.4 6977.5];

%CPU_KellerMiksis_julia.jl
julia_cpu.x = [4        8       16      32      64      128     256 	1024        2048        4096     10000];
julia_cpu.z = [3.898	6.911	13.139	25.709	50.941	102.874	201.296	771.332     1562.726	3161.66  NaN]; 

%GPU_KellerMiksis_julia.jl
julia_gpu.x = [4	8	16	32	64	128	256	512 1024];
julia_gpu.z = [3149.76, 3960.32, 5799.04, 8655.04, 12756.8, 19045.4, 27221.8, 46577.3, 69702.7]; 

fig=figure(1); hold on;
p1=plot(mpgos.x,mpgos.z);
p2=plot(julia_cpu.x,julia_cpu.z);
p3=plot(julia_gpu.x,julia_gpu.z);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e1 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2);
set(p2,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);
set(p3,'Color',[0.5 0 0.5], 'LineStyle','-','LineWidth',2);

xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('runtime (s)','FontSize',20,'FontWeight','bold','FontName','Times');