function KellerMiksisRuntimeComparison

%y: transient only, z: full run

%mpgos...
mpgos.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
%mpgos.y = [32.1 32.0 31.9 31.7 31.8 31.8 32.1 34.0 49.2 73.8 98.6 122.5 146.6 195.4 294.1 493.2 1221.5 6557.4];
mpgos.z = [34.1 34.0 33.9 33.7 33.8 33.8 34.1 36.2 52.4 78.6 105.0 130.3 156.0 208.1 313.1 524.7 1299.4 6977.5];

%CPU_KellerMiksis_julia.jl
julia_cpu.x = [4	8	16	32	64	128	256	1024 2048 4096 10000];
%julia_cpu.y = [4.114	7.288	13.803	26.862	53.259	105.801	222.343	853.743 1709.705	3384.723	8231.672863];
julia_cpu.z = [4.428	7.796	14.827	29.002	57.240	113.616	226.456	914.129392 1823.728	3652.775 NaN]; 

%GPU_KellerMiksis_julia.jl
julia_gpu.x = [4	8	16	32	64	128	256	512 1024];
julia_gpu.z = [3628.48 4645.76 6697.18 10074.9 15063.4 22195.2 31552 57664 80327]; 

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