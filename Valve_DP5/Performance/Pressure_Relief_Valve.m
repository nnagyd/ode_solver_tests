function PressureReliefValveRuntimeComparison

%y: transient only   z: full problem

%mpgos
mpgos.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
%mpgos.y = [3.2 3.1 3.2 3.2 3.2 3.2 3.3 4.1 6.6 9.2 11.9 14.6 17.4 23.2 35.1 59.0 147.2 796.8];
mpgos.z = [3.3 3.3 3.4 3.4 3.5 3.6 3.8 4.9 8.1 11.8 15.0 18.8 22.0 29.3 44.1 74.0 185.0 1000.7];

%CPU_Valve_julia_ensemble.jl
julia_cpu1.x = [8        16      32      48      64      128      256     512     1024        2048        4096        10000];
%julia_cpu.y = [4.431	10.451	20.501	30.972	42.123	93.255	194.194	386.703	774.20434	1538.097223	3065.774351	NaN]; 
julia_cpu1.z = [4.826	10.853	22.482	33.579	44.409	89.889	197.314	411.064	890.36208	1884.667078	3976.813194	NaN];

%CPU_Valve_julia_loop.jl
julia_cpu2.x = [16       48      128     256     512     1024        2048        16384];
julia_cpu2.z = [1.98	 11.35	 26.06	 59.6	 117.88	 239.3		 451.8       3601];


%GPU_Valve_julia.jl   -    incorrect results
julia_gpu.x = [16	128     1024    8192];
julia_gpu.z = [313 1009 1474 1667];



fig=figure(1); hold on;
p1=plot(mpgos.x,mpgos.z);
p2=plot(julia_cpu1.x,julia_cpu1.z);
p3=plot(julia_cpu2.x,julia_cpu2.z);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e-0 1e4],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2);
set(p2,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);
set(p3,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);


xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('runtime of PRV (s)','FontSize',20,'FontWeight','bold','FontName','Times');