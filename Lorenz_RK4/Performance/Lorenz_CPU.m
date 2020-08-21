function LorenzRuntimeComparison

AhnertCPUv2.x = [254 507 1013 2039 4034 8121 16351 32634 64578 128895 257267 517946 1024907 2081296 4154152];
AhnertCPUv2.y = [0.00209 0.00420 0.00827 0.01678 0.03304 0.06638 0.13604 0.26260 0.53284 1.07044 2.12911 4.27722 8.50749 17.09075 33.99354];

%CPU_Lorenz_odeint_vcl.cpp
odeint_vcl_CPU.x = [256      512      1024    2048    4096    7680     15360    46080   92160   184320   307200   768000   4147200 ];
odeint_vcl_CPU.y = [.003     .006     .011    .025    .039    .057     .137     .286    .553    1.109    1.818    4.298    22.701];

%CPU_Lorenz_Julia_ensemble.jl
julia_CPU.x = [256      1024	2048	4096	7680	15360	46080	92160	184320	307200	768000	4147200];
julia_CPU.y = [0.002	0.009	0.018	0.036	0.068	0.136	0.409	0.819	1.639	2.744	6.857	37.203];

%CPU_Lorenz_Julia.jl
julia_CPU_noEnsemble.x = [256	1024	2048	4096	7680	15360	46080	92160	184320	307200	768000	4147200];
julia_CPU_noEnsemble.y = [0.002	0.009	0.018	0.037	0.068	0.138	0.417	0.826	1.654	2.754	6.891	37.220];

%CPU_Lorenz_simple_Julia.jl
julia_CPU_simple.x = [256   1024	2048	4096	7680	15360	46080	92160	184320	307200	768000	4147200];
julia_CPU_simple.y = [0.002 0.007	0.015	0.030	0.057	0.114	0.343	0.692	1.380	2.298	5.746	31.036];

%CPU_Lorenz_hand_tuned.cpp
VCL_hand_tuned_CPU.x = [256      512         1024        2048        4096        7680        15360       46080       92160       184320      307200      768000      4147200];
VCL_hand_tuned_CPU.y = [0.00067  0.001339    0.002689    0.005364    0.01077     0.020154    0.040287    0.120895    0.241721    0.483811    0.805893    2.01435     10.8739];

%CPU_Lorenz_odeint.cpp
odeint_CPU_loop.x = [256      512         1024        2048        4096        7680        15360       46080       92160       184320      307200      768000      4147200];
odeint_CPU_loop.y = [7        5           11          22           44         83           167         503         971        1928         3214        8032       43358]/1000;

fig=figure(1); hold on;
p1=plot(AhnertCPUv2.x,AhnertCPUv2.y);
%p2=plot(julia_CPU.x,julia_CPU.y);
p3=plot(odeint_vcl_CPU.x,odeint_vcl_CPU.y);
p4=plot(julia_CPU_simple.x,julia_CPU_simple.y);
%p5=plot(julia_CPU_noEnsemble.x,julia_CPU_noEnsemble.y);
p6=plot(VCL_hand_tuned_CPU.x,VCL_hand_tuned_CPU.y);
p7=plot(odeint_CPU_loop.x,odeint_CPU_loop.y);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e-3 1e2],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[0 0 0], 'LineStyle','-','LineWidth',2);
%set(p2,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);
set(p3,'Color',[0 0.6 0], 'LineStyle','-','LineWidth',2);
set(p4,'Color',[0.7 0.0 1], 'LineStyle','-','LineWidth',2);
%set(p5,'Color',[0 0.7 1], 'LineStyle','-','LineWidth',2);
set(p6,'Color',[1 0.1 0.5], 'LineStyle','-','LineWidth',2);
set(p7,'Color',[0.5 0.5 0], 'LineStyle','-','LineWidth',2);
xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('1000\timesRK4 step (s)','FontSize',20,'FontWeight','bold','FontName','Times');