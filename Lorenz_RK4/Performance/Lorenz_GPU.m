function LorenzRuntimeComparison

%MPGOS
MPGOS_REG.x = [256 768 1536 3072 3840 5120 7680 15360 30720 46080 61440 76800 92160 122880 184320 307200 768000 4147200];
MPGOS_REG.y = [0.000479 0.000474 0.000475 0.000490 0.000494 0.000546 0.000648 0.001130 0.002170 0.003225 0.004278 0.005342 0.006382 0.008511 0.012742 0.021189 0.054579 0.305186];

%HAND-Tuned GPU
Spec_V1.x = [256 768 1536 3072 3840 5120 7680 15360 46080 92160 184320 307200 768000 4147200];
Spec_V1.y = [0.000239 0.000240 0.000242 0.000262 0.000263 0.000332 0.000427 0.000846 0.002538 0.005067 0.009288 0.016163 0.043672 0.246159];

AhnertVexCLv2.x = [255 503 1004 2021 4066 8041 16180 32560 64987 131954 263353 530177 1031136 2076033 4179501];
AhnertVexCLv2.y = [0.0131 0.0128 0.0128 0.0126 0.0128 0.0127 0.0132 0.0156 0.0307 0.0526 0.0957 0.1828 0.3564 0.7527 1.4675];

AhnertThrust.x = [255 504 1014 2023 4071 8191 16204 32619 64541 131057 261578 522086 1042054 2079883 4151397];
AhnertThrust.y = [0.0484 0.0504 0.0545 0.0545 0.0555 0.0607 0.0834 0.1458 0.2650 0.4916 0.9584 1.8687 3.7168 7.3928 15.0007];

%GPU_Lorenz_julia.jl
julia_GPU.x = [256	1024	2048	4096	7680	15360	46080	92160	184320	307200	768000	4147200];
julia_GPU.y = [0.160	0.188	0.239	0.281	0.220	0.257	0.483	0.859	1.537	2.377	6.090	31.611];

%GPU_Lorenz_thrust_odeint.cu
odeint_GPU_Thrust.x = [256      512         1024        2048        4096        7680        15360       46080       92160       184320      307200      768000      4147200];
odeint_GPU_Thrust.y = [97        39           38          38         42          43            59          136         250        470          760        1863        9920]/1000;

fig=figure(1); hold on;
p1=plot(MPGOS_REG.x,MPGOS_REG.y);
p2=plot(AhnertVexCLv2.x,AhnertVexCLv2.y);
p3=plot(AhnertThrust.x,AhnertThrust.y);
p4=plot(Spec_V1.x,Spec_V1.y);
p5=plot(julia_GPU.x,julia_GPU.y);
p6=plot(odeint_GPU_Thrust.x,odeint_GPU_Thrust.y);

set(gca,'XScale','log','YScale','log','XLim',[1e2 1e7],'YLim',[1e-4 1e3],'XGrid','on','YGrid','on','Box','on','FontSize',14,'FontWeight','bold','FontName','Times');
set(p1,'Color',[1 0 0], 'LineStyle','-','LineWidth',2);
set(p2,'Color',[0 0 0], 'LineStyle','-','LineWidth',2);
set(p3,'Color',[0.6 0.5 0], 'LineStyle','-','LineWidth',2);
set(p4,'Color',[1 0.3 0.3], 'LineStyle','-','LineWidth',2);
set(p5,'Color',[0 0 1], 'LineStyle','-','LineWidth',2);
set(p6,'Color',[0 0.7 0.1], 'LineStyle','-','LineWidth',2);
xlabel('N','FontSize',20,'FontWeight','bold','FontName','Times');
ylabel('1000\timesRK4 step (s)','FontSize',20,'FontWeight','bold','FontName','Times');