clc;
clear all;
close all;

%% Parameters PEPE

% m_b=0.51;
% m_h=0.005;
% 
% I_by=0.00003;
% I_hy=0.00005;
% 
% g=9.81;
% 
% l_b=0.1;
% l_h=0.05;

%% Parameters

m_b=0.011;
m_h=0.005;

I_by=0.00001;
I_hy=0.00001;

g=9.81;

l_b=0.05;
l_h=0.04;


%% Controlador
Kp_theta = 10;
Kp_gamma = 10;

Kv_theta = 1;
Kv_gamma = 1;

Kp = [Kp_theta 0;0 Kp_gamma];
Kv = [Kv_theta 0;0 Kv_gamma];

%% Simullation parameters

tsim=10;
dt=0.002;

%% Launch Simulation

out = sim('BirdHeadModel3',tsim);

figure(1)
out.theta.plot
hold on
out.gamma.plot 
out.theta_dot.plot 
out.gamma_dot.plot
legend('theta','gamma','theta dot','gamma dot')

grid minor

%% Animation
for i = 1:5:length(out.theta.Data)
    hold off
    PlotBird(out.theta.Data(i),out.gamma.Data(i),l_h,l_b);
    text(-0.05,0.2,"Timer: "+num2str(i*dt,2));
    pause(0.01) 
end
