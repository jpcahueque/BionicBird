clc;
clear all;
close all;

%% Parameters

m_b=0.011;
m_h=0.005;

I_by=0.00001;
I_hy=0.00001;

g=9.81;

l_b=0.05;
l_h=0.04;

l_ec = 0.020; 
l_ej = 0.030;

syms theta(t) gamma(t)

%% Simullation parameters

tsim=15;
dt=0.002;

%% Launch Simulation
% 9 iterations
Kp_theta = [10;20;20];
Kp_gamma = [10;20;20];

Kv_theta = [5;5;10];
Kv_gamma = [5;5;10];
iK = 3;

%initial conditions in degrees

init_cond = deg2rad([0 0;-30 30;-15 15]);
iInitCond = 1;

%frequencies in Hz
frequencies = [0.05 0.15 0.25]; %hz
iFreq = 3;

amplitude_theta = deg2rad(10);
amplitude_gamma = deg2rad(-10);

theta = timeseries(zeros(9,1));
theta_d = theta;
theta_dot = theta;
theta_ddot = theta;
gamma = theta;
gamma_d = theta;
gamma_dot = theta;
gamma_ddot = theta;

for k = 1:iFreq
    for j=1:iInitCond
        for i=1:1:iK
            f = frequencies(k);
            theta_init = init_cond(j,1);
            gamma_init = init_cond(j,2);
            Kp = [Kp_theta(i) 0;0 Kp_gamma(i)];
            Kv = [Kv_theta(i) 0;0 Kv_gamma(i)];
            out = sim('BirdHeadModelControlled_MovementControl',tsim);
            theta((k-1)*iFreq+i)=out.theta;
            theta_dot((k-1)*iFreq+i)=out.theta_dot;
            theta_ddot((k-1)*iFreq+i)=out.theta_ddot;
            theta_d((k-1)*iFreq+i)=out.theta_d;
            gamma((k-1)*iFreq+i)=out.gamma;
            gamma_d((k-1)*iFreq+i)=out.gamma_d;
            gamma_dot((k-1)*iFreq+i)=out.gamma_dot;
            gamma_ddot((k-1)*iFreq+i)=out.gamma_ddot;
            tau_theta((k-1)*iFreq+i)=out.tau_theta;
            tau_gamma((k-1)*iFreq+i)=out.tau_gamma;
            t = out.tout;
        end
    end
end
%% Plot three figures with frequency 0.25
clc
close all
figure
plot(t,theta_d(1).Data,'-')
hold on
for i=1:3
    plot(t,theta(i).Data,'--')
    hold on  
end
grid minor

legend('$\theta_d$','$K_{p2}=10 K_{v2}=5$','$ K_{p2}=20 K_{v2}=5$','$ K_{p2}=20 K_{v2}=10$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('Radians','Interpreter','latex');
title('Controlled $\theta $, Freq = 0.05Hz','Interpreter','latex')

figure
plot(t,theta_d(4).Data,'-')
hold on
for i=4:6
    plot(t,theta(i).Data,'--')
    hold on
end
grid minor

legend('$\theta_d$','$K_{p2}=10 K_{v2}=5$','$ K_{p2}=20 K_{v2}=5$','$ K_{p2}=20 K_{p2}=10$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('Radians','Interpreter','latex');
title('Controlled $\theta $ Freq = 0.15Hz','Interpreter','latex')

figure
plot(t,rad2deg(theta_d(7).Data),'-')
hold on
for i=7:9
    plot(t,rad2deg(theta(i).Data),'--')
    hold on
end
grid minor

legend('$\theta_d$','$K_{p2}=10 K_{v2}=5$','$ K_{p2}=20 K_{v2}=5$','$ K_{p2}=20 K_{v2}=10$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled $\theta $ Freq = 0.25Hz','Interpreter','latex')
ylim([-25,25])
set(gca,'FontSize',18)

%% Plot Position Velocities and Force
figure
plot(t,rad2deg(theta(7).Data))
hold on
plot(t,rad2deg(gamma(7).Data))
plot(t,rad2deg(out.theta_d.Data))
plot(t,rad2deg(out.gamma_d.Data))
legend('$\theta$','$\gamma$','$\theta_d$','$\gamma_d$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Freq = 0.25Hz','Interpreter','latex')
set(gca,'FontSize',18)
grid minor

figure
plot(t,rad2deg(theta_dot(7).Data))
hold on
plot(t,rad2deg(gamma_dot(7).Data))
legend('$\dot{\theta}$','$\dot{\gamma}$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('$deg/s$','Interpreter','latex');
set(gca,'FontSize',18)

grid minor

% figure
% plot(t,out.F_bz.Data)
% fontsize(18,"points")
% hline = findobj(gcf, 'type', 'line');
% legend('$F_{bz}$','Interpreter','latex')
% grid minor
% 
% figure
% plot(t,out.F_mz.Data,'r')
% legend('$F_m$','Interpreter','latex')
% xlabel('Time (s)','Interpreter','latex');
% ylabel('g','Interpreter','latex');
% set(gca,'FontSize',18)
% findobj(gcf)
% hline = findobj(gcf, 'type', 'line');
% set(hline,'LineWidth',2)
% grid minor

figure
plot(t,out.F_bz.Data)
ylabel('g','Interpreter','latex');
yyaxis right
plot(t,out.F_mz.Data,'r')
legend('$F_{bz}$','$F_m$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('g','Interpreter','latex');
set(gca,'FontSize',18)
findobj(gcf)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
grid minor
