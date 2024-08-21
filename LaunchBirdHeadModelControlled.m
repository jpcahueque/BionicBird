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
l_a=0.02;

alpha = deg2rad(0);

syms theta(t) gamma(t)

%% Simullation parameters

tsim=15;
dt=0.002;

%% Controlador
Kp_theta = 20; %8
Kp_gamma = 10; %2 for the real life implementation

Kv_theta = 20;
Kv_gamma = 10;

Kp = [Kp_theta 0;0 Kp_gamma];
Kv = [Kv_theta 0;0 Kv_gamma];

%% Launch Simulation
% Only one iteration
out = sim('BirdHeadModelControlled',tsim);
figure(1)
out.theta.plot
hold on
out.gamma.plot 
out.theta_dot.plot 
out.gamma_dot.plot
legend('theta','gamma','theta dot','gamma dot')

grid minor

%% Launch Simulation
% 9 iterations
Kp_theta = [1;1;1;0.5;1;2;10;20;20];
Kp_gamma = [1;1;1;0.5;1;2;10;20;20];

Kv_theta = [0.5;1;2;1;1;1;5;5;10];
Kv_gamma = [0.5;1;2;1;1;1;5;5;10];

theta = timeseries(zeros(1,9));
theta_d = theta;
theta_dot = theta;
theta_ddot = theta;
gamma = theta;
gamma_d = theta;
gamma_dot = theta;
gamma_ddot = theta;

for i=1:1:9
    Kp = [Kp_theta(i) 0;0 Kp_gamma(i)];
    Kv = [Kv_theta(i) 0;0 Kv_gamma(i)];
    out = sim('BirdHeadModelControlledResponses',tsim);
    theta(i)=out.theta;
    theta_dot(i)=out.theta_dot;
    theta_ddot(i)=out.theta_ddot;
    theta_d(i)=out.theta_d;
    gamma(i)=out.gamma;
    gamma_d(i)=out.gamma_d;
    gamma_dot(i)=out.gamma_dot;
    gamma_ddot(i)=out.gamma_ddot;
    tau_theta(i)=out.tau_theta;
    tau_gamma(i)=out.tau_gamma;
    t = out.tout;
end
%% Plot three figures with diferent responses of theta
figure
plot(t,rad2deg(theta_d(1).Data),'-')
hold on
for i=1:3
    plot(t,rad2deg(theta(i).Data),'--')
    hold on
end
grid minor
set(gca,'FontSize',18)

legend('$\theta_d$','$K_{p2}=1 K_{v2}=0.5$','$ K_{p2}=1 K_{v2}=1$','$ K_{p2}=1 K_{v2}=2$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\theta $','Interpreter','latex')

figure
plot(t,rad2deg(theta_d(1).Data),'-')
hold on
for i=4:6
    plot(t,rad2deg(theta(i).Data),'--')
    hold on
end
grid minor
set(gca,'FontSize',18)

legend('$\theta_d$','$K_{p2}=0.5 K_{v2}=1$','$ K_{p2}=1 K_{v2}=1$','$ K_{p2}=2 K_{v2}=1$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\theta $','Interpreter','latex')

figure
plot(t,rad2deg(theta_d(1).Data),'-')
hold on
for i=7:9
    plot(t,rad2deg(theta(i).Data),'--')
    hold on
end
grid minor
set(gca,'FontSize',18)

legend('$\theta_d$','$K_{p2}=10 K_{v2}=5$','$ K_{p2}=20 K_{v2}=5$','$ K_{p2}=20 K_{v2}=10$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\theta $','Interpreter','latex')

%% Plot three figures with diferent responses of gamma
figure
plot(t,rad2deg(gamma_d(1).Data),'-')
hold on
for i=1:3
    plot(t,rad2deg(gamma(i).Data),'--')
    hold on
end
ylim([5,35])
grid minor
set(gca,'FontSize',18)

legend('$\gamma_d$','$K_{p2}=1 K_{v2}=0.5$','$ K_{p2}=1 K_{v2}=1$','$ K_{p2}=1 K_{v2}=2$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\gamma $','Interpreter','latex')

figure
plot(t,rad2deg(gamma_d(1).Data),'-')
hold on
for i=4:6
    plot(t,rad2deg(gamma(i).Data),'--')
    hold on
end
ylim([5,35])
grid minor
set(gca,'FontSize',18)

legend('$\gamma_d$','$K_{p2}=0.5 K_{v2}=1$','$ K_{p2}=1 K_{v2}=1$','$ K_{p2}=2 K_{v2}=1$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\gamma $','Interpreter','latex')

figure
plot(t,rad2deg(gamma_d(1).Data),'-')
hold on
for i=7:9
    plot(t,rad2deg(gamma(i).Data),'--')
    hold on
end
ylim([10,35])
grid minor
set(gca,'FontSize',18)

legend('$\gamma_d$','$K_{p2}=10 K_{v2}=5$','$ K_{p2}=20 K_{v2}=5$','$ K_{p2}=20 K_{v2}=10$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
title('Controlled response variation of $\gamma $','Interpreter','latex')

%% Forces and velocities
figure
plot(t,rad2deg(out.theta.Data))
hold on
plot(t,rad2deg(out.gamma.Data))
plot(t,rad2deg(out.theta_d.Data))
plot(t,rad2deg(out.gamma_d.Data))
grid minor


legend('$\theta$','$\gamma$','$\theta_d$','$\gamma_d$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
ylim([-35,35])
% title('Controlled response variation of $\gamma $','Interpreter','latex')
set(gca,'FontSize',18)

figure
plot(t,out.F_bz.Data)
fontsize(18,"points")
hline = findobj(gcf, 'type', 'line');
legend('$F_{bz}$','Interpreter','latex')
figure
% plot(t,out.F_D.Data)
plot(t,out.F_mz.Data)
grid minor

% legend('$F_{bz}$','$F_m$','Interpreter','latex')
legend('$F_m$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('g','Interpreter','latex');
set(gca,'FontSize',18)
findobj(gcf)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)


figure
plot(t,rad2deg(out.theta_dot.Data))
hold on
plot(t,rad2deg(out.gamma_dot.Data))
grid minor


legend('$\dot{\theta}$','$\dot{\gamma}$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('$deg/s$','Interpreter','latex');
% title('Velocities','Interpreter','latex')
set(gca,'FontSize',18)



%% Launch Simulation
% Do an example that covers the angles that the bird can sweep
out = sim('BirdHeadModelControlledAngleRange',tsim);
figure
out.theta.plot
hold on
out.gamma.plot 
out.theta_dot.plot 
out.gamma_dot.plot
legend('theta','gamma','theta dot','gamma dot')

grid minor

%% save workspace
save("CondIni-1010Kp201Kv101.mat","out","Kp","Kv")

%% Animation
for i = 1:5:length(out.theta.Data)
    hold off
    PlotBird(out.theta.Data(i),out.gamma.Data(i),l_h,l_b);
    text(-0.05,0.15,"Timer: "+num2str(i*dt,2));
    pause(0.01) 
end

%% Plot Responses
% Required: Import the respective workspaces

load("Kp11Kv11.mat");
load("Kp11Kv21.mat");
load("Kp11Kv051.mat");
load("Kp21Kv11.mat");
load("Kp101Kv51.mat");
load("Kp201Kv51.mat");
load("Kp201Kv101.mat");

figure(1)

i = 1:75:outKp11Kv11.theta.Length;

plot(outKp11Kv11.theta_d.Time(i),outKp11Kv11.theta_d.Data(i))
hold on
% outKp11Kv11.gamma_d.plot ;


a = plot(outKp11Kv11.theta_d.Time(i),outKp11Kv11.theta.Data(i),'LineStyle','none','Marker','o','MarkerSize',4);
% b = outKp11Kv11.gamma.plot; 


a = plot(outKp11Kv11.theta_d.Time(i),outKp21Kv11.theta.Data(i),'LineStyle','none','Marker','+');
% b = outKp21Kv11.gamma.plot ;

a = plot(outKp11Kv11.theta_d.Time(i),outKp11Kv21.theta.Data(i),'LineStyle','none','Marker','<');
% b = outKp11Kv21.gamma.plot ;

a = plot(outKp11Kv11.theta_d.Time(i),outKp11Kv051.theta.Data(i),'LineStyle','none','Marker','x');
% b = outKp11Kv21.gamma.plot ;

a = plot(outKp11Kv11.theta_d.Time(i),outKp101Kv51.theta.Data(i),'LineStyle','none','Marker','diamond');
% b = outKp11Kv21.gamma.plot ;

a = plot(outKp11Kv11.theta_d.Time(i),outKp201Kv51.theta.Data(i),'LineStyle','none','Marker','>');
% b = outKp11Kv21.gamma.plot ;

a = plot(outKp11Kv11.theta_d.Time(i),outKp201Kv101.theta.Data(i),'LineStyle','none','Marker','>');
% b = outKp11Kv21.gamma.plot ;

grid minor

legend('$\theta_d$','$Kp=1 Kv=1$','$ Kp=2 Kv=1$','$ Kp=1 Kv=2$','$ Kp=1 Kv=0.5$','$Kp=10 Kv=5$','$Kp=20 Kv=5$','$Kp=20 Kv=10$','Interpreter','latex')
% legend('$\theta_d$','$\gamma_d$','$\theta Kp=1 Kv=1$','$\gamma Kp=1 Kv=1$','$\theta Kp=2 Kv=1$','$\gamma Kp=1 Kv=1$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('Radians','Interpreter','latex');
title('Controlled response variation of $\theta $','Interpreter','latex')

%% Plot the response of the ang. positions when the init. cond. are -10 and
% 10; The desired angle is 10 degrees for both
% Required: 
load("CondIni-1010Kp201Kv101.mat")
figure
out.theta.plot ;
hold on
out.gamma.plot ;
grid minor

legend('$\theta$','$\gamma$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('Radians','Interpreter','latex');
title('Controlled response of $\theta $','Interpreter','latex')

figure
out.tau_theta.plot ;
hold on
out.tau_gamma.plot ;
grid minor

legend('$\tau_\theta$','$\tau_\gamma$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('N*mm','Interpreter','latex');
title('Response of $\tau $','Interpreter','latex')
