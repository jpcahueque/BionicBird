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

l_a = 0.020; 
l_ej = 0.030;
alpha = 0;
syms theta(t) gamma(t)

%% Simullation parameters

tsim=15;
dt=0.002;

%% Controlador
Kp_theta = 20;
Kp_gamma = 20;

Kv_theta = 10;
Kv_gamma = 10;

Kp = [Kp_theta 0;0 Kp_gamma];
Kv = [Kv_theta 0;0 Kv_gamma];

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

%% Animation
f5 = figure('Name','PositionControlAnimation','Color','white','Position', [100 100 1000 600],'MenuBar','None');
pause(5)
for i = 1:50:length(out.theta.Data)
    subplot(3,7,[1,2,3,8,9,10,15,16,17]);
    hold off
    PlotBird(out.theta.Data(i),out.gamma.Data(i),l_h,l_b);
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$z$','Interpreter', 'latex')
    set(gca,'FontSize',16)
    % text(-0.05,0.15,"Timer: "+num2str(i*dt,2));

    subplot(3,7,[4,5,6,7]);
    plot(out.theta.Time(1:i),rad2deg(out.theta.Data(1:i)),'r-','LineWidth',2)
    hold on
    plot(out.gamma.Time(1:i),rad2deg(out.gamma.Data(1:i)),'b-','LineWidth',2)
    % grid minor
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('deg','Interpreter', 'latex')
    legend('$\theta$','$\gamma$','Location','southeast','Interpreter', 'latex')
    xlim([0,15])
    ylim([-40,40])
    set(gca,'FontSize',16)
    
    subplot(3,7,[11,12,13,14]);
    plot(out.theta.Time(1:i),rad2deg(out.theta_dot.Data(1:i)),'r-','LineWidth',2)
    hold on
    plot(out.gamma.Time(1:i),rad2deg(out.gamma_dot.Data(1:i)),'b-','LineWidth',2)
    % grid minor
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('deg','Interpreter', 'latex')
    legend('$\dot{\theta}$','$\dot{\gamma}$','Location','southeast','Interpreter', 'latex')
    xlim([0,15])
    ylim([-30,30])

    set(gca,'FontSize',16)
    
    subplot(3,7,[18,19,20,21]);
    yyaxis left
    
    plot(out.gamma.Time(1:i),out.F_mz.Data(1:i),'b-','LineWidth',2)
    % legend('$F_{mz}$','Location','southeast','Interpreter', 'latex')
    ylim([4.5,5.5])
    ylabel('g','Interpreter', 'latex')
    yyaxis right
    plot(out.theta.Time(1:i),out.F_bz.Data(1:i),'r-','LineWidth',2)
    ylim([-0.5,0.5])
    % grid minor
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    
    legend('$F_{mz}$','$F_{bz}$','Location','southeast','Interpreter', 'latex')
    xlim([0,15])
    % ylim([-0.1,5.5])

    set(gca,'FontSize',16)

    pause(0.01) 
end