
figure
t = out.tout;
plot(t,rad2deg(out.theta_d.Data))
hold on
plot(out.Bird_Pitch.Time,(out.Bird_Pitch.Data))
legend('$theta_d$','$theta$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
ylim([-90,90])
set(gca,'FontSize',18)

grid minor


figure
t = out.tout;
plot(t,rad2deg(out.gamma_d.Data))
hold on
plot(out.Bird_Pitch.Time,(out.Head_Pitch.Data))
legend('$gamma_d$','$gamma$','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');
ylabel('deg','Interpreter','latex');
ylim([-90,90])
set(gca,'FontSize',18)

grid minor