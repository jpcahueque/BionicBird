function PlotBird(theta,gamma,l_h,l_b)
    body=[-0.05 0.05; 0 0; 0 0; 1 1];
    head = [0.04; 0; 0; 1];
    motor = [-0.006 0.006 0.006 -0.006 -0.006; 0 0 0 0 0; 0 0 0.01 0.01 0;1 1 1 1 1];
    propeller = [0 0 -0.023 0.023; 0 0 0 0; 0.01 0.02 0.02 0.02; 1 1 1 1];
    base = [0 0 -0.01 0.01;0 0 0 0;0 -0.15 -0.15 -0.15];

    TB=[cos(theta) 0 sin(theta) 0;0 1 0 0;-sin(theta) 0 cos(theta) 0;0 0 0 1];
    bodyR=TB*body;

    TH = [cos(theta+gamma) 0 sin(theta+gamma) bodyR(1,2);0 1 0 0;-sin(theta+gamma) 0 cos(theta+gamma) bodyR(3,2); 0 0 0 1];
    headR = TH*head;
    
    TM = [cos(theta+gamma) 0 sin(theta+gamma) headR(1);0 1 0 0;-sin(theta+gamma) 0 cos(theta+gamma) headR(3); 0 0 0 1];
    motorR = TM*motor; 

    propellerR = TM*propeller;

    plot([bodyR(1,1) bodyR(1,2)],[bodyR(3,1) bodyR(3,2)],'r-', 'LineWidth', 2);
    hold on
    plot([bodyR(1,2) headR(1)],[bodyR(3,2) headR(3)],'b-', 'LineWidth', 2);
    plot(motorR(1,:), motorR(3,:), 'b-', 'LineWidth', 1);
    plot(propellerR(1,:), propellerR(3,:), 'k-', 'LineWidth', 1);
    plot(base(1,:), base(3,:), 'k-', 'LineWidth', 1);
    hold off
    xlim([-0.1 0.15])
    ylim([-0.2 0.2])
    grid minor
    
end