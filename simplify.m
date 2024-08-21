syms l_h l_b m_h m_b x g
syms I_bx I_by I_bz I_hx I_hy I_hz
syms theta(t) gamma(t) theta_dot(t) gamma_dot(t) theta_ddot(t) gamma_ddot(t)
syms tau_theta tau_gamma
R_eta = [cos(theta), 0, sin(theta);
                  0, 1, 0;    
        -sin(theta), 0, cos(theta)];
ex = [1; 0; 0];
ey = [0; 1; 0];
ez = [0; 0; 1];
R_gamma_y = [cos(gamma), 0, sin(gamma);
                      0, 1, 0;    
            -sin(gamma), 0, cos(gamma)];
W_n = [1, 0, -sin(theta);
       0, 1, 0;    
       0, 0, cos(theta)];
Ib = diag([I_bx,I_by,I_bz]);
Ih = diag([I_hx,I_hy,I_hz]);

q = [0; 0; 0; 0; theta; 0; gamma];
q_dot = [0; 0; 0; 0; theta_dot; 0; gamma_dot];
q_ddot = [0; 0; 0; 0; theta_ddot; 0; gamma_ddot];

%% M_xi
M_xi = diag((m_b+m_h)*ones(3,1));

%% M_eta    
x = ex*l_b+R_gamma_y*ex*l_h;
xbody = formula(x);
S_x = [0, -xbody(3), xbody(2); 
       xbody(3), 0, -xbody(1); 
       -xbody(2), xbody(1), 0 ];
M_eta = W_n.'*(Ib - m_h*S_x*S_x+R_gamma_y*Ih*R_gamma_y.')*W_n;

%% M_gamma
m_gamma = m_h*l_h^2 + I_hy;

%% M_xi_eta
x = (ex*l_b + R_gamma_y*ex*l_h);
xbody = formula(x);
S_x = [0, -xbody(3), xbody(2); 
       xbody(3), 0, -xbody(1); 
       -xbody(2), xbody(1), 0 ];
M_xi_eta = -m_h*R_eta*S_x*W_n;

%% M_xi_gamma
M_xi_gamma = -m_h*R_eta*R_gamma_y*ez*l_h;

%% M_eta_gamma
M_eta_gamma = W_n.'*(ey*I_hy-m_h*S_x*R_gamma_y*ez*l_h);

%% M_q
M_q = [M_xi M_xi_eta M_xi_gamma; M_xi_eta.' M_eta M_eta_gamma; M_xi_gamma.' M_eta_gamma.' m_gamma];

%% Aceleraciones
A = M_q*q_ddot;

%% Coriolis
y = q_dot.'*M_q*q_dot;
J = jacobian(y,[gamma]);
Jbody = formula(J);
C = diff(M_q,t)*q_dot - 0.5*[0;0;0;0;0;0;Jbody(1)];
Cbody = formula(C)
%% G_q
G_q = [M_xi M_xi_eta M_xi_gamma].'*[0; 0; g];
Gbody = formula(G_q);
%% 

Tau = A + C + G_q;
Taubody = formula(Tau);
eqn = Tau == [0;0;0;0;tau_theta;0;tau_gamma];

simplify(Taubody(5))
% simplify(Taubody(7));

%% 
% eqn1 = Taubody(7) == tau_gamma;
% solve(eqn1,gamma,'PrincipalValue',true);
% [y1] = solve(eqn,gamma,Real=true)
% simplify(Taubody(7));
theta_dot = diff(theta,t);
theta_ddot = diff(theta_dot,t);
gamma_dot = diff(gamma,t);
gamma_ddot = diff(gamma_dot,t);
I_hy = 0.0001;
I_by = 0.0001;
l_h = 0.035;
l_b = 0.030;
m_h = 0.0049;
g = -9.8;
Kp = 100;
Kv = 10;

% function theta_ddot = odefun(t,y)
% dydt = 5*y-3;
% end
% ode1 = theta_ddot == (tau_gamma - ((I_hy + l_h^2*m_h)*gamma_ddot + (g*l_h*m_h*cos(gamma + theta)) + l_b*l_h*m_h*sin(gamma_dot)*theta_dot^2 ))/ (I_hy+l_h^2*m_h*sin(2*gamma) + l_b*l_h*m_h*sin(gamma));
% ode2 = gamma_ddot == (tau_theta - ((I_hy+I_by+(l_b^2+l_h^2*m_h+ 2*l_b*l_h*m_h*cos(gamma))*theta_ddot + (g*(l_h*m_h*sin(gamma - theta)-l_b*m_h*sin(theta))) + (2*l_h^2*m_h*cos(2*gamma) + l_b*l_h*m_h*cos(gamma))*gamma^2 - 2*l_b*l_h *m_h*sin(gamma)*theta_dot*gamma_dot)))/(I_hy+l_h^2*m_h*sin(2*gamma)+l_b*l_h*m_h*sin(gamma));
% 
% odes = [ode1; ode2];
% 
% S = dsolve(odes)

%%

m11 = I_hy+I_by+(l_b^2+l_h^2)*m_h+2*l_b*l_h*m_h*cos(gamma);
m12 = I_hy+l_h^2*m_h*sin(2*gamma)+l_b*l_h*m_h*sin(gamma);
m22 = I_hy+l_h^2*m_h;

M=[m11 m12;
   m12 m22];


c11 = -2*l_b*l_h*m_h*sin(gamma)*diff(gamma,t);
c12 = (2*l_h^2*m_h*cos(2*gamma)+l_b*l_h*m_h*cos(gamma))*diff(gamma,t);
c21 = l_b*l_h*m_h*sin(gamma)*diff(theta,t);
c22 = 0;

C=[c11 c12; 
   c21 c22];
M_dot = diff(M,t)
C+C.'