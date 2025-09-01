addpath("Robotics2_script/")

syms M B K D Fd Fq g delta_t delta_q delta_dt delta_dq delta_ddt delta_ddq delta_u real
eq1 = B*delta_ddt + K*(delta_t - delta_q) + D*(delta_dt - delta_dq) + Fd*delta_dt == delta_u;
eq2 = M*delta_ddq + K*(delta_q - delta_t) + D*(delta_dq - delta_dt) + Fq*delta_dq -M*g == 0;

states = [delta_t; delta_q; delta_dt; delta_dq];
derivs = [delta_ddt, delta_ddq];

f = equations_to_state([eq1, eq2], states, derivs, delta_u);
disp(f)

[A, B] = linearize_system(f, states, delta_u, [0;0;0;0], 0);
disp(A); 
disp(B);

tf_and_rootlocus(A, B, [0;1;0;0], 0);

[Co, detCo, is_ctrl] = check_controllability(A, B, 4); 

% 
% syms x(t) q(t) u(t) a b c d k real
% eqs = { diff(x,2) + a*diff(x,1) + b*x + k*(x - q) == u, ...
%         diff(q,1) == -c*q + d*x };

% sys = sym_ode_to_state(eqs, {u(t)}, {x, q});   % y = [x; q]
% 
% % Jacobiane simboliche:
% sys.A_sym, sys.B_sym, sys.C_sym, sys.D_sym

% % Linearizzazione all'origine con parametri fissati:
% eqp = struct('x0', sym([0;0;0]), 'u0', sym(0), 'a',1, 'b',2, 'c',0.5, 'd',1, 'k',3);
% sys2 = sym_ode_to_state(eqs, {u(t)}, {x, q}, eqp);
% sys2.A, sys2.B

%% metodo per Lyapunov: 
syms theta q M g theta_d Kp K real 
V = 0.5*K*(theta-q)^2 + 0.5*Kp*(theta_d - theta)^2 - M*g*(theta_d+q-theta); %la Lyapunov candidate calcolata in theta_dot e q_dot nulli
g = gradient(V, [theta,q]);
H = hessian(V, [theta, q]); 

sol = solve(g==0, [theta, q], 'Real', true); 
theta_star = sol.theta;
q_star = sol.q; 
disp(theta_star); 
disp(q_star); 
disp(det(H));


%% dinamica 2R con Payload aggiunto: 

syms q1 q2 l1 dc1 dc2 dq1 dq2 Ic1 Ic2 m1 m2 ddq1 ddq2 real
A1 = DHmatrix(0, l1, 0, q1, true);
A2 = DHmatrix(0, dc2, 0, q2, true);
disp(simplify(A1*A2))
T = simplify(A1*A2);
pc2_dot = [-l1*sin(q1)*dq1 - sin(q1+q2)*dc2*(dq1+dq2); 
    l1*cos(q1)*dq1 + cos(q1+q2)*dc2*(dq1+dq2)]; 
T2 = (pc2_dot.'*pc2_dot); 
T2 = simplify(T2);
disp(T2);

%pc3, usa stesso sistema di riferimento del secondo COM ma in posizione
%"avanzata", quindi: 

A1 = DHmatrix(0, l1, 0, q1, true);
A2 = DHmatrix(0, l2+dc3, 0, q2, true);
disp(simplify(A1*A2))
T = simplify(A1*A2);
pc3_dot = [-sin(q1+q2)*(dc3+l2)*(dq1+dq2)-l1*cos(q1)*dq1; 
    cos(q1+q2)*(dc3+l2)*(dq1+dq2)+l1*cos(q1)*dq1]; 

T3 = collect(expand(pc3_dot.'*pc3_dot), [dq1^2, dq2^2, dq1*dq2]); 
T1 = 0.5*(Ic1+m1*dc1^2)*dq1^2; 
T2 = 0.5*(dq1^2*(Ic2+m2*dc2^2+l1^2+2*cos(q2)*dc2*l1) + dq2^2*(Ic2+m2*dc2^2)+dq1*dq2*(2*dc2^2+2*cos(q2)*dc2*l1+Ic2));
T_kin = T1+T2; 
M = inertia_matrix_from_kinetic_energy(T_kin, [dq1;dq2]);
disp(M)
[c, C] = inertia_matrix_to_coriolis(M, [dq1;dq2], [ddq1;ddq2]);


syms w t real
a = int(sin(w*t),0,t); 
disp(a)
