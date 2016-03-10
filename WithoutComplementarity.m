clear all
close all
clc

import casadi.*
tt = 0.05; % integration time (sampling time)
ni = 3;
% Declare model variables
theta = MX.sym('theta'); % theta1
thetad = MX.sym('thetad'); % theta1_dot
x = MX.sym('x');  % base x coordinate
xd = MX.sym('xd'); % base x dot
z = MX.sym('z'); % base y coordinate
zd = MX.sym('zd'); % base y dot

% state vector
state = [theta; thetad; x; xd; z; zd];
nv = size(state,1);
q = [theta; x; z];
dq = [thetad; xd; zd];
ddq = MX.sym('ddq',size(q,1)); % theta1_dot_dot

u = MX.sym('u',ni);  % u(1) = f (GRF), u(2) = M (Momentum of the flywheel)
fn = u(1);
ft = u(2);
M = u(3);

grav = -9.81;
m = 1; % mass of the pendulum
len = 1;

% Model equations
% xdot = [x2;     - fr*x2 + grav/ll*sin(x1) - u];
v_sq = xd^2 + zd^2;
E = 0.5*m*len*v_sq + 0.5*m*len^2*thetad^2;
V = m*grav*z;

% Lagrangian
Lag = E - V;


% Equation of motion
% eq = jtimes(gradient(Lag,dq),q,dq) - gradient(Lag,q);
eq = jacobian(gradient(Lag,dq),q)*dq - gradient(Lag,q);
xdot = [thetad; eq(1) - M - ft*len*cos(theta) - fn*len*sin(theta); xd; eq(2) - ft; zd; eq(3) + fn];
% Objective term
L = theta^2 + thetad^2 + (x+0.5)^2 + xd^2 + (z-1)^2 + zd^2;


% Continuous time dynamics
f = Function('f', {state, u}, {xdot, L});

% Control discretization
N = 100; % number of control intervals
IntStep = 4; % RK4 steps per interval
T = tt * N; % Time horizon
DT = T/N/IntStep;
X0 = MX.sym('X0', nv);
U = MX.sym('U',ni);
X = X0;
Q = 0;
tau = 0.1;
% Runge Kutta 4 integrator
for j=1:IntStep
    [k1, k1_q] = easycall(f, X, U);
    [k2, k2_q] = easycall(f, X + DT/2 * k1, U);
    [k3, k3_q] = easycall(f, X + DT/2 * k2, U);
    [k4, k4_q] = easycall(f, X + DT * k3, U);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = Function('F', {X0, U}, {X, Q});

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0', nv);
w = {w{:}, X0};
lbw = [lbw; 0.1; 0; 1; 0; 2; 0];
ubw = [ubw; 0.1; 0; 1; 0; 2; 0];
w0 = [w0; 0; 0; 0; 0; 0; 0];

% Formulate the NLP
Xk = X0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], ni);
    w = {w{:}, Uk};
    lbw = [lbw; 0; -inf; -inf];   % normal ground reaction force u(1) can only be positive
    ubw = [ubw; inf; inf; inf];
    w0 = [w0; 0; 0; 0];

    % Integrate till the end of the interval
    [Xk_end, Jk] = easycall(F, Xk, Uk);
    J=J+Jk;
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nv);
    w = {w{:}, Xk};
    if k == N-1
        lbw = [lbw; -inf; -inf;  -inf;  -inf;  -inf;  -inf];  % z coordinate can only be positive
        ubw = [ubw; inf; inf;  inf;  inf; inf;  inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    else
        lbw = [lbw; -inf; -inf;  -inf;  -inf;  -inf;  -inf];  % z coordinate can only be positive
        ubw = [ubw;  inf;  inf;  inf;  inf;  inf;  inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    end
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; 0; 0; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0; 0; 0];
    % the contact point must always be positive
    g = {g{:}, Xk(5)-len*cos(Xk(1))};
    lbg = [lbg; 0];
    ubg = [ubg; inf];
%     % add complementarity constraint
%     g = {g{:}, Uk(2)*Xk(5)};
%     lbg = [lbg; 0];
%     ubg = [ubg; tau];
%     g = {g{:}, Uk(1)*Xk(5)};
%     lbg = [lbg; 0];
%     ubg = [ubg; tau];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
arg = struct('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
sol = solver(arg);
w_opt = full(sol.x);

% Plot the solution
close all;
x1_opt = w_opt(1:9:end);
x2_opt = w_opt(2:9:end);
x3_opt = w_opt(3:9:end);
x4_opt = w_opt(4:9:end);
x5_opt = w_opt(5:9:end);
x6_opt = w_opt(6:9:end);
u1_opt = w_opt(7:9:end);
u2_opt = w_opt(8:9:end);
u3_opt = w_opt(9:9:end);

tgrid = linspace(0, T, N+1);
clf;
subplot(2,1,1),hold on
plot(tgrid, x1_opt, 'k--')
% plot(tgrid, x2_opt, 'r--')
% plot(tgrid, x3_opt, 'k^')
% plot(tgrid, x4_opt, 'r*')
plot(tgrid, x3_opt, 'g')
plot(tgrid, x5_opt, 'b')

% stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('theta','x','z')
subplot(2,1,2);
% handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan],[u3_opt; nan]]);hold on;
stairs(tgrid, [u1_opt;nan],'r'), hold on;
stairs(tgrid, [u2_opt;nan],'k')
stairs(tgrid, [u3_opt;nan],'g')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('fn', 'ft', 'M');
hold off;

figure(2)
n = size(x1_opt,1);
for k=1:n
    P0z = x5_opt(k) - len*cos(x1_opt(k));    
    P0x = x3_opt(k) - len*sin(x1_opt(k)); 
    % end effector coordinates
    zz = [x5_opt(k), P0z];
    xx = [x3_opt(k), P0x];
    % base coordinates
    floatz = [0, 0, P0z];
    floatx = [0, P0x, P0x];
    floorx = [-1 1];
    floory = [0 0];
    plot(xx,zz,'k', floatx, floatz, 'r', x3_opt(k), x5_opt(k),'ro', floorx, floory, 'k--')
%     axis([-1 1 -0.5 1.5])
    % Store the frame
    drawnow
    pause(0.1)
end
