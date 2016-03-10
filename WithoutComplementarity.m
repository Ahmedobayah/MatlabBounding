clear all
close all
clc

import casadi.*
tt = 0.1; % integration time (sampling time)
ni = 2;
% Declare model variables
theta = MX.sym('theta'); % theta1
thetad = MX.sym('thetad'); % theta1_dot
bx = MX.sym('bx');  % base x coordinate
bxd = MX.sym('bxd'); % base x dot
by = MX.sym('by'); % base y coordinate
byd = MX.sym('byd'); % base y dot

% state vector
x = [theta; thetad; bx; bxd; by; byd];
nv = size(x,1);
q = [theta; bx; by];
dq = [thetad; bxd; byd];
ddq = MX.sym('ddq',size(q,1)); % theta1_dot_dot

u = MX.sym('u',ni);  % u(1) = horizontal GRF, u(2) = normal to the ground GRF
grav = -9.81;
m = 1; % mass of the pendulum
M = 0.05; % mass of the floating base
len = 1;

% Model equations
% xdot = [x2;     - fr*x2 + grav/ll*sin(x1) - u];
E1 = 0.5*m*len*bxd^2 + 0.5*m*len*byd^2 + 0.5*m*(len^2*cos(theta)^2*thetad^2 + len^2*sin(theta)^2*thetad^2 - 2*len*cos(theta)*bxd*thetad - 2*len*sin(theta)*byd*thetad);
E2 = 0.5*M*len*bxd^2 + 0.5*M*len*byd^2;
E = E1 + E2;
V = m*grav*(cos(theta)*len + by) + M*grav*by;
Lag = E - V;

% Equation of motion
% eq = jtimes(gradient(Lag,dq),q,dq) - gradient(Lag,q);
eq = jacobian(gradient(Lag,dq),q)*dq - gradient(Lag,q);
xdot = [thetad; eq(1); bxd; eq(2)+u(1); byd; eq(3)+u(2)];
% Objective term
L = (bx-2)^2 + by^2 + theta^2;


% Continuous time dynamics
f = Function('f', {x, u}, {xdot, L});

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
lbw = [lbw; 0.1; 0; 0; 0; 0; 0];
ubw = [ubw; 0.1; 0; 0; 0; 0; 0];
w0 = [w0; 0; 0; 0; 0; 0; 0];

% Formulate the NLP
Xk = X0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], ni);
    w = {w{:}, Uk};
    lbw = [lbw; -inf; 0];   % normal ground reaction forces u(2) can only be positive
    ubw = [ubw; inf; inf];
    w0 = [w0; 0; 0];

    % Integrate till the end of the interval
    [Xk_end, Jk] = easycall(F, Xk, Uk);
    J=J+Jk;
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nv);
    w = {w{:}, Xk};
    if k == N-1
        lbw = [lbw; -inf; -inf;  -inf;  -inf;  0;  -inf];  % y base coordinate can only be positive
        ubw = [ubw; inf; inf;  inf;  inf; inf;  inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    else
        lbw = [lbw; -inf; -inf;  -inf;  -inf;  0;  -inf];  % y base coordinate can only be positive
        ubw = [ubw;  inf;  inf;  inf;  inf;  inf;  inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    end
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; 0; 0; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0; 0; 0];
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
x1_opt = w_opt(1:8:end);
x2_opt = w_opt(2:8:end);
x3_opt = w_opt(3:8:end);
x4_opt = w_opt(4:8:end);
x5_opt = w_opt(5:8:end);
x6_opt = w_opt(6:8:end);
u1_opt = w_opt(7:8:end);
u2_opt = w_opt(8:8:end);

tgrid = linspace(0, T, N+1);
clf;
subplot(2,1,1),hold on
plot(tgrid, x1_opt, 'k--')
plot(tgrid, x2_opt, 'r--')
% plot(tgrid, x3_opt, 'k^')
% plot(tgrid, x4_opt, 'r*')
plot(tgrid, x3_opt, 'g')
plot(tgrid, x5_opt, 'b')

% stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('theta','theta dot','base x','base y')
subplot(2,1,2);
handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan]]);hold on;
plot(tgrid, x3_opt,'r^')
plot(tgrid, x5_opt,'k*')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('u1', 'u2', 'base x', 'base y');
hold off;

figure(2)
n = size(x1_opt,1);
for k=1:n
    P1y = len*cos(x1_opt(k));    
    P1x = len*sin(x1_opt(k)); 
    % end effector coordinates
    yy = [x5_opt(k), P1y + x5_opt(k)];
    xx = [x3_opt(k), P1x+x3_opt(k)];
    % CoM of the pendulum
    mx = (x3_opt(k)*M + (P1x + x3_opt(k)*m)/(m+M));
    my = (x5_opt(k)*M + (P1y + x5_opt(k)*m)/(m+M));
    % base coordinates
    floatby = [0, 0, x5_opt(k)];
    floatbx = [0, x3_opt(k), x3_opt(k)];
    plot(xx,yy,'k',floatbx,floatby,'r',mx,my,'ro')
    axis([-1 x3_opt(k)+2 -1 x5_opt(k)+2])
    % Store the frame
    drawnow
    pause(0.1)
end
