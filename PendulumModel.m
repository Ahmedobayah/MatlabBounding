clear all
close all
clc

import casadi.*
tt = 0.1; % integration time (sampling time)

% Declare model variables
theta = MX.sym('theta'); % theta1
thetad = MX.sym('thetad'); % theta1_dot
x = MX.sym('x'); % x
xd = MX.sym('xd'); % x dot
z = MX.sym('z'); % z
zd = MX.sym('zd'); % z dot
X = [theta; thetad; x; xd; z; zd];
nv = size(X,1);
q = [theta; x; z];
dq = [thetad; xd; zd];
ddq = MX.sym('ddq'); % theta1_dot_dot
u = MX.sym('u');
grav = 9.81;
m = 1;
len = 0.7;
I = len^2*m;
fr = 0.5; % friction coeff
% Model equations
% xdot = [x2;     - fr*x2 + grav/ll*sin(x1) - u];
E = 0.5*I*thetad^2 + 0.5 *m * xd^2;
U = - m*grav*cos(theta)*len;
Lag = E - U;

% Equation of motion
eq = jtimes(gradient(Lag,dq),q,dq) - gradient(Lag,q);
xdot = [thetad; eq(1) - u - fr* thetad ; xd; eq(2); zd; eq(3)];
% Objective term
% L = x1^2 + x2^2;
L = cos(theta);

% Continuous time dynamics
f = Function('f', {X, u}, {xdot, L});

% Control discretization
N = 100; % number of control intervals
M = 4; % RK4 steps per interval
T = tt * N; % Time horizon
DT = T/N/M;
X0 = MX.sym('X0', nv);
U = MX.sym('U');
X = X0;
Q = 0;

for j=1:M
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
theta0 = 10/180*pi;
lbw = [lbw; theta0; 0; -sin(theta0); 0; cos(theta0); 0];
ubw = [ubw; theta0; 0; -sin(theta0); 0; cos(theta0); 0];
w0 = [w0; 0; 0; 0; 0; 0; 0];

% Formulate the NLP
Xk = X0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw; 0];
    ubw = [ubw; 0];
    w0 = [w0;  0];

    % Integrate till the end of the interval
    [Xk_end, Jk] = easycall(F, Xk, Uk);
    J=J+Jk;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nv);
    w = {w{:}, Xk};
    if k ==N-1
        lbw = [lbw; -inf; -inf; -inf; -inf; -inf; -inf];
        ubw = [ubw; inf; inf; inf; inf; inf; inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    else
        lbw = [lbw; -inf; -inf; -inf; -inf; -inf; -inf];
        ubw = [ubw;  inf;  inf; inf; inf; inf; inf];
        w0 = [w0; 0; 0; 0; 0; 0; 0];
    end
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; 0; 0; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0; 0; 0];
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
x1_opt = w_opt(1:7:end);
x2_opt = w_opt(2:7:end);
x3_opt = w_opt(3:7:end);
x4_opt = w_opt(4:7:end);
x5_opt = w_opt(5:7:end);
x6_opt = w_opt(6:7:end);
u_opt = w_opt(7:7:end);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x1_opt, 'k--')
plot(tgrid, x2_opt, 'r--')
stairs(tgrid, [u_opt;nan], '-.')
xlabel('t')
legend('theta','theta dot','u')

figure(2)
base = [0,0];
n = size(x1_opt,1);
for k=1:n
    Pz = len*cos(x1_opt(k));    
    Px = x3_opt(k); 
    xfoot = Px - len*sin(x1_opt(k));
    zfoot = Pz - len*cos(x1_opt(k));
    zz = [zfoot, Pz];
    xx = [xfoot, Px];
    plot(xx,zz)
    axis([-1 9 -5 5])
%     % Store the frame
    drawnow
    pause(0.05)
end