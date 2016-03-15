% Generalized Capture Point
% Given a CP find the suitable CoM parabolic trajectory 

CP = [3,0];


clear all
close all
clc

import casadi.*
Tst = 0.15; % stance time
Tsw = 1.0; % swing time
T = Tst + Tsw; % period corresponding to 1 step
tt = 0.005; % integration time (sampling time)

ni = 7;
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
% fn = (3*u(4)*(1-t)^2*t - 3*u(5)*(1-t)*t^2);
% ft = (3*u(6)*(1-t)^2*t - 3*u(7)*(1-t)*t^2);
fn = u(1);
ft = u(2);
M = u(3);
% add force/torque limits
u1max = 150;
u1min = 0;
u2max = 10;
u2min = -10;
u3max = 80;
u3min = -80;

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
xdot = [thetad; eq(1) - M - ft*len*cos(theta) - fn*len*sin(theta); xd; eq(2) - ft/m; zd; eq(3) + fn/m];
% Objective term
L = theta^2 + 100*u'*u;

% Continuous time dynamics
f = Function('f', {state, u}, {xdot, L});

% Control discretization
N = T/tt; % number of control intervals
IntStep = 4; % RK4 steps per interval
DT = T/N/IntStep;
X0 = MX.sym('X0', nv);
U = MX.sym('U',ni);
X = X0;
Q = 0;
tau = 0.5;
Jfinal = 0;

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
lbw = [lbw; 0; 0; 0; 0; 1; 0];
ubw = [ubw; 0; 0; 0; 0; 1; 0];
w0 = [w0; 0; 0; 0; 0; 0; 0];

% Formulate the NLP
Xk = X0;

for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], ni);
    w = {w{:}, Uk};
    lbw = [lbw; u1min; u2min; u3min; -1000; -1000; -1000; -1000];   % normal ground reaction force u(1) can only be positive
    ubw = [ubw; u1max; u2max; u3max; 1000; 1000; 1000; 1000];
    w0 = [w0; 0; 0; 0; 0; 0; 0; 0];

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
        J = J + Jfinal;
%     elseif k == 100
%         lbw = [lbw; -inf; -inf; -inf; -inf; 0;  -inf];  % z coordinate can only be positive
%         ubw = [ubw;  inf; inf; inf; inf; inf; inf];
%         w0 = [w0; 0; 0; 0; 0; 0; 0];
    else
        lbw = [lbw; -inf; -inf;  -inf;  -inf;  0;  -inf];  % z coordinate can only be positive
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
    % add complementarity constraint
%     g = {g{:}, Uk(3)*(Xk(5)-len*cos(Xk(1)))};
%     lbg = [lbg; 0];
%     ubg = [ubg; tau];
%     g = {g{:}, Uk(2)*(Xk(5)-len*cos(Xk(1)))};
%     lbg = [lbg; 0];
%     ubg = [ubg; tau];
%     g = {g{:}, Uk(1)*(Xk(5)-len*cos(Xk(1)))};
%     lbg = [lbg; 0];
%     ubg = [ubg; tau];
    % parametrize fn
    t = k*tt/Tst;
    g = {g{:}, Uk(1) - (3*Uk(4)*(1-t)^2*t - 3*Uk(5)*(1-t)*t^2)};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
    g = {g{:}, Uk(2) - (3*Uk(6)*(1-t)^2*t - 3*Uk(7)*(1-t)*t^2)};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
    if (k> 0)&&(k < Tst/tt)
        g = {g{:}, Uk(4) - tmpUk4, Uk(5) - tmpUk5};
        lbg = [lbg; 0; 0];
        ubg = [ubg; 0; 0];
        g = {g{:}, Uk(6) - tmpUk6, Uk(7) - tmpUk7};
        lbg = [lbg; 0; 0];
        ubg = [ubg; 0; 0];
    elseif k ==0
        tmpUk4 = Uk(4);
        tmpUk5 = Uk(5);
        tmpUk6 = Uk(6);
        tmpUk7 = Uk(7);
    else
    g = {g{:}, Uk(4), Uk(5)};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];
    g = {g{:}, Uk(6), Uk(7)};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];
    end
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
x1_opt = w_opt(1:13:end);
x2_opt = w_opt(2:13:end);
x3_opt = w_opt(3:13:end);
x4_opt = w_opt(4:13:end);
x5_opt = w_opt(5:13:end);
x6_opt = w_opt(6:13:end);
u1_opt = w_opt(7:13:end);
u2_opt = w_opt(8:13:end);
u3_opt = w_opt(9:13:end);
u4_opt = w_opt(10:13:end);
u5_opt = w_opt(11:13:end);
u6_opt = w_opt(12:13:end);
u7_opt = w_opt(13:13:end);

tgrid = 0:tt:T-tt;
clf;
pend = [x5_opt - len*cos(x1_opt),x5_opt];
subplot(3,1,1),hold on
plot(tgrid, x1_opt, 'k--')
plot(tgrid, x5_opt - len*cos(x1_opt), 'g')
plot(tgrid, x5_opt, 'b')
xlabel('t')
legend('theta [rad]','z foot [m]','z [m]')
xlabel('time [s]');
subplot(3,1,2),
plot(tgrid, x4_opt,'k'), hold on;
plot(tgrid, x6_opt,'r'), plot(tgrid, x2_opt,'g');
% stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x dot [m/s]','z dot [m/s]', 'theta dot [rad/s]')
xlabel('time [s]');
subplot(3,1,3);
% handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan],[u3_opt; nan]]);hold on;
stairs(tgrid, [u1_opt;nan],'r'), hold on;
stairs(tgrid, [u2_opt;nan],'k')
stairs(tgrid, [u3_opt;nan],'g')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('fn [N]', 'ft [N]', 'M [Nm]');
xlabel('time [s]');
hold off;

figure(2)
n = size(x1_opt,1);
Concatzz = []; Concatxx = [];
for k=1:n
    P0z = x5_opt(k) - len*cos(x1_opt(k));    
    P0x = x3_opt(k) - len*sin(x1_opt(k)); 
    % end effector coordinates
    zz = [x5_opt(k), P0z];
    Concatzz = [Concatzz; zz];
    xx = [x3_opt(k), P0x];
    Concatxx = [Concatxx; xx];
    % base coordinates
    floatz = [0, 0, P0z];
    floatx = [0, P0x, P0x];
    floorx = [-1 1];
    floory = [0 0];
    plot(xx,zz,'k', floatx, floatz, 'r', x3_opt(k), x5_opt(k),'ro', floorx, floory, 'k--', Concatxx', Concatzz','k')
%     axis([-2 2 -0.5 5.5])
    % Store the frame
    xlabel('x'); ylabel('z');
    title('L = \theta^2 + x^2 + (z-2)^2  (\tau = 0.5)');
    drawnow
    pause(0.1)
end
