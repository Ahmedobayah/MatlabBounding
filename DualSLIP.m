clear all
close all
clc

import casadi.*

m = 1;
len = 0.7;
lenTr = 0.7;
grav = - 9.81;
% Gap lenght [m]
GapLenght = 0.3;
% Obstacle height
obstacle_height = 0.3;

StabMarX = 0;
StabMarZ = 0;
% Desired step lenght
StepLenght = GapLenght + 2*StabMarX;
z_flight = obstacle_height + StabMarZ;

Tlanding = (ceil(sqrt(- z_flight/0.5/grav)*100)/100);
Tsw = 2*Tlanding;
% Horizontal speed of the pendulum
HorSpeedFinal = StepLenght/Tsw;
SpringCompr = obstacle_height/2;
K = ceil( 2*(- m*grav*(SpringCompr + z_flight)+0.5*m*HorSpeedFinal^2)/SpringCompr^2);
omega = sqrt(K/m);
freq = omega/2/pi;
Tst = 1/freq*0.5; % stance time is half of the spring period
tt = 0.005; % integration time (sampling time t = 0.005 works)
Tst = ceil(Tst*100)/100; % round the landing instant to the upper centi-second (so that the counter is integer)
T = Tst + Tsw; % period corresponding to 1 step
Duty = Tst/T;
Tliftoff = Tlanding + Tst;
% number of control inputs
ni = 4;
% Declare model variables
theta1 = MX.sym('theta1'); % theta1
thetad1 = MX.sym('thetad1'); % theta1_dot
theta2 = MX.sym('theta2'); % theta1
thetad2 = MX.sym('thetad2'); % theta1_dot
theta3 = MX.sym('theta3'); % theta1
thetad3 = MX.sym('thetad3'); % theta1_dot

x = MX.sym('x');  % base x coordinate
xd = MX.sym('xd'); % base x dot
z = MX.sym('z'); % base y coordinate
zd = MX.sym('zd'); % base y dot

% state vector
state = [theta1; thetad1; x; xd; z; zd; theta2; thetad2; theta3; thetad3];
nv = size(state,1);
q = [theta1; x; z; theta2; theta3];
dq = [thetad1; xd; zd; thetad2; thetad3];
ddq = MX.sym('ddq',size(q,1)); % theta1_dot_dot

u = MX.sym('u',ni);
% Control inputs
fnFr = u(1)*sin(theta2);
ftFr = u(1)*cos(theta2);
MFr = u(2);
fnHi = u(3)*sin(theta3);
ftHi = u(3)*cos(theta3);
MHi = u(4);
fn1 = u(1)*cos(theta1);
ft1 = u(1)*sin(theta1);
fn3 = u(1)*cos(theta3);
ft3 = u(1)*sin(theta3);

% add force/torque limits
u1max = inf;
u1min = 0;
u2max = inf;
u2min = -inf;
% Model equations
% xdot = [x2;     - fr*x2 + grav/ll*sin(x1) - u];
v_sq = xd^2 + zd^2;
E = 0.5*m*len*v_sq + 0.5*m*len^2*(thetad1^2 + thetad2^2);
V = m*grav*z;

% Lagrangian
Lag = E - V;

% Equation of motion
% eq = jtimes(gradient(Lag,dq),q,dq) - gradient(Lag,q);
eq = jacobian(gradient(Lag,dq),q)*dq - gradient(Lag,q);
% xdot = [thetad; eq(1) - M - ft*len*cos(theta) - fn*len*sin(theta); xd; eq(2) - ft/m; zd; eq(3) + fn/m];
xdot = [thetad1; eq(1) - MHi - fn1*len*sin(theta1) - ft1*len*cos(theta1); xd; eq(2) + ft1/m + ft3/m; zd; eq(3) + fn1/m + fn3/m; thetad2; eq(4) - fnHi*lenTr/2*sin(theta1 + theta2) + fnFr*lenTr/2*sin(theta1 + theta2); thetad3; eq(5) - MHi + fn3*len*sin(theta3)];
% Objective term
L = theta1^2 + u'*u;


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
% Number of steps
StepNum = 1;
% initial position
XInit = [0, StepLenght, 2*StepLenght, 3*StepLenght];
% INCREASE in the main control quantities
HorSpeedFinal =[StepLenght/T; 0; 0; 0];
HorSpeedInit = [0; StepLenght/T; 0; 0];
CaptPoint = zeros(1,StepNum);
Xfinal = [StepLenght, StepLenght, StepLenght, StepLenght];

x1_opt = []; x2_opt = []; x3_opt = []; x4_opt = [];
x5_opt = []; x6_opt = []; x7_opt = []; x8_opt = []; 
x9_opt = []; x10_opt = [];
u1_opt = []; u2_opt = []; u3_opt = []; u4_opt = [];
for ii = 1:StepNum   
    
    CaptPoint(ii) = XInit(ii) + HorSpeedFinal(ii)/sqrt(-grav/len);
    
    Jfinal = 0;
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
    lbw = [lbw; 0; 0; XInit(ii); HorSpeedInit(ii); z_flight+len; 0; 0; 0; 0; 0;];
    ubw = [ubw; 0; 0; XInit(ii); HorSpeedInit(ii); z_flight+len; 0; 0; 0; 0; 0;];
    w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
    
    % Formulate the NLP
    Xk = X0;
    
    for k=0:N-1
        % New NLP variable for the control
        Uk = MX.sym(['U_' num2str(k)], ni);
        w = {w{:}, Uk};
        lbw = [lbw; u1min; u2min; u1min; u2min];
        ubw = [ubw; u1max; u2max; u1max; u2max];
        w0 = [w0; 0; 0; 0; 0];
        
        % Integrate till the end of the interval
        [Xk_end, Jk] = easycall(F, Xk, Uk);
        J=J+Jk;
        % New NLP variable for state at end of interval
        Xk = MX.sym(['X_' num2str(k+1)], nv);
        w = {w{:}, Xk};
        if k == N-1
            lbw = [lbw; -inf; -inf;  -inf;  -inf; -inf; -inf;  -inf;  -inf; -inf; -inf];  % z coordinate can only be positive
            ubw = [ubw; inf; inf;  inf;  inf; inf; inf;  inf;  inf; inf; inf];
            w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            J = J + Jfinal;
        else
            lbw = [lbw; -inf; -inf;  -inf;  -inf;  -inf;  -inf;  -inf;  -inf;  -inf;  -inf];  % z coordinate can only be positive
            ubw = [ubw;  inf;  inf;  inf;  inf;  inf;  inf;  inf;  inf;  inf;  inf];
            w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        end
        % Add equality constraint
        g = {g{:}, Xk_end-Xk};
        lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        if k ==0
            Xinit = Xk;
            g = {g{:}, Uk(1)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
        elseif (k >= Tlanding/tt)&&(k <= Tliftoff/tt)
            if (k == Tlanding/tt)
                %     x,z coordinates of the foot at touch down
                Xtouchd = Xk(3) - len*sin(Xk(1));
                Ztouchd = Xk(5) - len*cos(Xk(1));
            end
            % current lenght of the pendulum at each iteration
            lx = (Xk(3) - Xtouchd);
            lz = (Xk(5) - Ztouchd);
            l = sqrt(lx^2+lz^2);
            % imposing the spring force to be proportional to the
            % compression of the pendulum
            g = {g{:}, Uk(1) - K*(len - l)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            % impose the contact point to be fixed during stance phase
            g = {g{:}, Xk(3) - l*sin(Xk(1)) - Xtouchd};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            g = {g{:}, Xk(5) - l*cos(Xk(1)) - Ztouchd};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            % impose the torque M to be null during stance phase
%             g = {g{:}, Uk(2)};
%             lbg = [lbg; 0];
%             ubg = [ubg; 0];
        elseif k == N-1
            g = {g{:}, Xk - Xinit};
            lbg = [lbg; 0; 0; Xfinal(ii); HorSpeedFinal(ii); 0; 0; 0; 0; 0; 0];
            ubg = [ubg; 0; 0; Xfinal(ii); HorSpeedFinal(ii); 0; 0; 0; 0; 0; 0];
            g = {g{:}, Uk(1)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        else
            g = {g{:}, Uk(1)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
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
    % allocating decision variables for plotting
    x1_opt = [x1_opt; w_opt(1:14:end)];
    x2_opt = [x2_opt; w_opt(2:14:end)];
    x3_opt = [x3_opt; w_opt(3:14:end)];
    x4_opt = [x4_opt; w_opt(4:14:end)];
    x5_opt = [x5_opt; w_opt(5:14:end)];
    x6_opt = [x6_opt; w_opt(6:14:end)];
    x7_opt = [x7_opt; w_opt(7:14:end)];
    X8_opt = [x8_opt; w_opt(8:14:end)];
    x9_opt = [x9_opt; w_opt(9:14:end)];
    x10_opt = [x10_opt; w_opt(10:14:end)];
    u1_opt = [u1_opt; w_opt(11:14:end); nan];
    u2_opt = [u2_opt; w_opt(12:14:end); nan];
    u3_opt = [u3_opt; w_opt(13:14:end); nan];
    u4_opt = [u4_opt; w_opt(14:14:end); nan];
    % Optimized Capture Point
    OptCP(ii) = x3_opt(end) + x4_opt(end)/sqrt(-grav/len);
end

%% Plot the solution
close all;

tgrid = 0:tt:StepNum*T+(StepNum-1)*tt;
clf;
pend = [[x5_opt] - len*cos(x1_opt),[x5_opt]];
subplot(3,1,1),hold on
plot(tgrid, x1_opt, 'k--')
plot(tgrid, [x5_opt] - len*cos(x1_opt), 'g')
plot(tgrid, [x5_opt], 'b')
xlabel('t')
legend('theta [rad]','z foot [m]','z [m]')
xlabel('time [s]');
subplot(3,1,2),
plot(tgrid, x4_opt,'k'), hold on;
plot(tgrid, [x6_opt],'r');
% plot(tgrid, x2_opt,'g');
% stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x dot [m/s]','z dot [m/s]')
xlabel('time [s]');
subplot(3,1,3);
Fn = [u1_opt].*cos(x1_opt);
Ft = [u1_opt].*sin(x1_opt);
% handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan],[u3_opt; nan]]);hold on;
stairs(tgrid, Fn,'r'), hold on;
stairs(tgrid, Ft,'b'),
stairs(tgrid, [u2_opt],'k')
% stairs(tgrid, [u3_opt;nan],'g')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('fn [N]', 'ft [N]', 'M [Nm]');
xlabel('time [s]');
hold off;

figure(2)
n = size(x1_opt,1);
ConcatzzStance = []; ConcatxxStance = [];
ConcatzzFly = []; ConcatxxFly = [];
for k=1:n
    P0z = x5_opt(k) - len*cos(x1_opt(k));
    P0x = x3_opt(k) - len*sin(x1_opt(k));
    % end effector coordinates
    zz = [x5_opt(k), P0z];
    xx = [x3_opt(k), P0x];
    % base coordinates
    %     floatz = [0, 0, P0z];
    %     floatx = [0, P0x, P0x];
    floorx = [-0.3 StepLenght+0.3];
    floory = [0 0];
    if (k>=Tlanding/tt)&&(k <= Tliftoff/tt) % stance phase
        if (k == Tlanding/tt)
            Xt = P0x; Zt = P0z;
        end
        Xtd = P0x-Xt;
        Ztd = P0z-Zt;
        l = sqrt(Xtd^2 + Ztd^2);
        P0z = x5_opt(k) - l*cos(x1_opt(k));
        P0x = x3_opt(k) - l*sin(x1_opt(k));
        ConcatzzStance = [ConcatzzStance; zz];
        ConcatxxStance = [ConcatxxStance; xx];
        plot(CaptPoint,0,'bo',xx,zz,'r', x3_opt(k), x5_opt(k),'ro', floorx, floory, 'k--', ConcatxxStance', ConcatzzStance','r', ConcatxxFly', ConcatzzFly','k', Xt,Zt,'ko')
    else % flight phase
        ConcatzzFly = [ConcatzzFly; zz];
        ConcatxxFly = [ConcatxxFly; xx];
        plot(OptCP,0,'bo',xx,zz,'k', x3_opt(k), x5_opt(k),'ro', floorx, floory, 'k--', ConcatxxFly', ConcatzzFly','k',ConcatxxStance', ConcatzzStance','r')
    end
    % Store the frame
    xlabel('x'); ylabel('z');
    title('L = \theta^2 + u^2; step lenght = 0.3 + 0.3 + 0.3 + 0.3 [m]');
    drawnow
    pause(0.01)
end