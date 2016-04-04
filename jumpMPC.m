clear all
close all
clc

import casadi.*

m = 1;
len = 0.7;
grav = - 9.81;
I = m*len^2;
% Gap lenght [m]
GapLenght = 1;
% Obstacle height
obstacle_height = 0.5;

StabMarX = 0;
StabMarZ = 0;
% Desired step lenght
StepLenght = GapLenght + 2*StabMarX;
z_flight = obstacle_height + StabMarZ;

Tlanding1 = (ceil(sqrt(- z_flight/0.5/grav)*100)/100);
Tsw = 2*Tlanding1;
% Horizontal speed of the pendulum
HorSpeedFinal = StepLenght/Tsw;
SpringCompr = obstacle_height*0.5;
K = ceil( 2*(- m*grav*(SpringCompr + z_flight)+0.5*m*HorSpeedFinal^2)/SpringCompr^2);
omega = sqrt(K/m);
freq = omega/2/pi;
Tst = 1/freq*0.5; % stance time is half of the spring period
tt = 0.01; % integration time (sampling time t = 0.005 works)
Tst = ceil(Tst*100)/100; % round the landing instant to the upper centi-second (so that the counter is integer)
T1 = Tst + Tsw; % period corresponding to 1 step
Duty = Tst/T1;
Tliftoff1 = Tlanding1 + Tst;
T = 2*T1;
Tlanding2 = Tlanding1 + T1;
Tliftoff2 = Tliftoff1 + T1;
% number of control inputs
ni = 2;
% Declare model variables
theta = MX.sym('theta'); % theta
thetad = MX.sym('thetad'); % theta_dot
x = MX.sym('x');  % base x coordinate
xd = MX.sym('xd'); % base x dot
z = MX.sym('z'); % base y coordinate
zd = MX.sym('zd'); % base y dot

% state vector
state = [theta; thetad; x; xd; z; zd];
nv = size(state,1);
q = [theta; x; z];
dq = [thetad; xd; zd];
ddq = MX.sym('ddq',size(q,1));

u = MX.sym('u',ni);
% Control inputs
fn = u(1)*cos(theta);
ft = u(1)*sin(theta);
M = u(2);
% add force/torque limits
u1max = inf;
u1min = 0;
u2max = inf;
u2min = -inf;
% Model equations
v_sq = xd^2 + zd^2;
E = 0.5*m*len*v_sq + 0.5*m*len^2*thetad^2;
V = m*grav*z;

% Lagrangian
Lag = E - V;

% Equation of motion
% eq = jtimes(gradient(Lag,dq),q,dq) - gradient(Lag,q);
eq = jacobian(gradient(Lag,dq),q)*dq - gradient(Lag,q);
% xdot = [thetad; eq(1) - M - ft*len*cos(theta) - fn*len*sin(theta); xd; eq(2) - ft/m; zd; eq(3) + fn/m];
xdot = [thetad;  -grav/len*(theta) - M/I - fn*len*(theta)/I - ft*len/I; xd; eq(2) + ft/m; zd; eq(3) + fn/m];
% Objective term
L = u'*u;

% Continuous time dynamics
f = Function('f', {state, u}, {xdot, L});

% Control discretization
N = T/tt; % number of control intervals


X0 = MX.sym('X0', nv);
U = MX.sym('U',ni);
X = X0;
Q = 0;
% Runge Kutta Integrator
[F,X,Q] = RK4(f,X,U,T,N,Q,X0);
% Number of steps
NumIter = 31;
% initial position
XInit = 0;
Xfinal = StepLenght;
% Schedule of the bounding step by step
HorSpeedFinal =StepLenght/T;  % speed xdot at the start of the opt
HorSpeedInit = 0;  % speed xdot at the end of the opt
StateInit = [0;0;XInit;HorSpeedInit;z_flight+len;0];
StateFinal = [0;0;Xfinal;StepLenght/T;0;0];
CaptPoint = zeros(1,NumIter);

% initialize the vectors of the solution
% x1_opt = []; x2_opt = []; x3_opt = []; x4_opt = [];
% x5_opt = []; x6_opt = []; u1_opt = []; u2_opt = [];
for ii = 1:5:NumIter
    display(ii);
    % initialize the vectors of the solution
    x1_opt = []; x2_opt = []; x3_opt = []; x4_opt = [];
    x5_opt = []; x6_opt = []; u1_opt = []; u2_opt = [];
    CaptPoint(ii) = XInit + HorSpeedFinal/sqrt(-grav/len);
    
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
    lbw = [lbw;StateInit];
    ubw = [ubw;StateInit];
    w0 = [w0; 0; 0; 0; 0; 0; 0];
    
    % Formulate the NLP
    Xk = X0;
    
    for k=0:N-1
        % New NLP variable for the control
        Uk = MX.sym(['U_' num2str(k)], ni);
        w = {w{:}, Uk};
        lbw = [lbw; u1min; u2min];
        ubw = [ubw; u1max; u2max];
        w0 = [w0; 0; 0];
        
        % Integrate till the end of the interval
        [Xk_end, Jk] = easycall(F, Xk, Uk);
        J=J+Jk;
        % New NLP variable for state at end of interval
        Xk = MX.sym(['X_' num2str(k+1)], nv);
        w = {w{:}, Xk};
        if k == N-1
            lbw = [lbw; -inf; -inf;  -inf;  -inf; -inf; -inf];  % z coordinate can only be positive
            ubw = [ubw; inf; inf;  inf;  inf; inf; inf];
            w0 = [w0; 0; 0; 0; 0; 0; 0];
            J = J + Jfinal;
        else
            lbw = [lbw; -inf; -inf;  -inf;  -inf;  -inf;  -inf];  % z coordinate can only be positive
            ubw = [ubw;  inf;  inf;  inf;  inf;  inf;  inf];
            w0 = [w0; 0; 0; 0; 0; 0; 0];
        end
        % Add equality constraint
        g = {g{:}, Xk_end-Xk};
        lbg = [lbg; 0; 0; 0; 0; 0; 0];
        ubg = [ubg; 0; 0; 0; 0; 0; 0];
        if k ==0
            Xinit = Xk;
            g = {g{:}, Uk(1)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            
        elseif ((k >= floor(Tlanding1/tt))&&(k <= ceil(Tliftoff1/tt)))||((k >= floor(Tlanding2/tt))&&(k <= ceil(Tliftoff2/tt)))
            if (k == floor(Tlanding1/tt))||(k == floor(Tlanding2/tt))
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
            g = {g{:}, Xk(3) - l*sin(Xk(1)) - Xtouchd} ;
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            g = {g{:}, Xk(5) - l*cos(Xk(1)) - Ztouchd};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            % impose the torque M to be null during stance phase
            %             g = {g{:}, Uk(2)};
            %             lbg = [lbg; 0];
            %             ubg = [ubg; 0];
        elseif k == N-1-ii+1
            g = {g{:}, Xk - Xinit};
            lbg = [lbg; StateFinal];
            ubg = [ubg; StateFinal];
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
    % options.ipopt = struct('max_iter',100,'acceptable_tol',10e+100);
    solver = nlpsol('solver', 'ipopt', prob);
    
    % Solve the NLP
    arg = struct('x0', w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg);
    sol = solver(arg);
    w_opt = full(sol.x);
    % allocating decision variables for plotting
    x1_opt = [w_opt(1:8:end)];
    x2_opt = [w_opt(2:8:end)];
    x3_opt = [w_opt(3:8:end)];
    x4_opt = [w_opt(4:8:end)];
    x5_opt = [w_opt(5:8:end)];
    x6_opt = [w_opt(6:8:end)];
    u1_opt = [w_opt(7:8:end); nan];
    u2_opt = [w_opt(8:8:end); nan];
    
    % Optimized Capture Point
    OptCP(ii) = x3_opt(end) + x4_opt(end)/sqrt(-grav/len);
    
    
    % Plot the optimal solution
    Fn = [u1_opt].*cos(x1_opt);
    Ft = [u1_opt].*sin(x1_opt);
    M = u2_opt;
    
    Plot(x1_opt,x2_opt,x3_opt,x4_opt, x5_opt, x6_opt, Fn, Ft, M, tt, T,NumIter,len,StepLenght,Tlanding1,Tlanding2,Tliftoff1,Tliftoff2,OptCP,CaptPoint);
    % pause()
    % % Simulate the solution and plot
    % [Theta, Thetad, Xx, Xxd, Zz, Zzd] = Simulation(Fn,Ft,M,N,x1_opt(1),x2_opt(1),x3_opt(1),x4_opt(1),x5_opt(1),x6_opt(1),T);
    % Plot(Theta', Thetad', Xx', Xxd', Zz', Zzd', Fn, Ft, M, tt, T,NumIter,len,StepLenght,Tlanding1,Tlanding2,Tliftoff1,Tliftoff2,OptCP,CaptPoint);
   
    % Current state equal to the new initial state
    StateInit = [x1_opt(ii);x2_opt(ii);x3_opt(ii);x4_opt(ii);x5_opt(ii);x6_opt(ii)];
    StateFinal = [0;0;StepLenght;StepLenght/T;0;0];
    Tlanding1 = Tlanding1 - (ii+1)*tt
    Tlanding2 = Tlanding2 - (ii+1)*tt
    Tliftoff1 = Tliftoff1 - (ii+1)*tt
    Tliftoff2 = Tliftoff2 - (ii+1)*tt 
end