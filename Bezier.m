close all

% explicit cubic Bezier curve:
t = 0:0.01:1;
m = 1;
grav = -9.81;
% t = Tstance(1:floor(end/2))';
Tstance = 0.15;
T = 2;

a0 = [0];
a1 = [1];
a2 = [-1];
a3 = [0];
alpha = 2;
Bez = alpha*((1-t).^3*a0 + 3*(1-t).^2.*t*a1 + (1-t).*t.^2*a2 + t.^3*a3);
figure, stairs(t,Bez);hold on;
% t2 = 1:0.01:2;
% a0 = [0];
% a1 = [0.5];
% a2 = [1];
% a3 = [0];
% Bez2 = (1-t2).^3*a0 + 3*(1-t2).^2.*t2*a1 + (1-t2).*t2.^2*a2 + t2.^3*a3;
% stairs(t2,Bez2,'k')