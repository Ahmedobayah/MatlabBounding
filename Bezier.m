close all

% explicit cubic Bezier curve:
Tstance = 0:0.01:1;
t = Tstance(1:floor(end/2));
a0 = 0;
a1 = 0.1;
a2 = 0.1;
a3 = 0;
Bez = (1-t).^3.*a0 + 3*(1-t).^2.*t*a1 + (1-t).*t.^2*a2 + t.^3*a3;
plot(t,Bez)