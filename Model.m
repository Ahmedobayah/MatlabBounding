function xdot_real = Model(X,U,Thetadd,Xxdd,Zzdd)
grav = -9.81;
m = 1;
len = 0.7;
I = m*len^2;

Theta = X(1); Thetad = X(2); Xx = X(3); Xxd = X(4); Zz = X(5); Zzd = X(6);
Fn = U(1); Ft = U(2); M = U(3);

% xdot_real = [Thetad;  -grav/len*sin(Theta) - M - Fn*len*sin(Theta) - Ft*len*cos(Theta); Xxd; Xxdd + Ft/m; Zzd; Zzdd + m*grav + Fn/m];
xdot_real = [Thetad;  Thetadd - M/I - Fn*len*sin(Theta)/I - Ft*len*cos(Theta)/I; Xxd; Xxdd + Ft/m; Zzd; Zzdd + m*grav + Fn/m];
end