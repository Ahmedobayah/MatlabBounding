function [Theta, Thetad, Xx, Xxd, Zz, Zzd] = Simulation(Fn,Ft,M,tt,len,m,grav,N,x1_opt,x2_opt,x3_opt,x4_opt,x5_opt,x6_opt)
% Simulation on the real robot
Theta = zeros(1,N); Thetad = zeros(1,N); Thetadd = zeros(1,N);
Xx = zeros(1,N); Xxd = zeros(1,N); Xxdd = zeros(1,N);
Zz = zeros(1,N); Zzd = zeros(1,N); Zzdd = zeros(1,N);

Xx(1) = x3_opt(1); Xxd(1) = x4_opt(1); Zz(1) = x5_opt(1);
xnew = [x1_opt(1); x2_opt(2); x3_opt(3); x4_opt(1); x5_opt(1); x6_opt(1)];
for k = 1:N
    xdot_real = [Thetad(k);  -grav/len*sin(Theta(k)) - M(k) - Fn(k)*len*sin(Theta(k)) - Ft(k)*len*cos(Theta(k)); Xxd(k); Xxdd(k) + Ft(k)/m; Zzd(k); Zzdd(k) + m*grav + Fn(k)/m];
    xnew = [xnew, xnew(:,end)+xdot_real.*tt];

    Theta(k+1) = xnew(1,k);
    Thetad(k+1) = xnew(2,k);
    Xx(k+1) = xnew(3,k);
    Xxd(k+1) = xnew(4,k);
    Zz(k+1) = xnew(5,k);
    Zzd(k+1) = xnew(6,k);
   
    if k>1
        Thetadd(k+1) = -(xnew(2,k)-xnew(2,k-1))/tt;
        Xxdd(k+1) = -(xnew(4,k)-xnew(4,k-1))/tt;
        Zzdd(k+1) = -(xnew(6,k)-xnew(6,k-1))/tt;
    else
        Thetadd(k+1) = 0;
        Xxdd(k+1) = 0;
        Zzdd(k+1) = 0;
    end
    
end