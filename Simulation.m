function [Theta, Thetad, Xx, Xxd, Zz, Zzd] = Simulation(Fn,Ft,M,N,x1_opt,x2_opt,x3_opt,x4_opt,x5_opt,x6_opt,T)
% Simulation on the real robot
Theta = zeros(1,N); Thetad = zeros(1,N); Thetadd = zeros(1,N);
Xx = zeros(1,N); Xxd = zeros(1,N); Xxdd = zeros(1,N);
Zz = zeros(1,N); Zzd = zeros(1,N); Zzdd = zeros(1,N);
Theta(1) = x1_opt; Thetad(1) = x2_opt; Xx(1) = x3_opt; Xxd(1) = x4_opt; Zz(1) = x5_opt; Zzd(1) = x6_opt;

IntStep = 4; % RK4 steps per interval
tt = T/N;
DT = T/N/IntStep;
for k = 1:N
   % Runge Kutta 4
   X(:,k) = [Theta(k);Thetad(k);Xx(k);Xxd(k);Zz(k);Zzd(k)];
   U = [Fn(k),Ft(k),M(k)];
   if k > 1
        thetadd = X(1,k-1);
        xdd = X(3,k-1);
        zdd = X(3,k-1);
   else
       thetadd = 0;
       xdd = 0;
       zdd = 0;
   end
   
   for ii = 1:IntStep
        k1 = Model(X(:,k),U,thetadd,xdd,zdd);
        [k2] = Model(X(:,k)+DT/2*k1,U,thetadd,xdd,zdd);
        [k3] = Model(X(:,k)+DT/2*k2,U,thetadd,xdd,zdd);
        [k4] = Model(X(:,k)+DT*k3,U,thetadd,xdd,zdd);
        X(:,k)=X(:,k)+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    Theta(k+1) = X(1,k);
    Thetad(k+1) = X(2,k);
    Xx(k+1) = X(3,k);
    Xxd(k+1) = X(4,k);
    Zz(k+1) = X(5,k);
    Zzd(k+1) = X(6,k);
   
%     if k>1
%         Thetadd(k+1) = -(X(2,k)-X(2,k-1))/tt;
%         Xxdd(k+1) = -(X(4,k)-X(4,k-1))/tt;
%         Zzdd(k+1) = -(X(6,k)-X(6,k-1))/tt;
%     end

    
end