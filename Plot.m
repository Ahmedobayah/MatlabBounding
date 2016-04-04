%% Plot the solution
function [] = Plot(x1_opt,x2_opt,x3_opt,x4_opt, x5_opt, x6_opt, Fn, Ft, M,tt,T,StepNum,len,StepLenght,Tlanding1,Tlanding2,Tliftoff1,Tliftoff2,OptCP,CaptPoint)

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

% handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan],[u3_opt; nan]]);hold on;
stairs(tgrid, Fn,'r'), hold on;
stairs(tgrid, Ft,'b'),
stairs(tgrid, [M],'k')
% stairs(tgrid, [u3_opt;nan],'g')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('fn [N]', 'ft [N]', 'M [Nm]');
xlabel('time [s]');
hold off;

figure,
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
    if (k>=floor(Tlanding1/tt))&&(k <= ceil(Tliftoff1/tt))||(k>=floor(Tlanding2/tt))&&(k <= ceil(Tliftoff2/tt)) % stance phase
        if (k == floor(Tlanding1/tt))
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
    title('L = \theta^2 + u^2');
    drawnow
    pause(0.05)
end

figure,
subplot(3,1,1),plot(tgrid,Fn),legend('normal GRF'), title('Ground Reaction Forces and Torques');
subplot(3,1,2),plot(tgrid,Ft), legend('tangent GRF');
subplot(3,1,3),plot(tgrid,M),legend('hip torque M');
