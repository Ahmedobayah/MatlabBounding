% Plot the solution
close all;

N = 150;
tt = 0.01;
T = tt*N;
len = 1;

x1_opt = w_opt(1:9:end);
x2_opt = w_opt(2:9:end);
x3_opt = w_opt(3:9:end);
x4_opt = w_opt(4:9:end);
x5_opt = w_opt(5:9:end);
x6_opt = w_opt(6:9:end);
u1_opt = w_opt(7:9:end);
u2_opt = w_opt(8:9:end);
u3_opt = w_opt(9:9:end);

tgrid = linspace(0, T, N+1);
clf;
subplot(2,1,1),hold on
plot(tgrid, x1_opt, 'k--')
% plot(tgrid, x2_opt, 'r--')
% plot(tgrid, x3_opt, 'k^')
% plot(tgrid, x4_opt, 'r*')
plot(tgrid, x5_opt - len*cos(x1_opt), 'g')
plot(tgrid, x5_opt, 'b')

% stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('theta','z base','z')
subplot(2,1,2);
% handle = stairs(tgrid,[[u1_opt; nan],[u2_opt; nan],[u3_opt; nan]]);hold on;
stairs(tgrid, [u1_opt;nan],'r'), hold on;
stairs(tgrid, [u2_opt;nan],'k')
stairs(tgrid, [u3_opt;nan],'g')
% handle(1).Marker = 'o'; handle(2).Marker = '*';
legend('fn', 'ft', 'M');
hold off;

figure(2)
n = size(x1_opt,1);
for k=1:n
    P0z = x5_opt(k) - len*cos(x1_opt(k));    
    P0x = x3_opt(k) - len*sin(x1_opt(k)); 
    % end effector coordinates
    zz = [x5_opt(k), P0z];
    xx = [x3_opt(k), P0x];
    % base coordinates
    floatz = [0, 0, P0z];
    floatx = [0, P0x, P0x];
    floorx = [-1 1];
    floory = [0 0];
    plot(xx,zz,'k', floatx, floatz, 'r', x3_opt(k), x5_opt(k),'ro', floorx, floory, 'k--')
%     axis([-2 2 -0.5 5.5])
    % Store the frame
    drawnow
    pause(0.1)
end


figure(3)
for k =1:n 
plot(x1_opt,x2_opt,'b',x1_opt(k),x2_opt(k),'r*')
xlabel('theta'),ylabel('theta dot')
drawnow
pause(0.05)
end