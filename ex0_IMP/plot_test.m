format long ;
hold on;
f1 = @(t,y) [0.5*y(2)^4-y(1); 2*y(1)-y(2)^4];
tspan = [0:0.001:0.4];
I0 = [8, 0.001];
[t,y] = ode23s(f1, tspan, I0);
t = t(2:end,:);
y = y(2:end,:);
plot(t,y(:,2));
DAT1 = load('IMP_P0_dt0.001.dat');
DAT2 = load('IMP_P0_dt0.002.dat');
DAT3 = load('IMP_P0_dt0.004.dat');
DAT4 = load('IMP_P0_dt0.008.dat');
DAT5 = load('IMP_P0_dt0.016.dat');
DAT6 = load('IMP_P0_dt0.400.dat');
plot(DAT1(:,1), DAT1(:,2), '-');%'o'
plot(DAT2(:,1), DAT2(:,2), '-');%'<'
plot(DAT3(:,1), DAT3(:,2), '-');%'>'
plot(DAT4(:,1), DAT4(:,2), '-');%'^'
plot(DAT5(:,1), DAT5(:,2), '-');%'v'
plot(DAT6(:,1), DAT6(:,2), '-');%'*'
axis([0.1,0.4,1.4,2]);
legend('T', 'dt=0.001', 'dt=0.002', 'dt=0.004', 'dt=0.008', 'dt=0.016', 'dt=0.1');
%err(1) = norm(DAT1(:,2)-y(1:1:end,2), 2);
%err(2) = norm(DAT2(:,2)-y(1:2:end,2), 2);
%err(3) = norm(DAT3(:,2)-y(1:4:end,2), 2);
%err(4) = norm(DAT4(:,2)-y(1:8:end,2), 2);
%err(5) = norm(DAT5(:,2)-y(1:16:end,2), 2);
%err(6) = norm(DAT6(:,2)-y(1:200:end,2), 2);
%err
%[I T] = dsolve('DI==T^4/2-I', 'DT==2*I-T^4', 'I(0)==8', 'T(0)==0.001')

