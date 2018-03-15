function plot_ex2(K, n, limiter, minmod);

format long;
hold on;

numer0 = load('ref/Marshak_NT_Nx1000_P0_0.23_dt5.75e-5.dat');
numer1 = load(['ex2_Nx',num2str(n),'_K',num2str(K) ,'_PP',num2str(limiter), '_MD', num2str(minmod), '.dat']);
%numer2 = load('NT1.dat');
%numer3 = load('NTF.dat');
x0 = numer0(:,1); y0 = numer0(:,3);
x1 = numer1(:,1); y1 = numer1(:,3);
%x2 = numer2(:,1); y2 = numer2(:,3);
%x3 = numer3(:,1); y3 = numer3(:,3);
plot(x0, y0, '-k');
plot(x1, y1, '-o');
%plot(x2, y2, '-o');
%plot(x3, y3, '-*');
%plot(x2, y2, '--');
legend('ref0.23', 'Num');

