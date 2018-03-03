function plot_ex1(K, n, limiter);

format long;
hold on;

numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K) ,'_PP',num2str(limiter),'.dat']);
numer2 = load('reflex.dat');
x1 = numer1(:,1); y1 = numer1(:,3);
x2 = numer2(:,1); y2 = numer2(:,3);
plot(x1, y1, 'r-o');
plot(x2, y2, 'b-o');
%axis([0,0.02,0,1]);
%axis([0,0.055,0,1]);

