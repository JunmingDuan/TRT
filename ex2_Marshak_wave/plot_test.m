function plot_ex2(K, n, limiter, minmod);

format long;
hold on;

numer1 = load(['ex2_Nx',num2str(n),'_K',num2str(K) ,'_PP',num2str(limiter), '_MD', num2str(minmod), '.dat']);
size(numer1)
x1 = numer1(:,1); y1 = numer1(:,3);
plot(x1, y1, 'r-o');
%axis([0,0.02,0,1]);

