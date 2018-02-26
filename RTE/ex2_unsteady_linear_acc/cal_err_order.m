format long ;
%f = @(x) (1-exp(4*x))/2;
%f = @(x) (1-exp(-2*x))/2;
%f = @(x) x;
f = @(x,mu) mu^2*cos(pi*(x+0.1)).^4;

poly(1) = @(x) 1;
poly(2) = @(x) x;
poly(3) = @(x) (3*x.^2-1)/2;
poly(4) = @(x) (5*x.^3-3*x)/2;
poly(5) = @(x) (35*x.^4-30*x.^2+3)/8;
transform = @(x) 

u =[-0.960289856497536;
-0.796666477413627;
-0.525532409916329;
-0.183434642495650;
0.183434642495650;
0.525532409916329;
0.796666477413627;
0.960289856497536];
wl =[ 0.101228536290376;
0.222381034453374;
0.313706645877887;
0.362683783378362;
0.362683783378362;
0.313706645877888;
0.222381034453374;
0.101228536290376];

LG2 = [-sqrt(3)/3,sqrt(3)/3];
LG3 = [-sqrt(15)/5, 0, sqrt(15)/5];
LG4 = [-0.8611363115940520, -0.3399810435848560, 0.3399810435848560, 0.8611363115940520];
LG8 = [-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650, ...
       0.183434642495650,  0.525532409916329,  0.796666477413627, 0.960289856497536];

x0 = linspace(0, 1, 1000);
Nx = 80;
numer1 = load(['ex2_Nx', num2str(Nx),'_K',num2str(2),'_PP0.dat']);
x1 = numer1(:,1);
w1 = numer1(:,2);
y1 = numer1(:,3:end);
if(size(y1,2) == 2)
  hold on;
  plot(x0, f(x0,LG2(1)), '-k', x1, y1(:,1), 'ro');
  plot(x0, f(x0,LG2(2)), '-k', x1, y1(:,2), 'ro');
  min(min(y1))
  hold off;
elseif(size(y1,2) == 3)
  hold on;
  plot(x0, f(x0,LG3(1)), '-k', x1, y1(:,1), 'ro');
  plot(x0, f(x0,LG3(2)), '-k', x1, y1(:,2), 'ro');
  plot(x0, f(x0,LG3(3)), '-k', x1, y1(:,3), 'ro');
  min(min(y1))
  hold off;
elseif(size(y1,2) == 4)
  hold on;
  plot(x0, f(x0,LG4(1)), '-k', x1, y1(:,1), 'ro');
  plot(x0, f(x0,LG4(2)), '-k', x1, y1(:,2), 'ro');
  plot(x0, f(x0,LG4(3)), '-k', x1, y1(:,3), 'ro');
  plot(x0, f(x0,LG4(4)), '-k', x1, y1(:,4), 'ro');
  min(min(y1))
  hold off;
elseif(size(y1,2) == 8)
  hold on;
  plot(x0, f(x0,LG8(1)), '-k', x1, y1(:,1), 'ro');
  plot(x0, f(x0,LG8(2)), '-k', x1, y1(:,2), 'ro');
  plot(x0, f(x0,LG8(3)), '-k', x1, y1(:,3), 'ro');
  plot(x0, f(x0,LG8(4)), '-k', x1, y1(:,4), 'ro');
  plot(x0, f(x0,LG8(5)), '-k', x1, y1(:,5), 'ro');
  plot(x0, f(x0,LG8(6)), '-k', x1, y1(:,6), 'ro');
  plot(x0, f(x0,LG8(7)), '-k', x1, y1(:,7), 'ro');
  plot(x0, f(x0,LG8(8)), '-k', x1, y1(:,8), 'ro');
  hold off;
  %cal_err
  %err = [
  %f(x1,LG8(1)) - y1(:,1);
  %f(x1,LG8(2)) - y1(:,2);
  %f(x1,LG8(3)) - y1(:,3);
  %f(x1,LG8(4)) - y1(:,4);
  %f(x1,LG8(5)) - y1(:,5);
  %f(x1,LG8(6)) - y1(:,6);
  %f(x1,LG8(7)) - y1(:,7);
  %f(x1,LG8(8)) - y1(:,8)];
  %fprintf('%e, %e, %e\n',norm(err,2)/sqrt(length(err)),norm(err,1)/length(err),norm(err,'inf'));
  %cal_err
  err1 = ...
  dot(abs(f(x1,LG8(1)) - y1(:,1)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(2)) - y1(:,2)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(3)) - y1(:,3)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(4)) - y1(:,4)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(5)) - y1(:,5)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(6)) - y1(:,6)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(7)) - y1(:,7)), w1)*1/Nx/2 + ...
  dot(abs(f(x1,LG8(8)) - y1(:,8)), w1)*1/Nx/2;
  err2 = ...
  dot((f(x1,LG8(1)) - y1(:,1)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(2)) - y1(:,2)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(3)) - y1(:,3)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(4)) - y1(:,4)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(5)) - y1(:,5)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(6)) - y1(:,6)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(7)) - y1(:,7)).^2, w1)*1/Nx/2 + ...
  dot((f(x1,LG8(8)) - y1(:,8)).^2, w1)*1/Nx/2;
  errf = max([ ...
  abs(f(x1,LG8(1)) - y1(:,1));
  abs(f(x1,LG8(2)) - y1(:,2));
  abs(f(x1,LG8(3)) - y1(:,3));
  abs(f(x1,LG8(4)) - y1(:,4));
  abs(f(x1,LG8(5)) - y1(:,5));
  abs(f(x1,LG8(6)) - y1(:,6));
  abs(f(x1,LG8(7)) - y1(:,7));
  abs(f(x1,LG8(8)) - y1(:,8))]);


  %fprintf('%e, %e, %e\n',norm(err,2)/sqrt(Nx),norm(err,1)/Nx,norm(err,'inf'));
  fprintf('%e, %e, %e\n', err1, err2, errf);

end

%P4
err4 = [
];
err = err4;
%log2(err(1,:)./err(2,:))
%log2(err(2,:)./err(3,:))
%log2(err(3,:)./err(4,:))
%log2(err(4,:)./err(5,:))


