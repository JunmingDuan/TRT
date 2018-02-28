f1 = @(t,y) [0.5*y(2)^4-y(1); 2*y(1)-y(2)^4];
tspan = [0 1.4];
I0 = [8, 0.001];
[t,y] = ode23s(f1, tspan, I0);
plot(t, y);
legend('I', 'T');
plot(t, y(:,2));
axis([0.1,0.4,1.4,2]);
y(end,1)-0.5*y(end,2)^4
%[I T] = dsolve('DI==T^4/2-I', 'DT==2*I-T^4', 'I(0)==8', 'T(0)==0.001')

