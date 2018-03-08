%/**
%* @file quartic_eqn.m
%* @brief Solve A*T^4 + T = C
%* @author Duan Junming, duanjm@pku.edu.cn
%* @version 1.0
%* @date 2018-03-08
%*/
function T = quartic_eqn(A, C)
hold on;
x0 = linspace(-2, 1, 1e3);
plot(x0, A*x0.^4+x0-C, '-r', x0, zeros(size(x0)), '-k');
plot(x0, -A*x0.^4+C, '-b');
str = '$$ \hat{u} $$';
text(0,2, str, 'Interpreter', 'latex');
xlabel(str, 'Interpreter', 'latex');
T0 = roots([A,0,0,1,-C]);
label = find(isreal(T0));
for i = 1:length(label);
  if(label(i)==1 && T>0 )
    T = T0(label);
    return;
  end
end
end

