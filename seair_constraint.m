%% 程序功能：条件约束函数，限制初始种群的选值（A+I<E,A<I）
function [c, ceq] = seair_constraint(x)
c = [x(1)+x(2)-1]; %c是小于0，ceq是等于0，即-x(1)+x(2)+x(3)<0，,x(2)-x(3)<0
ceq = [];