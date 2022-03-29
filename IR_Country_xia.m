%% 程序功能：将估计的四个国家参数带入自己国家模拟个国家疫情初始状态
function [cc, y] = IR_Country_xia(C,a,pop) 
t = size(C,1); 
%给各人群变量分配内存
global no2;
I = zeros(1,t);di = zeros(1,t);R = zeros(1,t);lambda = zeros(1,t);
I(1) = pop(1); %有症状感染人群数（未检测到的）
R(1) = pop(2); %康复人群数
di(1) = pop(3);
beta0 = a(1); %传染率
Yi = pop(5); %康复率
beta1 = a(2);
for i = 1:1:t-1
    lambda(i) = beta0 + beta1/(1+exp(a(3)*(i-a(4))));
    eta = a(5)+a(6)*(no2-C(i,1)); %no2均值 美国 11 巴西  印度14 英国8  俄罗斯
    di(i+1) = lambda(i)*I(i)-eta*I(i) - Yi*I(i);
    I(i+1) = I(i)+lambda(i)*I(i)-eta*I(i) - Yi*I(i);
    R(i+1) = R(i)+Yi*I(i);
end
y = [di' lambda'];
cc = [I(end)  R(end)  di(end) a(1) pop(5)];
end