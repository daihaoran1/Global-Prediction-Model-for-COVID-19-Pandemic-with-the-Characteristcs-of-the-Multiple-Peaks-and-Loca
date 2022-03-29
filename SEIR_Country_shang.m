%% 程序功能：将估计的四个国家参数带入自己国家模拟个国家疫情初始状态
function [cc, y] = SEIR_Country_shang(C,a,pop) 
%% 固定参数
global di;global did;global sigma;global v;global Ya;global Yi;global Yid;global no2;
%% 将各各国估计参数带入自己国家
t = size(C,1); 
%给各人群变量分配内存
S = zeros(1,t); E = zeros(1,t);A = zeros(1,t);I = zeros(1,t);R = zeros(1,t);ID = zeros(1,t); EdConcount = zeros(1,t); lambda = zeros(1,t);
         
E(1) = pop(2); %潜伏者人群数
A(1) = pop(3); %无症状感染人群数（未检测到的）
I(1) = pop(4); %有症状感染人群数（未检测到的）
ID(1) = pop(5); %通过检测具有感染性的人群数，即正在医院接受治疗的（包括了无症状感染人群和有症状感染人群）
R(1) = pop(6); %康复人群数
S(1) = pop(1)+pop(6); %易感染人群数
EdConcount(1) = pop(7); %每天新增确诊人数
beta0 = a(1); %传染率
beta1 = a(2);
theta = pop(8); %无症状感染者的检测率
fa = pop(9); %有症状感染者的检测率
for i = 1:1:t-1
    lambda(i) = beta0 + beta1/(1+exp(a(3)*(a(4)-i)));  %易感染人群变为潜伏人群的概率
    eta =  a(5)+a(6)*(no2-C(i,1));
    S(i+1) = S(i)-lambda(i)*S(i)+eta*I(i)+eta*A(i); 
    E(i+1) = E(i)+lambda(i)*S(i)-sigma*E(i);
    A(i+1) = A(i)+v*sigma*E(i)-(Ya)*A(i)-eta*A(i);
    I(i+1) = I(i)+(1-v)*sigma*E(i)-(Yi+di)*I(i)-eta*I(i);
    ID(i+1) = ID(i)+theta*A(i)+fa*I(i)-(Yid+did)*ID(i);
    R(i+1) = R(i)+Yid*ID(i)+Yi*I(i)+Ya*A(i);
    EdConcount(i+1) = theta*A(i)+fa*I(i);
end
y = [EdConcount' lambda'];
cc = [S(end) E(end) A(end) I(end) ID(end) R(end) EdConcount(end) theta fa];
end