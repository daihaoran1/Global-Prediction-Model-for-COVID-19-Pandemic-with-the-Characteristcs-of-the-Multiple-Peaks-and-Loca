%% 程序功能：为遗传算法提供自适应函数（SEAIIDR模型），返回simerr作为评价指标，即每次进化会选出最小的simerr值
function simerr = seaircal_IR(a)
%% 定义参数
global C;  %总人数
global t; %时间
global I;global di; global R;global pop;
I(1) = pop(1); %有症状感染人群数（未检测到的）
R(1) = pop(2); %康复人群数
di(1) = pop(3);
%% 循环以定步长1天求解传染病模型
for i = 1:1:t-1
    lambda = a(1) + a(2)/(1+exp(a(3)*(a(4)-i))); 
    di(i+1) = lambda*I(i)- a(5)*I(i); 
    I(i+1) = I(i)+lambda*I(i) - a(5)*I(i); 
    R(i+1) = R(i)+a(5)*I(i);
end
y = di';
%% 计算R2值
%计算误差评定指数：模型参数条件下SARS 传播过程与真实过程的拟合度值
simerr = -FITNESS(C(:,4),y); %因为ga()函数默认输出的目标函数的结果是最小值，但是R2的结果要最大才对，所以加负号