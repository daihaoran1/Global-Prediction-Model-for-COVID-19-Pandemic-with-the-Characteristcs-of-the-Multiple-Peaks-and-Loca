%% 程序功能：为遗传算法提供自适应函数（SEAIIDR模型），返回simerr作为评价指标，即每次进化会选出最小的simerr值
function simerr = seaircal_IR_shang(a)
%% 定义参数
global C;  %总人数
global t; %时间
global I;global di; global R;global pop;global no2;

I(1) = pop(1); %有症状感染人群数（未检测到的）
R(1) = pop(2); %康复人群数
di(1) = pop(3);
beta0 = a(1); %传染率
Yi = pop(5); %康复率
beta1 = a(2);
%% 循环以定步长1天求解传染病模型
for i = 1:1:t-1
    lambda = beta0 + beta1/(1+exp(a(3)*(a(4)-i))); 
    eta = a(5)+a(6)*(no2-C(i,1)); %no2均值 美国 11 巴西  印度14 英国8  俄罗斯
    di(i+1) = lambda*I(i)-eta*I(i) - Yi*I(i);
    I(i+1) = I(i)+lambda*I(i)-eta*I(i) - Yi*I(i);
    R(i+1) = R(i)+Yi*I(i);
end
y = di';
%% 计算R2值
%计算误差评定指数：模型参数条件下SARS 传播过程与真实过程的拟合度值
simerr = -FITNESS(C(:,4),y); %因为ga()函数默认输出的目标函数的结果是最小值，但是R2的结果要最大才对，所以加负号