%% 程序功能：为遗传算法提供自适应函数（SEAIIDR模型），返回simerr作为评价指标，即每次进化会选出最小的simerr值
function simerr = seaircal_SEIR(a)
%% 定义参数
global C; global t; global S; global E; global A; global I; global ID; global R; global EdConcount;global Nh; global di;
global did; global sigma; global v; global Ya; global Yi; global Yid; global pop;
E(1) = pop(1); %潜伏者人群数
A(1) = pop(2); %无症状感染人群数（未检测到的）
I(1) = pop(3); %有症状感染人群数（未检测到的）
ID(1) = pop(4); %通过检测具有感染性的人群数，即正在医院接受治疗的（包括了无症状感染人群和有症状感染人群）
R(1) = pop(5); %康复人群数
S(1) = Nh - pop(1) - pop(2) - pop(3)- pop(4) -pop(5); %易感染人群数
EdConcount(1) = pop(6); %每天新增确诊人数
theta = a(5); %无症状感染者的检测率
fa = a(6); %有症状感染者的检测率
%% 循环以定步长1天求解传染病模型
for i = 1:1:t-1
    lambda = a(1) + a(2)/(1+exp(a(3)*(a(4)-i))); %易感染人群变为潜伏人群的概率
    S(i+1) = S(i)-lambda*S(i); 
    E(i+1) = E(i)+lambda*S(i)-sigma*E(i);
    A(i+1) = A(i)+v*sigma*E(i)-(Ya)*A(i);
    I(i+1) = I(i)+(1-v)*sigma*E(i)-(Yi+di)*I(i);
    ID(i+1) = ID(i)+theta*A(i)+fa*I(i)-(Yid+did)*ID(i);
    R(i+1) = R(i)+Yid*ID(i)+Yi*I(i)+Ya*A(i);
    EdConcount(i+1) = theta*A(i)+fa*I(i);
end
y = EdConcount';
%% 计算R2值
%计算误差评定指数：模型参数条件下SARS 传播过程与真实过程的拟合度值
simerr = -FITNESS(C(:,4),y); %因为ga()函数默认输出的目标函数的结果是最小值，但是R2的结果要最大才对，所以加负号