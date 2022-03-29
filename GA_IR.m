%% 求解第一阶段无人为干预措施的IR参数
%% 初始化程序
clc;clear;
%% 固定参数
global C;  %总人数
global t; %时间
global num; %遗传算法
global I;global di; global R;global pop;

C = xlsread('德国-疫情数据.xlsx','0'); %只需改动数据来源即可
pop = [1 0 1];
LB = [0 0 0 0 0]; 
UB = [1 1 5 100 0.1];
nvars = 5; %待估计的模型参数数量
num = 50; %给定遗传算法重复的次数
t = size(C,1);
result1 = zeros(num,nvars+1);
I = zeros(1,t);di = zeros(1,t); R = zeros(1,t);

%定义待调用的两个子程序seaircal()和seair_constraint()
ObjectiveFunction = @seaircal_IR;
ConstraintFunction = @seair_constraint;
    
%循环ga()遗传算法
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x是待估变量，fval是返回的自适应值
    result1(i,:) = [x,-fval];
end
%将计算结果导出为*.xls 格式的电子表格文档
result = sortrows(result1,nvars+1,'descend');  

[c, y] = IR_Country(C,result(1,:),pop); %返回的y是就是每日新增感染人数，c就是截止到当天的所有人群数和一些参数，用于带入到下一步计算中p