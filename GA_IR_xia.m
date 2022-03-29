%% 估计疫情曲线下降期的IR模型参数
%% 初始化程序
clc;clear;
%% 固定参数
global C;  %总人数
global t; %时间
global num; %遗传算法
global I;global di; global R;global pop;global no2;

C = xlsread('印度-疫情数据.xlsx','5'); %数据
pop = [32509617.71	422484047.8	36183.04367	0.62157152	0.064329007];
no2 = 15; %no2均值 美国 11 巴西6  印度15 英国9  俄罗斯14 德国10

LB = [0 0 0 0 -1 -1]; 
UB = [1 1 1 200 1 1];
% 英国的限制参数为
% LB = [0 0 0 0 0 0]; 
% UB = [2 2 2 200 1 1];

nvars = 6; %待估计的模型参数数量
num = 30; %给定遗传算法重复的次数
t = size(C,1);
%% 循环每个国家的参数
result1 = zeros(num,nvars+1);
%给各种群变量分配内存
I = zeros(1,t);di = zeros(1,t); R = zeros(1,t);
%定义待调用的两个子程序seaircal()和seair_constraint()
ObjectiveFunction = @seaircal_IR_xia;
ConstraintFunction = @seair_constraint;

%循环ga()遗传算法
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x是待估变量，fval是返回的自适应值
    result1(i,:) = [x,-fval];
end

%将计算结果导出为*.xls 格式的电子表格文档
result = sortrows(result1,nvars+1,'descend');  

[c, y] = IR_Country_xia(C,result(1,:),pop); %返回的y是就是每日新增感染人数，c就是截止到当天的所有人群数和一些参数，用于带入到下一步计算中pop