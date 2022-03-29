%% 程序的功能：利用遗传算法估计公布案例开始时实际中已存在的潜伏者人群数E(1)、无症状感染人群数（未检测）A(1)和
%% 初始化程序
clc;clear;
%% 要改动的参数
%限定待估计的模型参数的约束空间，LB为下限，UB为上限
LB = [0 0 0 0 -1 -1]; 
UB = [1 1 1 200 1 1];
%% 固定参数
global C;  %总人数
global t; %时间
global num; %遗传算法
global S; global E; global A; global I; global ID; global R; global EdConcount; global di; global pop;global no2;
global did; global sigma; global v; global alfa; global Ya; global Yi; global Yid; global Nh; %总人口

di = 0.015; %有症状感染者的死亡率
did = 0.015; %通过检测具有感染性的人群的死亡率
sigma = 1/5.2; %从潜伏状态到感染状态（包括又症状和无症状感染者）的概率
v = 0.5; %无症状感染者的比例，1-v就是有症状感染者的比例
alfa = 0.5; %无症状感染者与有症状感染者相比，其感染率的修改参数
Ya = 0.13978; %无症状感染者的康复率
Yi = 0.13978; %有症状感染者的康复率
Yid = 1/15; %通过检测具有感染性的人群的康复率                   各地区不同
nvars = 6; %待估计的模型参数数量[p1 p2 p3 p4 eta1 eta2]
num = 30; %给定遗传算法重复的次数
%% 循环每个国家的参数

result1 = zeros(num,nvars+1);
Nh = 66573504; %总人口  美国2020年人口326766748；英国66573504；德国82293457；西班牙46397452；印度1354051854；俄罗斯人口143964709；墨西哥130759074；巴西210867954；

%给各种群变量分配内存
S = zeros(1,t);E = zeros(1,t);A = zeros(1,t);I = zeros(1,t);ID = zeros(1,t);R = zeros(1,t);EdConcount = zeros(1,t);

C = xlsread('英国-疫情数据.xlsx','1'); %数据   重新跑一遍
pop = [62014020.24	1640283.417	658798.9202	621293.875	32748.6208	1575926.318	5712.95119	2.47E-07	0.01013113];
no2 = 9;   %no2均值 美国 11 巴西6  印度15 英国9  俄罗斯14 德国10
t = size(C,1);

%定义待调用的两个子程序seaircal()和seair_constraint()
ObjectiveFunction = @seaircal_SEIR_shang;
ConstraintFunction = @seair_constraint;

%循环ga()遗传算法
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x是待估变量，fval是返回的自适应值
    result1(i,:) = [x,-fval];
end

%将计算结果导出为*.xls 格式的电子表格文档
result = sortrows(result1,nvars+1,'descend');  
[c, y] = SEIR_Country_shang(C,result(1,:),pop); %返回的y是就是每日新增感染人数，c就是截止到当天的所有人群数和一些参数，用于带入到下一步计算中pop   
