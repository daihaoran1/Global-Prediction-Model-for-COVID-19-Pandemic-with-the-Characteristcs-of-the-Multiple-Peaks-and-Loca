%% 程序的功能：利用遗传算法估计公布案例开始时实际中已存在的潜伏者人群数E(1)、无症状感染人群数（未检测）A(1)和
clc;clear;
%% 固定参数
global C; global t; global S; global E; global A; global I; global ID; global R; global EdConcount; global di;
global did; global sigma; global v; global Ya; global Yi; global Yid; global Nh; global pop;

di = 0.015; %有症状感染者的死亡率
did = 0.015; %通过检测具有感染性的人群的死亡率
sigma = 1/5.2; %从潜伏状态到感染状态（包括又症状和无症状感染者）的概率
v = 0.5; %无症状感染者的比例，1-v就是有症状感染者的比例
Ya = 0.13978; %无症状感染者的康复率
Yi = 0.13978; %有症状感染者的康复率
Yid = 1/15; %通过检测具有感染性的人群的康复率
nvars = 6; %待估计的模型参数数量
num = 30; %给定遗传算法重复的次数

C = xlsread('印度-疫情数据.xlsx','0'); %数据
Nh = 1354051854; %总人口  美国2020年人口326766748；英国66573504；德国82293457；印度1354051854；俄罗斯人口143964709；巴西210867954
pop = [0, 0, 1, 0, 0, 1]; % E A I ID R EDveryday

LB = [0 0 0 0 0 0.01]; 
UB = [1 1 1 200 0.0001 0.1];
t = size(C,1); result1 = zeros(num,nvars+1);S = zeros(1,t); E = zeros(1,t); A = zeros(1,t);I = zeros(1,t); ID = zeros(1,t);R = zeros(1,t); EdConcount = zeros(1,t);
ObjectiveFunction = @seaircal_SEIR; FUN = @seair_constraint;
%循环ga()遗传算法
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,FUN); %x是待估变量，fval是返回的自适应值
    result1(i,:) = [x,-fval];
end
result = sortrows(result1,nvars+1,'descend');     

[c, y] = SEIR_Country(C, result(1,:),pop);