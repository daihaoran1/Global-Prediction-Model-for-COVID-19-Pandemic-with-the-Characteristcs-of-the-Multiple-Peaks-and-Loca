%% ����Ĺ��ܣ������Ŵ��㷨���ƹ���������ʼʱʵ�����Ѵ��ڵ�Ǳ������Ⱥ��E(1)����֢״��Ⱦ��Ⱥ����δ��⣩A(1)��
%% ��ʼ������
clc;clear;
%% Ҫ�Ķ��Ĳ���
%�޶������Ƶ�ģ�Ͳ�����Լ���ռ䣬LBΪ���ޣ�UBΪ����
LB = [0 0 0 0 -1 -1]; 
UB = [1 1 1 100 1 1];
%% �̶�����
global C;  %������
global t; %ʱ��
global num; %�Ŵ��㷨
global S; global E; global A; global I; global ID; global R; global EdConcount; global di; global pop;global no2;
global did; global sigma; global v; global alfa; global Ya; global Yi; global Yid; global Nh; %���˿�

di = 0.015; %��֢״��Ⱦ�ߵ�������
did = 0.015; %ͨ�������и�Ⱦ�Ե���Ⱥ��������
sigma = 1/5.2; %��Ǳ��״̬����Ⱦ״̬��������֢״����֢״��Ⱦ�ߣ��ĸ���
v = 0.5; %��֢״��Ⱦ�ߵı�����1-v������֢״��Ⱦ�ߵı���
alfa = 0.5; %��֢״��Ⱦ������֢״��Ⱦ����ȣ����Ⱦ�ʵ��޸Ĳ���
Ya = 0.13978; %��֢״��Ⱦ�ߵĿ�����
Yi = 0.13978; %��֢״��Ⱦ�ߵĿ�����
Yid = 1/15; %ͨ�������и�Ⱦ�Ե���Ⱥ�Ŀ�����                   ��������ͬ
nvars = 6; %�����Ƶ�ģ�Ͳ�������[p1 p2 p3 p4 eta1 eta2]
num = 30; %�����Ŵ��㷨�ظ��Ĵ���
%% ѭ��ÿ�����ҵĲ���

result1 = zeros(num,nvars+1);
Nh = 210867954; %���˿�  ����2020���˿�326766748��Ӣ��66573504���¹�82293457��������46397452��ӡ��1354051854������˹�˿�143964709��ī����130759074������210867954��

%������Ⱥ���������ڴ�
S = zeros(1,t);E = zeros(1,t);A = zeros(1,t);I = zeros(1,t);ID = zeros(1,t);R = zeros(1,t);EdConcount = zeros(1,t);

C = xlsread('ӡ��-��������.xlsx','4'); %����   ������һ��
pop = [1172773633	113133128.4	8268818.761	8174891.082	4045426.588	123467923.5	450283.3251	9.98E-05	0.055737305];
no2 = 15;  %no2��ֵ ���� 11 ����6  ӡ��15 Ӣ��9  ����˹14 �¹�10

t = size(C,1);

%��������õ������ӳ���seaircal()��seair_constraint()
ObjectiveFunction = @seaircal_SEIR_xia;
ConstraintFunction = @seair_constraint;

%ѭ��ga()�Ŵ��㷨
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x�Ǵ���������fval�Ƿ��ص�����Ӧֵ
    result1(i,:) = [x,-fval];
end

%������������Ϊ*.xls ��ʽ�ĵ��ӱ���ĵ�
result = sortrows(result1,nvars+1,'descend');  
[c, y] = SEIR_Country_xia(C,result(1,:),pop); %���ص�y�Ǿ���ÿ��������Ⱦ������c���ǽ�ֹ�������������Ⱥ����һЩ���������ڴ��뵽��һ��������pop   
