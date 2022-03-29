%% ����Ĺ��ܣ������Ŵ��㷨���ƹ���������ʼʱʵ�����Ѵ��ڵ�Ǳ������Ⱥ��E(1)����֢״��Ⱦ��Ⱥ����δ��⣩A(1)��
clc;clear;
%% �̶�����
global C; global t; global S; global E; global A; global I; global ID; global R; global EdConcount; global di;
global did; global sigma; global v; global Ya; global Yi; global Yid; global Nh; global pop;

di = 0.015; %��֢״��Ⱦ�ߵ�������
did = 0.015; %ͨ�������и�Ⱦ�Ե���Ⱥ��������
sigma = 1/5.2; %��Ǳ��״̬����Ⱦ״̬��������֢״����֢״��Ⱦ�ߣ��ĸ���
v = 0.5; %��֢״��Ⱦ�ߵı�����1-v������֢״��Ⱦ�ߵı���
Ya = 0.13978; %��֢״��Ⱦ�ߵĿ�����
Yi = 0.13978; %��֢״��Ⱦ�ߵĿ�����
Yid = 1/15; %ͨ�������и�Ⱦ�Ե���Ⱥ�Ŀ�����
nvars = 6; %�����Ƶ�ģ�Ͳ�������
num = 30; %�����Ŵ��㷨�ظ��Ĵ���

C = xlsread('ӡ��-��������.xlsx','0'); %����
Nh = 1354051854; %���˿�  ����2020���˿�326766748��Ӣ��66573504���¹�82293457��ӡ��1354051854������˹�˿�143964709������210867954
pop = [0, 0, 1, 0, 0, 1]; % E A I ID R EDveryday

LB = [0 0 0 0 0 0.01]; 
UB = [1 1 1 200 0.0001 0.1];
t = size(C,1); result1 = zeros(num,nvars+1);S = zeros(1,t); E = zeros(1,t); A = zeros(1,t);I = zeros(1,t); ID = zeros(1,t);R = zeros(1,t); EdConcount = zeros(1,t);
ObjectiveFunction = @seaircal_SEIR; FUN = @seair_constraint;
%ѭ��ga()�Ŵ��㷨
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,FUN); %x�Ǵ���������fval�Ƿ��ص�����Ӧֵ
    result1(i,:) = [x,-fval];
end
result = sortrows(result1,nvars+1,'descend');     

[c, y] = SEIR_Country(C, result(1,:),pop);