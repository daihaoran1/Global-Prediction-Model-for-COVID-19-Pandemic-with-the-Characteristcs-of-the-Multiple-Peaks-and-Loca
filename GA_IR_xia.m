%% �������������½��ڵ�IRģ�Ͳ���
%% ��ʼ������
clc;clear;
%% �̶�����
global C;  %������
global t; %ʱ��
global num; %�Ŵ��㷨
global I;global di; global R;global pop;global no2;

C = xlsread('ӡ��-��������.xlsx','5'); %����
pop = [32509617.71	422484047.8	36183.04367	0.62157152	0.064329007];
no2 = 15; %no2��ֵ ���� 11 ����6  ӡ��15 Ӣ��9  ����˹14 �¹�10

LB = [0 0 0 0 -1 -1]; 
UB = [1 1 1 200 1 1];
% Ӣ�������Ʋ���Ϊ
% LB = [0 0 0 0 0 0]; 
% UB = [2 2 2 200 1 1];

nvars = 6; %�����Ƶ�ģ�Ͳ�������
num = 30; %�����Ŵ��㷨�ظ��Ĵ���
t = size(C,1);
%% ѭ��ÿ�����ҵĲ���
result1 = zeros(num,nvars+1);
%������Ⱥ���������ڴ�
I = zeros(1,t);di = zeros(1,t); R = zeros(1,t);
%��������õ������ӳ���seaircal()��seair_constraint()
ObjectiveFunction = @seaircal_IR_xia;
ConstraintFunction = @seair_constraint;

%ѭ��ga()�Ŵ��㷨
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x�Ǵ���������fval�Ƿ��ص�����Ӧֵ
    result1(i,:) = [x,-fval];
end

%������������Ϊ*.xls ��ʽ�ĵ��ӱ���ĵ�
result = sortrows(result1,nvars+1,'descend');  

[c, y] = IR_Country_xia(C,result(1,:),pop); %���ص�y�Ǿ���ÿ��������Ⱦ������c���ǽ�ֹ�������������Ⱥ����һЩ���������ڴ��뵽��һ��������pop