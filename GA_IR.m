%% ����һ�׶�����Ϊ��Ԥ��ʩ��IR����
%% ��ʼ������
clc;clear;
%% �̶�����
global C;  %������
global t; %ʱ��
global num; %�Ŵ��㷨
global I;global di; global R;global pop;

C = xlsread('�¹�-��������.xlsx','0'); %ֻ��Ķ�������Դ����
pop = [1 0 1];
LB = [0 0 0 0 0]; 
UB = [1 1 5 100 0.1];
nvars = 5; %�����Ƶ�ģ�Ͳ�������
num = 50; %�����Ŵ��㷨�ظ��Ĵ���
t = size(C,1);
result1 = zeros(num,nvars+1);
I = zeros(1,t);di = zeros(1,t); R = zeros(1,t);

%��������õ������ӳ���seaircal()��seair_constraint()
ObjectiveFunction = @seaircal_IR;
ConstraintFunction = @seair_constraint;
    
%ѭ��ga()�Ŵ��㷨
for i =1:num
    [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[]); %x�Ǵ���������fval�Ƿ��ص�����Ӧֵ
    result1(i,:) = [x,-fval];
end
%������������Ϊ*.xls ��ʽ�ĵ��ӱ���ĵ�
result = sortrows(result1,nvars+1,'descend');  

[c, y] = IR_Country(C,result(1,:),pop); %���ص�y�Ǿ���ÿ��������Ⱦ������c���ǽ�ֹ�������������Ⱥ����һЩ���������ڴ��뵽��һ��������p