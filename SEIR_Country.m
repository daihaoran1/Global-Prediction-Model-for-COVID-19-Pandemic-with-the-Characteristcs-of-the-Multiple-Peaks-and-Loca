%% �����ܣ������Ƶ��ĸ����Ҳ��������Լ�����ģ������������ʼ״̬
function [cc, y] = SEIR_Country(C,a,pop) 
%% �������
global Nh; global di;
global did; global sigma; global v; global Ya; global Yi; global Yid;
t = size(C,1); 
%������Ⱥ���������ڴ�
S = zeros(1,t); E = zeros(1,t);A = zeros(1,t);I = zeros(1,t);R = zeros(1,t);ID = zeros(1,t); EdConcount = zeros(1,t); lambda = zeros(1,t);
E(1) = pop(1); %Ǳ������Ⱥ��
A(1) = pop(2); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
I(1) = pop(3); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
ID(1) = pop(4); %ͨ�������и�Ⱦ�Ե���Ⱥ����������ҽԺ�������Ƶģ���������֢״��Ⱦ��Ⱥ����֢״��Ⱦ��Ⱥ��
R(1) = pop(5); %������Ⱥ��
S(1) = Nh - pop(1) - pop(2) - pop(3)- pop(4) -pop(5); %�׸�Ⱦ��Ⱥ��
EdConcount(1) = pop(6); %ÿ������ȷ������
lambda = zeros(1,t);
theta = a(5); %��֢״��Ⱦ�ߵļ����
fa = a(6); %��֢״��Ⱦ�ߵļ����
%% ѭ���Զ�����1����⴫Ⱦ��ģ��
for i = 1:1:t-1
    lambda(i) = a(1) + a(2)/(1+exp(a(3)*(a(4)-i))); %�׸�Ⱦ��Ⱥ��ΪǱ����Ⱥ�ĸ���
    S(i+1) = S(i)-lambda(i)*S(i); 
    E(i+1) = E(i)+lambda(i)*S(i)-sigma*E(i);
    A(i+1) = A(i)+v*sigma*E(i)-(Ya)*A(i);
    I(i+1) = I(i)+(1-v)*sigma*E(i)-(Yi+di)*I(i);
    ID(i+1) = ID(i)+theta*A(i)+fa*I(i)-(Yid+did)*ID(i);
    R(i+1) = R(i)+Yid*ID(i)+Yi*I(i)+Ya*A(i);
    EdConcount(i+1) = theta*A(i)+fa*I(i);
end
y = [EdConcount' lambda'];
cc = [S(end) E(end) A(end) I(end) ID(end) R(end) EdConcount(end) theta fa];
end