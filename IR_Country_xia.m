%% �����ܣ������Ƶ��ĸ����Ҳ��������Լ�����ģ������������ʼ״̬
function [cc, y] = IR_Country_xia(C,a,pop) 
t = size(C,1); 
%������Ⱥ���������ڴ�
global no2;
I = zeros(1,t);di = zeros(1,t);R = zeros(1,t);lambda = zeros(1,t);
I(1) = pop(1); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
R(1) = pop(2); %������Ⱥ��
di(1) = pop(3);
beta0 = a(1); %��Ⱦ��
Yi = pop(5); %������
beta1 = a(2);
for i = 1:1:t-1
    lambda(i) = beta0 + beta1/(1+exp(a(3)*(i-a(4))));
    eta = a(5)+a(6)*(no2-C(i,1)); %no2��ֵ ���� 11 ����  ӡ��14 Ӣ��8  ����˹
    di(i+1) = lambda(i)*I(i)-eta*I(i) - Yi*I(i);
    I(i+1) = I(i)+lambda(i)*I(i)-eta*I(i) - Yi*I(i);
    R(i+1) = R(i)+Yi*I(i);
end
y = [di' lambda'];
cc = [I(end)  R(end)  di(end) a(1) pop(5)];
end