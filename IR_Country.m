%% �����ܣ������Ƶ��ĸ����Ҳ��������Լ�����ģ������������ʼ״̬
function [cc, y] = IR_Country(C,a,pop) 
t = size(C,1); 
%������Ⱥ���������ڴ�
I = zeros(1,t);di = zeros(1,t);R = zeros(1,t);lambda = zeros(1,t);

I(1) = pop(1); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
R(1) = pop(2); %������Ⱥ��
di(1) = pop(3);
for i = 1:1:t-1
    lambda(i) = a(1) + a(2)/(1+exp(a(3)*(a(4)-i))); 
    di(i+1) = lambda(i)*I(i)- a(5)*I(i); 
    I(i+1) = I(i)+lambda(i)*I(i) - a(5)*I(i); 
    R(i+1) = R(i)+a(5)*I(i);
end
y = [di' lambda'];
cc = [I(end)  R(end)  di(end) a(1)+a(2) a(5)];
end