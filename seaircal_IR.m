%% �����ܣ�Ϊ�Ŵ��㷨�ṩ����Ӧ������SEAIIDRģ�ͣ�������simerr��Ϊ����ָ�꣬��ÿ�ν�����ѡ����С��simerrֵ
function simerr = seaircal_IR(a)
%% �������
global C;  %������
global t; %ʱ��
global I;global di; global R;global pop;
I(1) = pop(1); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
R(1) = pop(2); %������Ⱥ��
di(1) = pop(3);
%% ѭ���Զ�����1����⴫Ⱦ��ģ��
for i = 1:1:t-1
    lambda = a(1) + a(2)/(1+exp(a(3)*(a(4)-i))); 
    di(i+1) = lambda*I(i)- a(5)*I(i); 
    I(i+1) = I(i)+lambda*I(i) - a(5)*I(i); 
    R(i+1) = R(i)+a(5)*I(i);
end
y = di';
%% ����R2ֵ
%�����������ָ����ģ�Ͳ���������SARS ������������ʵ���̵���϶�ֵ
simerr = -FITNESS(C(:,4),y); %��Ϊga()����Ĭ�������Ŀ�꺯���Ľ������Сֵ������R2�Ľ��Ҫ���Ŷԣ����ԼӸ���