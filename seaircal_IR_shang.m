%% �����ܣ�Ϊ�Ŵ��㷨�ṩ����Ӧ������SEAIIDRģ�ͣ�������simerr��Ϊ����ָ�꣬��ÿ�ν�����ѡ����С��simerrֵ
function simerr = seaircal_IR_shang(a)
%% �������
global C;  %������
global t; %ʱ��
global I;global di; global R;global pop;global no2;

I(1) = pop(1); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
R(1) = pop(2); %������Ⱥ��
di(1) = pop(3);
beta0 = a(1); %��Ⱦ��
Yi = pop(5); %������
beta1 = a(2);
%% ѭ���Զ�����1����⴫Ⱦ��ģ��
for i = 1:1:t-1
    lambda = beta0 + beta1/(1+exp(a(3)*(a(4)-i))); 
    eta = a(5)+a(6)*(no2-C(i,1)); %no2��ֵ ���� 11 ����  ӡ��14 Ӣ��8  ����˹
    di(i+1) = lambda*I(i)-eta*I(i) - Yi*I(i);
    I(i+1) = I(i)+lambda*I(i)-eta*I(i) - Yi*I(i);
    R(i+1) = R(i)+Yi*I(i);
end
y = di';
%% ����R2ֵ
%�����������ָ����ģ�Ͳ���������SARS ������������ʵ���̵���϶�ֵ
simerr = -FITNESS(C(:,4),y); %��Ϊga()����Ĭ�������Ŀ�꺯���Ľ������Сֵ������R2�Ľ��Ҫ���Ŷԣ����ԼӸ���