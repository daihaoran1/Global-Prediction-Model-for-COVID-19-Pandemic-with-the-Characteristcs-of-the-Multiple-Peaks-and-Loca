%% �����ܣ�Ϊ�Ŵ��㷨�ṩ����Ӧ������SEAIIDRģ�ͣ�������simerr��Ϊ����ָ�꣬��ÿ�ν�����ѡ����С��simerrֵ
function simerr = seaircal_SEIR_shang(a)
%% �������
global C;global t;global S;global E;global A;global I;global ID;global R;global EdConcount;global di;global did;global sigma;global v;global Ya;global Yi;global Yid;
global pop;global no2;

E(1) = pop(2); %Ǳ������Ⱥ��
A(1) = pop(3); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
I(1) = pop(4); %��֢״��Ⱦ��Ⱥ����δ��⵽�ģ�
ID(1) = pop(5); %ͨ�������и�Ⱦ�Ե���Ⱥ����������ҽԺ�������Ƶģ���������֢״��Ⱦ��Ⱥ����֢״��Ⱦ��Ⱥ��
R(1) = pop(6); %������Ⱥ��
S(1) = pop(1)+pop(6); %�׸�Ⱦ��Ⱥ��
EdConcount(1) = pop(7); %ÿ������ȷ������
beta0 = a(1); %��Ⱦ��
beta1 = a(2);
theta = pop(8); %��֢״��Ⱦ�ߵļ����
fa = pop(9); %��֢״��Ⱦ�ߵļ����
%% ѭ���Զ�����1����⴫Ⱦ��ģ��
for i = 1:1:t-1
    lambda = beta0 + beta1/(1+exp(a(3)*(a(4)-i))); %�׸�Ⱦ��Ⱥ��ΪǱ����Ⱥ�ĸ���
    eta = a(5)+a(6)*(no2-C(i,1)); %a(5)+a(6)*(11-C(i,1))
    S(i+1) = S(i)-lambda*S(i)+eta*I(i)+eta*A(i); 
    E(i+1) = E(i)+lambda*S(i)-sigma*E(i);
    A(i+1) = A(i)+v*sigma*E(i)-(Ya)*A(i)-eta*A(i);
    I(i+1) = I(i)+(1-v)*sigma*E(i)-(Yi+di)*I(i)-eta*I(i);
    ID(i+1) = ID(i)+theta*A(i)+fa*I(i)-(Yid+did)*ID(i);
    R(i+1) = R(i)+Yid*ID(i)+Yi*I(i)+Ya*A(i);
    EdConcount(i+1) = theta*A(i)+fa*I(i);
end
y = EdConcount';
%% ����R2ֵ
%�����������ָ����ģ�Ͳ���������SARS ������������ʵ���̵���϶�ֵ
simerr = -FITNESS(C(:,4),y); %��Ϊga()����Ĭ�������Ŀ�꺯���Ľ������Сֵ������R2�Ľ��Ҫ���Ŷԣ����ԼӸ���