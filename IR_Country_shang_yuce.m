%% �����Ų���������㲢������ѵ���Ͻ��
no2 = 9; %no2��ֵ ���� 11 ����6  ӡ��15 Ӣ��9  ����˹14 �¹�10
pop = [7197293.644	43223444.56	27679.5603	0.611931266	0.032053263];
C = xlsread('Ӣ��-��������.xlsx','12'); %����


y_zhenshi = C(:,4);
y_zhenshi = y_zhenshi(~isnan(y_zhenshi));
result = xlsread('Ӣ��-��������.xlsx','Ԥ�����'); %����
t = size(C,1);y = zeros(t,size(result,1));r2= zeros(size(result,1),2);
for j = 1:size(result,1)
a = result(j,:);
 I = zeros(1,t);di = zeros(1,t);R = zeros(1,t);lambda = zeros(1,t);
I(1) = pop(1); R(1) = pop(2); di(1) = pop(3);
beta0 = a(1); %��Ⱦ��
Yi = pop(5); %������
beta1 = a(2);
%% ѭ���Զ�����1����⴫Ⱦ��ģ��
for i = 1:1:t-1
    lambda(i) = beta0 + beta1/(1+exp(a(3)*(a(4)-i))); 
    eta = a(5)+a(6)*(no2-C(i,1)); 
    di(i+1) = lambda(i)*I(i) - eta*I(i) - Yi*I(i);
    I(i+1) = I(i)+lambda(i)*I(i)-eta*I(i) - Yi*I(i);
    R(i+1) = R(i)+Yi*I(i);
end
y(:,j) = di';
r2(j,1) = j;
r2(j,2) = FITNESS(y_zhenshi,y(1:size(y_zhenshi,1),j));
end
r2 = sortrows(r2,2,'descend');
