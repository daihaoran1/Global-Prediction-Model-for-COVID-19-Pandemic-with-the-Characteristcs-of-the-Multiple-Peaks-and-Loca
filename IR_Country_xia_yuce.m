%% 将最优参数带入计算并返回最佳的拟合结果
no2 = 15; %no2均值 美国 11 巴西6  印度15 英国9  俄罗斯14  德国10
pop = [32509617.71	422484047.8	36183.04367	0.62157152	0.064329007];
C = xlsread('印度-疫情数据.xlsx','5'); %数据


y_zhenshi = C(:,4);
y_zhenshi = y_zhenshi(~isnan(y_zhenshi));
result = xlsread('印度-疫情数据.xlsx','预测参数'); %数据
t = size(C,1);y = zeros(t,size(result,1));r2= zeros(size(result,1),2);
for j = 1:size(result,1)
a = result(j,:);
 I = zeros(1,t);di = zeros(1,t);R = zeros(1,t);lambda = zeros(1,t);
I(1) = pop(1); R(1) = pop(2); di(1) = pop(3);
beta0 = a(1); %传染率
Yi = pop(5); %康复率
beta1 = a(2);
%% 循环以定步长1天求解传染病模型
for i = 1:1:t-1
    lambda(i) = beta0 + beta1/(1+exp(a(3)*(i-a(4)))); 
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
