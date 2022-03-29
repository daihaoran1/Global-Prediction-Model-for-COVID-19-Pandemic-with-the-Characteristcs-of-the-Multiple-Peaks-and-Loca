%% 函数功能：计算真实数据和估计参数之间的R2值，a为真实数据，b为估计值,a为n行1列
% data = xlsread('印度-疫情数据.xlsx','Sheet4');
% [a,b]=R2(data(:,1),data(:,2));
function [sr, rmse]=R2(a,b)
sizeA=size(a);
top = 0;
for i=1:1:sizeA(1)
    top = top + (a(i,1)-b(i,1))^2;
end
bot=0;
ma=mean(a);
for i=1:1:sizeA(1)
    bot = bot + (a(i,1)-ma)^2;
end
sr = 1 - top/bot;

rm = 0;
for i = 1:sizeA(1)
    rm = rm + abs(a(i)-b(i))/b(i);
end
rmse = rm/sizeA(1);