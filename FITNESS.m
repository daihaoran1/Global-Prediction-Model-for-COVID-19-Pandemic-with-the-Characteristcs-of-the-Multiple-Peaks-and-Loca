%% 函数功能：计算真实数据和估计参数之间的R2值，a为真实数据，b为估计值
function sr=FITNESS(a,b)
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