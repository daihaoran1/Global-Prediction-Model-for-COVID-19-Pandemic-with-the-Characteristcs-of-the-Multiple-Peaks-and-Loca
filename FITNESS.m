%% �������ܣ�������ʵ���ݺ͹��Ʋ���֮���R2ֵ��aΪ��ʵ���ݣ�bΪ����ֵ
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