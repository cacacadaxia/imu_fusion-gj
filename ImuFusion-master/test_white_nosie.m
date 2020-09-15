


std_ = 0.02;
noise = std_*randn(1,1e6);
ksdensity(noise)
ksdensity(cos(noise/2))
zeta = std_;
 [out1,out2 ]  = count_perset(sin(noise/2) ,zeta ,0);

list = pos(:,1:3)- x_r(1:3,:)';
figure;ksdensity(list(:,3));
 [out1,out2 ] = count_perset(list ,0.08 ,0)
function [out1,out2 ] = count_perset(l ,zeta ,miu)%%数列、标准差、均值
cunt = 0;
for i = 1:length(l)
    if l(i)>miu-zeta&&l(i)<miu+zeta
        cunt = cunt+1;
    end
end
out1 = cunt/length(l);
cunt = 0;
for i = 1:length(l)
    if l(i)>miu-3*zeta&&l(i)<miu+3*zeta
        cunt = cunt+1;
    end
end
out2 = cunt/length(l);
end