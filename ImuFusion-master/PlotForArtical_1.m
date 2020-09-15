% =========================================================================
%
%                  论文画图
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 4月22日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%--------------------------------------------------------------------------
close all;
l1 = pos_measure(:,2)-pos(:,2);
l2 = pos(:,2)- x_r(2,:)';
figure;plot(l1(1:10:end),'--k');
hold on ;plot(l2(1:10:end),'r','LineWidth',1)

l1 = pos_measure(:,3)-pos(:,3);
l2 = pos(:,3)- x_r(3,:)';
figure;plot(l1(1:10:end),'--k');
hold on ;plot(l2(1:10:end),'r','LineWidth',1)

figure1 = figure;
subplot1 = subplot(3,1,1);
l1 = pos_measure(:,2)-pos(:,2);
l2 = pos(:,2)- x_r(2,:)';
plot(time(1:10:end),l1(1:10:end),'--k');
hold on ;plot(time(1:10:end),l2(1:10:end),'r','LineWidth',1)
set(gca,'Fontname','Times New Roman','fontsize',16);
xlabel('t /s')
ylabel('y-error /mm')
hl = legend('振动补偿算法','KF算法')

legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.13151658767773 0.603734579313592 0.252369668246445 0.0501193317422435],...
    'Orientation','horizontal');

% 创建 axes
axes1 = axes('Parent',figure1,...
    'Position',[0.697772511848348 0.863961813842482 0.214549763033171 0.0835322195704074]);
hold(axes1,'on');

% 创建 plot
plot(time(1:10:end),l2(1:10:end),'Parent',axes1,'Color',[1 0 0]);
xlim(axes1,[10 100]);
box(axes1,'on');
% 创建 arrow
% annotation(figure1,'arrow',[0.699052132701424 0.654028436018959],...
%     [0.86873508353222 0.768496420047733],'Color',[1 0 0],'LineWidth',2);

%%
subplot(3,1,2)
l1 = pos_measure(:,3)-pos(:,3);
l2 = pos(:,3)- x_r(3,:)';
plot(time(1:10:end),l1(1:10:end),'--k');
hold on ;plot(time(1:10:end),l2(1:10:end),'r','LineWidth',1)
set(gca,'Fontname','Times New Roman','fontsize',16);
xlabel('t /s')
ylabel('z-error /mm')


%%
subplot3 = subplot(3,1,3)
l1 = (x_oula'- pos(:,4:6))/pi*180;

plot(time(1:10:end) , l1(1:10:end,1) ,'k');hold on
plot(time(1:10:end) , l1(1:10:end,2) ,'r');hold on
plot(time(1:10:end) , l1(1:10:end,3) ,'g')
hl = legend('x axis','y axis','z axis')
xlabel('t /s')
ylabel('\phi-error /deg')

legend1 = legend(subplot3,'show');
set(legend1,...
    'Position',[0.13151658767773 0.603734579313592 0.252369668246445 0.0501193317422435],...
    'Orientation','horizontal');
set(gca,'Fontname','Times New Roman','fontsize',16);
