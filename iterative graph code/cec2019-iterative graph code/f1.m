clear all
clc

load('IIGKSOme_huatu.mat','diedaihui_IIGKSOme');
load('GKSO_huatu.mat', 'diedaihui_GKSO');
load('TTHHO_huatu.mat', 'diedaihui_TTHHO');
load('MTVSCA_huatu.mat', 'diedaihui_MTVSCA');
load('MGKSO_huatu.mat', 'diedaihui_MGKSO');
load('MDBO_huatu.mat', 'diedaihui_MDBO');
load('IGKSO_huatu.mat', 'diedaihui_IGKSO');
load('IGKSO2_huatu.mat', 'diedaihui_IGKSO2');
% load('ESO_huatu.mat', 'diedaihui_ESO');
load('ALA_huatu.mat', 'diedaihui_ALA');
% %------------diedaihui_IIGKSOme----------------------%
for i=1:10
    figure
fname=i;
% fname=9
d1=reshape(diedaihui_IIGKSOme(fname,1:100000),8,[]);
D1=min(d1,[],1);
yyyy1=log(D1);
plot(yyyy1,'Color',[0.04 0.09 0.27],'lineWidth',1.5)
% plot(yyyyy1,'r','lineWidth',1)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
hold on
% % % -----------------------diedaihui_GKSO-------------------
d2=reshape(diedaihui_GKSO(fname,1:100000),8,[]);
D2=min(d2,[],1);
yyyyy2=log(D2);
plot(yyyyy2,'Color',[0.1 0.8 0.9],'lineWidth',1.5)
% plot(yyyyy2,'g','lineWidth',1)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
 hold on



% % %---------------------diedaihui_TTHHO------------
d3=reshape(diedaihui_TTHHO(fname,1:100000),8,[]);
D3=min(d3,[],1);
yyyyy3=log(D3);
plot(yyyyy3,'Color',[0.42 0.35 0.80],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% % %-------------------diedaihui_MTVSCA-------------------------
d4=reshape(diedaihui_MTVSCA(fname,1:100000),8,[]);
D4=min(d4,[],1);
yyyyy4=log(D4);
plot(yyyyy4,'Color',[0.00 0.79 0.34],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% % %-------------------diedaihui_MGKSO-----------------------
% d5=reshape(every_bestf5(28,67:666),3,[]);
d5=reshape(diedaihui_MGKSO(fname,1:100000),8,[]);
D5=min(d5,[],1);
yyyyy5=log(D5);
% yyyyy5(1)=[8.74107688677474];
plot(yyyyy5,'Color',[1 0 0],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %------------------------diedaihui_MDBO-------------------
d6=reshape(diedaihui_MDBO(fname,1:100000),8,[]);
D6=min(d6,[],1);
yyyyy6=log(D6);
plot(yyyyy6,'Color',[1 0.5 1],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %------------------------diedaihui_IGKSO-------------------
d7=reshape(diedaihui_IGKSO(fname,1:100000),8,[]);
D7=min(d7,[],1);
yyyyy7=log(D7);
plot(yyyyy7,'Color',[0 0.2 0.6],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %------------------------diedaihui_IGKSO2------------------
d8=reshape(diedaihui_IGKSO2(fname,1:100000),8,[]);
D8=min(d8,[],1);
yyyyy8=log(D8);
plot(yyyyy8,'Color',[1 0.5 0],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %------------------------diedaihui_ESO-------------------
% d9=reshape(diedaihui_ESO(fname,1:100000),8,[]);
% D9=min(d9,[],1);
% yyyyy9=log(D9);
% plot(yyyyy9,'Color',[1 0.8 1],'lineWidth',1.5)
% set(gca,'XTick',0:40:1000)
% set(gca,'XLim',[0 200])
% set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
% %set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
% hold on
%------------------------diedaihui_diedaihui_ALA-------------------
d10=reshape(diedaihui_ALA(fname,1:100000),8,[]);
D10=min(d10,[],1);
yyyyy10=log(D10);
plot(yyyyy10,'Color',[0.6 0.2 0.8],'lineWidth',1.5)
set(gca,'XTick',0:2000:100000)
% %  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 10000])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% 


legend('IGKSO','GKSO','TTHHO','MTVSCA','MGKSO','MDBO','IGKSO1','EGKSO','ALA');
xlabel('FEs(\times10^5)');
ylabel('Fitness Value(log)');
% title('F1');
title(['F', num2str(i)]);
% title('F25');
end
%%%%F5
% 在原图上用半透明矩形标记放大区域
% zoom_x = [4000, 8000]; % 需要放大的x范围（根据实际调整）
% zoom_y = [0,0.5];
% rectangle('Position', [4000, 0, 4000, 0.5], ...
%           'EdgeColor', 'k', 'LineWidth', 1, 'LineStyle', '--', 'FaceColor', [0.9 0.9 0.9 0.2]);
% % 在主图中创建放大的子图
% axes('Position', [0.4 0.4 0.2 0.2]); % 位置和大小：[左下x, 左下y, 宽, 高]
% box on; % 显示边框
% 
% % 重新绘制放大区域的曲线
% hold on;
% plot(zoom_x(1):zoom_x(2), yyyy1(zoom_x(1):zoom_x(2)), 'Color', [0.04 0.09 0.27], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy2(zoom_x(1):zoom_x(2)), 'Color', [0.1 0.8 0.9], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy3(zoom_x(1):zoom_x(2)), 'Color', [0.42 0.35 0.80], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy4(zoom_x(1):zoom_x(2)), 'Color', [0.00 0.79 0.34], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy5(zoom_x(1):zoom_x(2)), 'Color', [1 0 0], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy6(zoom_x(1):zoom_x(2)), 'Color', [1 0.5 1], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy7(zoom_x(1):zoom_x(2)), 'Color', [0 0.2 0.6], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy8(zoom_x(1):zoom_x(2)), 'Color', [1 0.5 0], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy10(zoom_x(1):zoom_x(2)), 'Color', [0.6 0.2 0.8], 'LineWidth', 1.5);
% xlim(zoom_x);
% % annotation('arrow', [0.6 0], [0.5 0.5], 'Color', 'k', 'LineWidth', 1.5, 'HeadWidth', 15);
% set(gca,'XTicklabel',{'0','0.4','0.8'})
% ylim([0 0.5]); % 设置y轴范围为[0,0.5]
% annotation('arrow', ...
%            [0.6 0.5], [0.17 0.36], ... % 从主图到放大图
%            'Color', 'k', ...
%            'LineWidth', 0.5, ...
%            'HeadWidth', 15, ...
%            'HeadLength', 12);
% annotation('arrow', [0.3 0.4], [0.7 0.8], 'Color', 'k', 'LineWidth', 1.5, 'HeadWidth', 15);
% %%%F6
% % 在原图上用半透明矩形标记放大区域
% zoom_x = [4000, 8000]; % 需要放大的x范围（根据实际调整）
% zoom_y = [0,0.001];
% rectangle('Position', [4000, 0, 4000, 0.001], ...
%           'EdgeColor', 'k', 'LineWidth', 1, 'LineStyle', '--', 'FaceColor', [0.9 0.9 0.9 0.2]);
% % 在主图中创建放大的子图
% axes('Position', [0.4 0.5 0.2 0.2]); % 位置和大小：[左下x, 左下y, 宽, 高]
% box on; % 显示边框
% 
% % 重新绘制放大区域的曲线
% hold on;
% plot(zoom_x(1):zoom_x(2), yyyy1(zoom_x(1):zoom_x(2)), 'Color', [0.04 0.09 0.27], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy2(zoom_x(1):zoom_x(2)), 'Color', [0.1 0.8 0.9], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy3(zoom_x(1):zoom_x(2)), 'Color', [0.42 0.35 0.80], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy4(zoom_x(1):zoom_x(2)), 'Color', [0.00 0.79 0.34], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy5(zoom_x(1):zoom_x(2)), 'Color', [1 0 0], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy6(zoom_x(1):zoom_x(2)), 'Color', [1 0.5 1], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy7(zoom_x(1):zoom_x(2)), 'Color', [0 0.2 0.6], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy8(zoom_x(1):zoom_x(2)), 'Color', [1 0.5 0], 'LineWidth', 1.5);
% plot(zoom_x(1):zoom_x(2), yyyyy10(zoom_x(1):zoom_x(2)), 'Color', [0.6 0.2 0.8], 'LineWidth', 1.5);
% xlim(zoom_x);
% % annotation('arrow', [0.6 0], [0.5 0.5], 'Color', 'k', 'LineWidth', 1.5, 'HeadWidth', 15);
% set(gca,'XTicklabel',{'0','0.4','0.8'})
% ylim([0 0.001]); % 设置y轴范围为[0,0.5]
% annotation('arrow', ...
%            [0.6 0.5], [0.13 0.45], ... % 从主图到放大图
%            'Color', 'k', ...
%            'LineWidth', 0.5, ...
%            'HeadWidth', 15, ...
%            'HeadLength', 12);
% % annotation('arrow', [0.3 0.4], [0.7 0.8], 'Color', 'k', 'LineWidth', 1.5, 'HeadWidth', 15);
