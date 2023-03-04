%% 调整颜色、线型、线宽、标记、标记大小等,根据你的需求适当取舍
LineWidth = 0.8;
MarkerSize = 5;
MarkerSet = {'+','x','o','s','d','^','>','p','h','*'};
% o	圆圈
% +	加号
% *	星号
% .	点
% x	叉号
% s	方形
% d	菱形
% ^	上三角
% v	下三角
% >	右三角
% <	左三角
% p	五角形
% h 六角形
for i=1:10
    set(pic(i),...%'Color',[1 0 1],...%颜色 'LineStyle','-.',...%线型
        'LineWidth',LineWidth,...%线宽
        'Marker',MarkerSet{i},...%标记
        'MarkerSize',MarkerSize,...%标记大小 'MarkerEdgeColor','k',...%标记边缘颜色 'MarkerFaceColor',[.49 1 .63],...%标记颜色
        'DisplayName',['\lambda^{' num2str(i) '}'])%图例名称
end




%% 添加和调整图例
lgd = legend;
set(lgd,'Interpreter','Tex','Location','best');
set(lgd,'NumColumns',2 ,'FontSize',8)%,...图例列数 'LineWidth',1,...%图例线宽 'FontSize',4,...%设置字体大小 'FontName','Arial',...%图例字体 'EdgeColor',[0.850980392156863 0.325490196078431 0.0980392156862745],...%图例边缘颜色 'Color',[0.0588235294117647 1 1]);%图例颜色
%% 添加 x y 标签
xlabel('$1/\Delta x$','Interpreter','LaTex','FontSize',10);
ylabel('Error','Interpreter','LaTex','FontSize',10)
%
%% 坐标刻度设置，是否显示小坐标和对数坐标
set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

%% 添加标准
%legend off;
y_std = exp(-2*k*log(one_dx)-0.5);
std_pic = plot(one_dx(1:end),y_std(1:end),'HandleVisibility','off');
set(std_pic,...%'Color',[1 0 1],...%颜色 'LineStyle','-.',...%线型
    'LineWidth',3,...
    'LineStyle','--')%图例名称
%% 添加文字
%[x,y] = ginput(1)
%text(gca,16,0,'ddd')

%k2 ISO_k stabilized
annotation(figure(1),'textbox',...
    [0.2 0.64 0.0 0.0],...
    'String',['slope=' num2str(2*k)],...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'FontSize',14,...
    'EdgeColor',[1 1 1]);

%k2 NISO_k
% annotation(figure(1),'textbox',...
%     [0.2 0.64 0.0 0.0],...
%     'String',['slope=' num2str(2*k)],...
%     'FitBoxToText','off',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'EdgeColor',[1 1 1]);

%k2 ISO_k
% annotation(figure(1),'textbox',...
%     [0.2 0.67 0.0 0.0],...
%     'String',['slope=' num2str(2*k)],...
%     'FitBoxToText','off',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'EdgeColor',[1 1 1]);

% 
%k1 NISO_k 
% annotation(figure(1),'textbox',...
%     [0.2 0.63 0.0 0.0],...
%     'String',['slope=' num2str(2*k)],...
%     'FitBoxToText','off',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'EdgeColor',[1 1 1]);

%k1 ISO_k
% annotation(figure(1),'textbox',...
%     [0.2 0.62 0.0 0.0],...
%     'String',['slope=' num2str(2*k)],...
%     'FitBoxToText','off',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'EdgeColor',[1 1 1]);