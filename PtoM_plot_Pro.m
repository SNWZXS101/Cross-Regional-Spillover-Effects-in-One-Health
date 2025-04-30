clear all
clc
% global data 

data=xlsread('Data_all_sm','广东广西2');
% data(data==0)=NaN;
% data=data(1:66,:);
h = 0.0236 ;
Tend = 236 ;tt=1:h:Tend;
tt = linspace(1,Tend,10000);
h = tt(2)-tt(1);

nSteps = length(tt);

filepath="     D:\GMDX\硕士文件202311\斑块交叉感染逸出模型\Codes\PINN_fitted\Result1_try\Training-results-25-02-11_035516-set1     ";
cd(filepath);%打开文件夹
temp=dir('*.txt');%扫描文件夹中的所有txt文件，得到结构体temp变量
[m,~]=size(temp);%提取结构体temp中的行数，为下面定义cell数组做准备
dataSet=cell(m,1);%定义一个空的cell数组，用来存放每一个txt文件的数据
%通过for循环，将上面每一个txt中的数据存放到cell数组dataSet中
for i=1:m
    dataSet{i,1}=load(temp(i).name);
end

filepath="  D:\GMDX\硕士文件202311\斑块交叉感染逸出模型\Codes\PINN_fitted\Result1_try   ";
cd(filepath);%打开文件夹

I1_py = dataSet{1}; I2_py = dataSet{4};
Iq1_py = dataSet{7}; Iq2_py = dataSet{8};
R1_py = dataSet{9}; R2_py = dataSet{10};
S1_py = dataSet{11}; S2_py = dataSet{12};
Sq1_py = dataSet{13}; Sq2_py = dataSet{14};

I_py = [I1_py,I2_py] ;
Iq_py = [Iq1_py,Iq2_py] ;
R_py = [R1_py,R2_py] ;
S_py = [S1_py,S2_py] ;
Sq_py = [Sq1_py,Sq2_py] ;

beta = dataSet{15} ;
c1_o = dataSet{16} ;
c2_o = dataSet{17} ;
d1_o = dataSet{18} ; 
d2_o = dataSet{19} ;
gamma = dataSet{20} ;
m1_o = dataSet{21} ;
m2_o = dataSet{22} ;
mu = dataSet{23} ;
q1_o = dataSet{24} ;
q2_o = dataSet{25} ;

Sfc = 1 ;
t_data = 1:Tend;

if length(q1_o) == 1
    q1_o = q1_o*ones(length(t_data),1);
    q2_o = q2_o*ones(length(t_data),1);
end

c1_t = interp1(t_data, c1_o, tt, 'linear')';
c2_t = interp1(t_data, c2_o, tt, 'linear')';
q1_t = interp1(t_data, q1_o, tt, 'linear')';
q2_t = interp1(t_data, q2_o, tt, 'linear')';
m1_t = interp1(t_data, m1_o, tt, 'linear')';
m2_t = interp1(t_data, m2_o, tt, 'linear')';
d1_t = interp1(t_data, d1_o, tt, 'linear')';
d2_t = interp1(t_data, d2_o, tt, 'linear')';

% m = [m1 m2];
% d = [d1 d2];
beta = [beta beta] ;

% State your initial values of compartments
N10 = 12700e4 * Sfc * 1;
N20 = 5037e4 * Sfc;
I10 = I1_py(1) ;
I20 = I2_py(1) ;
I0 = [I10 I20] ;
S0 = [N10 - I10, N20 - I20] ;
Sq0 = [Sq1_py(1) Sq2_py(1)] ;
Iq0 = [Iq1_py(1) Iq2_py(1)] ;
R0 = [R1_py(1) R2_py(1)] ;
X0 = [S0 Sq0 I0 Iq0 R0];

% par1 = [ d1 d2 m1 m2  beta mu gamma];
% par2 = [ c1, c2, q1, q2];

% 初始化解矩阵
X = zeros(nSteps, length(X0));
X(1, :) = X0;

% Runge-Kutta 4阶迭代
for i = 1 : nSteps-1
    tn = tt(i);
    xn = X(i, :);
    c = [c1_t(i), c2_t(i)];% in for .. it's single value of the certain time point
    q = [q1_t(i), q2_t(i)];
    m = [m1_t(i), m2_t(i)];
    d = [d1_t(i), d2_t(i)];
    x = xn ;
    % 计算k1, k2, k3, k4
    
     % 状态变量分解
    S = x(1:2); Sq = x(3:4); I = x(5:6); Iq = x(7:8); R = x(9:10);
    N = [sum(x(1:2:end)), sum(x(2:2:end))];
    lambda = [c(1)*I(1)*S(1)/N(1), c(2)*I(2)*S(2)/N(2)];
    
    % 模型方程
    dS = -lambda + (1-q).*(1-beta).*(1-m).*lambda + mu.*Sq + (1-flip(q)).*(1-flip(beta)).*flip(m).*flip(lambda);
    dSq = q.*(1-beta).*lambda - mu*Sq;
    dI = (1-q).*beta.*(1-m).*lambda - gamma.*I + (1-flip(q)).*flip(beta).*flip(m).*flip(lambda) - d.*I + flip(d).*flip(I);
    dIq = q.*beta.*lambda - gamma.*Iq;
    dR = gamma .* (I + Iq);
    
    dt = h;
    % 四阶 Runge-Kutta 方法的 k1, k2, k3, k4 计算
    k1 = dt * [dS, dSq, dI, dIq, dR];
    k2 = dt * [dS, dSq, dI, dIq, dR];
    k3 = dt * [dS, dSq, dI, dIq, dR];
    k4 = dt * [dS, dSq, dI, dIq, dR];

    % 更新解
    X(i+1, :) = x + (k1 + 2*k2 + 2*k3 + k4) / 6;
%     X(i+1, :) = x + dt * [dS, dSq, dI, dIq, dR];
end

% options=odeset('Reltol',1e-8,'AbsTol',1e-8);
% [T,X]=ode45(@DWYodeX,tt,X0,options,par1,par2);

S = X(:,1:2); Sq = X(:,3:4); I = X(:,5:6); Iq = X(:,7:8); R = X(:,9:10);
N1 = sum(X(:,1:2:end),2);
N2 = sum(X(:,2:2:end),2);
N = [N1 N2];

c = [c1_t, c2_t];% To calculate the new_I whose steps is 10000
q = [q1_t, q2_t];
m = [m1_t, m2_t];
d = [d1_t, d2_t];

lambda = c.*I.*S./N ;
new_I = (1-q).*beta.*(1-m).*lambda  +  (1-flip(q,2)).*flip(beta).*flip(m,2).*flip(lambda,2) + flip(d,2).*flip(I,2) ;
new_Iq = q.*beta.*lambda ; 
new_m = (1-flip(q,2)).*flip(beta).*flip(m,2).*flip(lambda,2)  ;
new_d = flip(d,2).*flip(I,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zihao = 12;
Timew = 130;

figure(1)

new_I = new_I * 1;

TT = 1:Tend;
c = [c1_o, c2_o];% To plot the original c(t)
q = [q1_o, q2_o];
m = [m1_o, m2_o];
d = [d1_o, d2_o];

subplot(211)
yyaxis left
h01 = plot(data(1:end,1),'k.','Markersize',18) ; hold on
fill([tt,fliplr(tt)],[zeros(1,length(tt)) ,fliplr(new_I(:,1)')],    [0.25 0.41 1]    ,'Edgecolor',  [0.2 0.63 0.79]     ) ;  hold on
% h1 = plot(tt,new_I(:,1), '-'  ,'LineWidth',2,          'HandleVisibility','on'); hold on 
alpha(0.5)
ylabel('Daily infections')

yyaxis right
h01c = plot(TT,c(:,1) ,'LineWidth',3,          'HandleVisibility','on'); hold on 

set(gca,'linewidth',1,'fontsize',zihao);
% set(gca,'yscale','log')
% set(gca,'XLim',[0 300])
legend([h01, h01c],'Data','c1')
title('回代GD')

subplot(212)
yyaxis left
h02 = plot(data(1:end,2),'r.','Markersize',18) ; hold on
fill([tt,fliplr(tt)],[zeros(1,length(tt)) ,fliplr(new_I(:,2)')],     [1 0.27 0]   ,'Edgecolor',[0.78 0.38 0.08]       ) ; hold on
% h2 = plot(tt,new_I(:,2), '--' ,'LineWidth',2,    'HandleVisibility','on'); hold on
alpha(0.5)
set(gca,'linewidth',1,'fontsize',zihao);
% set(gca,'yscale','log')
% set(gca,'XLim',[0 300])
ylabel('Daily infections')

yyaxis right
h02c = plot(TT,c(:,2) ,'LineWidth',3,          'HandleVisibility','on'); hold on 

legend([h02, h02c],'Data','c2')
title('回代GX')
box on




figure(2)

new_I = new_I * 1;

% TT = 1:Tend;
% c = [c1_o, c2_o];
% q = [q1_o, q2_o];

subplot(211)
yyaxis left
h01 = plot(data(1:end,1),'k.','Markersize',18) ; hold on
fill([tt,fliplr(tt)],[zeros(1,length(tt)) ,fliplr(new_I(:,1)')],    [0.25 0.41 1]    ,'Edgecolor',  [0.2 0.63 0.79]       ) ;  hold on
% h1 = plot(tt,new_I(:,1), '-'  ,'LineWidth',2,          'HandleVisibility','on'); hold on 
alpha(0.5)
ylabel('Daily infections')

yyaxis right
h01c = plot(TT,q(:,1) ,'LineWidth',3,          'HandleVisibility','on'); hold on 

set(gca,'linewidth',1,'fontsize',zihao);
% set(gca,'yscale','log')
% set(gca,'XLim',[0 300])

legend([h01, h01c],'Data','q1')
title('回代GD')

subplot(212)
yyaxis left
h02 = plot(data(1:end,2),'r.','Markersize',18) ; hold on
fill([tt,fliplr(tt)],[zeros(1,length(tt)) ,fliplr(new_I(:,2)')],     [1 0.27 0]   ,'Edgecolor',[0.78 0.38 0.08]       ) ; hold on
% h2 = plot(tt,new_I(:,2), '--' ,'LineWidth',2,    'HandleVisibility','on'); hold on
alpha(0.5)
set(gca,'linewidth',1,'fontsize',zihao);
% set(gca,'yscale','log')
% set(gca,'XLim',[0 300])
ylabel('Daily infections')

yyaxis right
h02c = plot(TT,q(:,2) ,'LineWidth',3,          'HandleVisibility','on'); hold on 

legend([h02, h02c],'Data','q1')
title('回代GX')
box on


figure(3)
subplot(211)
h1 = plot(m(1:end,1),'k.','Markersize',18 ,   'HandleVisibility','on') ; hold on
h2 = plot(m(1:end,2),'bo','Markersize',5 ,   'HandleVisibility','on') ; hold on
% h1 = plot(TT,m(:,1) ,'LineWidth',3,          'HandleVisibility','on'); hold on 
% h2 = plot(TT,m(:,2) ,'LineWidth',3,          'HandleVisibility','on'); hold on 
legend([h1,h2],'m1','m2')
set(gca,'linewidth',1,'fontsize',zihao);

subplot(212)
h1 = plot(d(1:end,1),'k.','Markersize',18 ,   'HandleVisibility','on') ; hold on
h2 = plot(d(1:end,2),'bo','Markersize',5 ,   'HandleVisibility','on') ; hold on
% h1 = plot(TT,d(:,1) ,'LineWidth',3,          'HandleVisibility','on'); hold on 
% h2 = plot(TT,d(:,2) ,'LineWidth',3,          'HandleVisibility','on'); hold on 
legend([h1,h2],'d1','d2')
set(gca,'linewidth',1,'fontsize',zihao);



data = {m(1:end,1), m(1:end,2), d(1:end,1), d(1:end,2)};

groupNames = {'m_1', 'm_2', 'd_1', 'd_2'}; % 组标签
colors = [0.2 0.6 1; 0.6 0.4 1; 0.6 0.2 0.8; 1 0.6 0.4]; % RGB颜色
colors = [0, 191, 255 ; 34, 139, 34; 255, 64, 64; 255, 165, 0] /255;

% 创建小提琴图
figure;
subplot(121)
hold on;
vp = violin(data(:,1:2), groupNames);

% 美化小提琴图
for i = 1:length(vp)
    vp(i).FaceColor = colors(i, :); % 设置小提琴颜色
    vp(i).FaceAlpha = 0.4; % 设置透明度
    vp(i).EdgeColor = [0 0 0]; % 设置边界颜色

    % 绘制数据点散点图（透明圆点）
    jitterX = (rand(size(data{i})) - 0.5) * 0.2; % 轻微扰动X坐标
    scatter(i + jitterX, data{i}, 20, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.6);
end

% 设置坐标轴
set(gca, 'FontSize', 12);
ylabel('Drought Severity/Duration', 'FontSize', 14);
title('Violin Plot Example', 'FontSize', 16);
grid on;
box on;
hold off;

subplot(122)
hold on;
vp = violin(data(:,3:4), groupNames);

% 美化小提琴图
for i = 1:length(vp)
    vp(i).FaceColor = colors(i+2, :); % 设置小提琴颜色
    vp(i).FaceAlpha = 0.4; % 设置透明度
    vp(i).EdgeColor = [0 0 0]; % 设置边界颜色

    % 绘制数据点散点图（透明圆点）
    jitterX = (rand(size(data{i})) - 0.5) * 0.2; % 轻微扰动X坐标
    scatter(i + jitterX, data{i}, 20, colors(i+2, :), 'filled', 'MarkerFaceAlpha', 0.6);
end

% 设置坐标轴
set(gca, 'FontSize', 12);
ylabel('Drought Severity/Duration', 'FontSize', 14);
title('Violin Plot Example', 'FontSize', 16);
grid on;
box on;
hold off;

% 创建图形
figure;
subplot(211)
hold on;
% 画第一条曲线
plot(tt, new_m(:,1), '-', 'LineWidth', 2.5, 'Color', [70, 130, 180]/255); % 蓝色曲线
% 画第二条曲线
plot(tt, new_m(:,2), '--', 'LineWidth', 2.5, 'Color', [255, 99, 71]/255); % 橙色曲线
% 美化图形
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Daily spillover', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
title('Time Series Comparison', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
% 设置图例
legend({'Caused by m_1', 'Caused by m_2'}, 'Location', 'northwest', 'FontSize', 12, 'FontName', 'Arial', 'Box', 'off');
% 美化轴
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
         'XColor', [0.3, 0.3, 0.3], 'YColor', [0.3, 0.3, 0.3], ...
         'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');

% 添加网格（降低透明度）
grid on;
box on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7]; % 灰色网格
ax.GridAlpha = 0.3; % 透明度调整

subplot(212)
hold on;
% 画第一条曲线
plot(tt, new_d(:,1), '-', 'LineWidth', 2.5, 'Color', [70, 130, 180]/255); % 蓝色曲线
% 画第二条曲线
plot(tt, new_d(:,2), '--', 'LineWidth', 2.5, 'Color', [255, 99, 71]/255); % 橙色曲线
% 美化图形
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Daily spillover', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
title('Time Series Comparison', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');

% 设置图例
legend({'Caused by d_1', 'Caused by d_2'}, 'Location', 'northwest', 'FontSize', 12, 'FontName', 'Arial', 'Box', 'off');

% 美化轴
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
         'XColor', [0.3, 0.3, 0.3], 'YColor', [0.3, 0.3, 0.3], ...
         'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');

% 添加网格（降低透明度）
grid on;
box on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7]; % 灰色网格
ax.GridAlpha = 0.3; % 透明度调整


display_time_series_stats(d2_o);   %   c1_o q1_o        data{1}


% display_epidemic_stats_2(q2_o);





