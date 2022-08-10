% 绘制分区示意图
d22 = xlsread('附件1.csv'); % 原顺序对应坐标
plot3(d22(:,1),d22(:,2),d22(:,3));

% clear,clc
% a0 = 0.00170; z0 = -300.8; % 初始解
% a = 0.6*a0:0.0001:1.3*a0;
% n = 130;
% l = zeros(n,1); % 促动器移动
% x = zeros(n,1);y = zeros(n,1);z = zeros(n,1);
% k = zeros(n,1);
% for i=1:length(a)% 遍历a
%     parfor j=1:n
%         l(j) = aa_f(a,z0,i,j); % 求aa参数
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
% [vmax,imax] = max(k);
% aa = a(imax); % 重新检验l
% 
% 
% %% a检验得l
% clear,clc
% aa = 0.00176; n = 130;z0 = -300.8; load data d1
% l = zeros(n,1);
% for i=1:length(aa)% 遍历a
%     parfor j=1:n
%         l(j) = aa_f(aa,z0,i,j); % 求aa参数
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
%    [vmax,imax] = max(k);
%    
% %% 遍历初始z0
% clear,clc
% z0 = -300.4;
% aa = 0.00176; n = 130;
% z00 = z0-0.5:0.1:z0+0.5;
% k = zeros(n,1);l = zeros(n,1);
% for i=1:length(z00)% 遍历z
%     parfor j=1:n
%         l(j)=z00_f(aa,z00,i,j);
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
% [vmax,imax] = max(k);
% zz0 = z00(imax);

% 计算约束，另一约束在q11
clear,clc
aa = 0.00181; n = 130;
zz0 = -300.8271;
l = zeros(n,1);
d = zeros(n,1);
for i=1
    parfor j=1:n
        l(j)=z00_f(aa,zz0,i,j);
    end
    k = sum(l>-0.6 & l<0.6);
end
[vmax,imax] = max(k);
figure,plot(1:n,l);
xlabel('侧剖面单区主索节点序号');
ylabel('促动器移动距离');
title('促动器约束情况')

%% 绘制解结果图
clear, clc
aa = 0.00181;zz0 = -300.8271;
s_plot1(aa,zz0)% 绘制平面图

%% 灵敏度分析a
clear,clc
n = 130;zz0 = -300.8271;
a_s = 0.00172:0.00002:0.00184;
k = zeros(length(a_s),1);
for i=1:length(a_s)
    for j=1:n
        aa = a_s(i);
        l(j)=aa_f(aa,zz0,j);
    end
    k(i) = sum(l>-0.6 & l<0.6);
end
figure,plot(a_s(1:length(a_s)),k);
title('参数a灵敏度检验')
xlabel('a参数值');ylabel('满足约束条件的主索节点数');
[s_vmax,s_imax] = max(k);

%% z0灵敏度
clear,clc
n = 130; aa = 0.00181;
zz_s = -301.0:0.02:-300.8;
k = zeros(length(zz_s),1);
for i=1:length(zz_s)
    for j=1:n
        zz0 = zz_s(i);
        l(j)=aa_f(aa,zz0,j);
    end
    k(i) = sum(l>-0.6 & l<0.6);
end
figure,plot(zz_s(1:length(zz_s)),k);
title('参数z0灵敏度检验')
xlabel('z0参数值');ylabel('满足约束条件的主索节点数');
[s_vmax,s_imax] = max(k);

%% 子函数部分
function s_plot(aa,zz0)
% 绘制解结果图
load data d22
x = -250:5:250; y = -250:5:250;
[X1,Y1] = meshgrid(x,y);
Z1 = (aa*(X1.^2+Y1.^2)+zz0);
Z1(X1<0 & Z1>-300.9 | Z1>-240) = NaN;
% 只绘制下部
figure,mesh(X1,Y1,Z1);  hold on; % 网格为抛物面
[x,y,z] = deal(d22(:,1),d22(:,2),d22(:,3)); % 原序
[X2,Y2] = meshgrid(x(1:5:end),y(1:5:end));
Z2 = griddata(x,y,z,X2,Y2);
Z2(X2<0 & Z2<-100) = NaN;
surf(X2,Y2,Z2);  % 只绘制半边基准球面
shading interp;

axis([-100 250 -250 250 -301 -160]);
title('基准球面与工作抛面')
end

function s_plot1(aa,zz0)
% 绘制解结果图
load data d22
x = -200:5:200; 
y = aa*x.^2+zz0;
plot(x,y,'r');hold on
y = -sqrt(300^2-x.^2);
plot(x,y,'b');
axis([-250 250 -301 -160]);
title('基准球面与工作抛面')
end

function [l,d]=z00_f(aa,z00,i,j)
% 求解方程
load data d1
n = 200;
x = zeros(n,1);y = zeros(n,1);z = zeros(n,1);
clear s
syms t
s = solve(d1(j,3).*t==aa*t^2.*(d1(j,1).^2+d1(j,2).^2)+z00(i),t);
s = min(abs(double(s)));
x(j) = d1(j,1).*double(s);
y(j) = d1(j,2).*double(s);
z(j) = d1(j,3).*double(s);
l = sqrt((d1(j,1)-x(j)).^2 +...
   (d1(j,2)-y(j)).^2+...
   (d1(j,3)-z(j)).^2);
if z(j)<d1(j,3)
   l=-l; 
end

end


function [l]=aa_f(a,zz0,j)
% 求aa参数
clear s
load data d1
n = 200;x = zeros(n,1);y = zeros(n,1);z = zeros(n,1);
syms t
s = solve(d1(j,3).*t==a*t^2.*(d1(j,1).^2+d1(j,2).^2)+zz0,t);
s = min(abs(double(s)));
x(j) = d1(j,1).*double(s);
y(j) = d1(j,2).*double(s);
z(j) = d1(j,3).*double(s);
l = sqrt((d1(j,1)-x(j)).^2 +...
   (d1(j,2)-y(j)).^2+...
   (d1(j,3)-z(j)).^2);
if z(j)<d1(j,3)
   l=-l; 
end
end



%   Z = (a0*(X.^2+Y.^2)+z0);
% end


% a = 0.8*a0:0.0001:1.2*z0;
% z = z0-2:0.001:z0+2;
% for i=1:length(a)% 遍历a
%     for j=1:length(z)% 遍历z
%         % 假设主索点都在理想抛物面上
%         % 求解反射情况
%         % 如果满足约束条件退出
%     end
% end


% [x,y,z] = deal(d22(:,1),d22(:,2),d22(:,3)); % 原序
% [X,Y] = meshgrid(x(1:5:end),y(1:5:end));
% Z=griddata(x,y,z,X,Y);
% surf(X,Y,Z);
% shading interp 













