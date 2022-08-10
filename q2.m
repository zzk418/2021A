clear,clc
% 输入alpha,beta,原坐标转化为新坐标
alpha = 36.795; beta = 78.169; % (度)

%% 计算促动器径向移动
% dxyz =  xyz_w - xyz_b;
% dls = sqrt(dxyz(:,1).^2+dxyz(:,2).^2+dxyz(:,3).^2);
% for i=1:541
%     if dxyz(i,3)<0
%         dls(i) = -dls(i);
%     end
% end
load q22 l
n = 1350;
figure,plot(1:1350,l); 
title('促动器径向约束');
xlabel('反射区域内主索节点序号');
ylabel('促动器径向移动距离');
% %% 绘制平面图
% aa = 0.00181; zz0 = -294.4367;
% s_plot1(aa,zz0);
    
%% 绘图
clear,clc
alpha = 36.795; beta = 78.169; % (度)
load data d5 % 主索节点坐标（原序);口径以内索点
aa = 0.00181; n = 541;
zz0 = -300.8271;
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
xlabel('x');ylabel('y');zlabel('z');
axis([-100 250 -250 250 -301 -160]);
title('基准球面与工作抛面')

l = zeros(n,1);d = zeros(n,1);
% 储存联立的焦点
xs = zeros(n,1);ys = zeros(n,1);zs = zeros(n,1);
for i=1
    parfor j=1:n
        [l(j),xs(j),ys(j),zs(j)]=z00_f1(aa,zz0,i,j);
    end
    k = sum(l>-0.6 & l<0.6);
end
% figure,plot(1:n,l);
% xlabel('侧剖面单区主索节点序号');ylabel('促动器移动距离');
% 绘制小口径，工作抛面变换
[xyz_w] = shift_f(alpha,beta,[xs ys zs]);xyz_w = xyz_w';
[Xs,Ys] = meshgrid(xyz_w(1:5:end,1),xyz_w(1:5:end,2));
Zs = griddata(xyz_w(1:5:end,1),xyz_w(1:5:end,2),xyz_w(1:5:end,3),Xs,Ys);
figure,surf(Xs,Ys,Zs);hold on; shading interp;
% 小口径对应基准球面变换
[xyz_b] = shift_f(alpha,beta,d5);xyz_b = xyz_b';

% 绘制基准球面
H = -sqrt(300^2 -(500/2)^2);
plot3(cosd(0:360)*500/2,sind(0:360)*500/2,...
    ones(1,361)*H,'color',[0.7 0.7 0.7])
load data d3 d21 d22
x_b=d22(:,1);y_b=d22(:,2);z_b=d22(:,3); % bais状态
[~,tri] = ismember(d3,d21(2:end,1));
plot3(x_b(tri(:,[1:3,1]))',y_b(tri(:,[1:3,1]))',z_b(tri(:,[1:3,1]))',...
    'color',[0.8,0.8,0.8]);
plot3(0,0,-300.4,'ro'); % 绘制原来中心
hold on
view([0,0,1]); axis equal
xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');

%% 输出结果
clear,clc
load q22.mat
load data2.mat
str = {};
num = [];
for i = 1:length(l)
    if abs(l(i))>0
        num  = [num;i];
    end
end
d = zeros(length(num),3);

for i = 1:length(num)
    str{i} = ch2index{num(i)};
    d(i,:) = s2(num(i),:);
    d(i,find(abs(d(i,:))<0.00000001)) = 0;
    ll(i) = l(num(i));
end
str = str';
ll = ll';

%% 函数部分
% function s_plot1(aa,zz0)
% % 绘制中间截面图
% load data d22
% x = -250:5:250; 
% z = aa*x.^2+zz0;
% [xs;zs] = 
% plot(x,z,'r');hold on
% y = -sqrt(300^2-x.^2);
% plot(x,y,'b');
% axis([-250 250 -301 -160]);
% title('基准球面与工作抛面')
% end


function [s2] = shift_f(alpha,beta,d22)
% 输入alpha,beta,原坐标转化为新坐标
r1 = atand(sind(alpha)/tand(beta));
r2 = asind(cosd(alpha)*cosd(beta)); % 空间角度关系
p1 = [cosd(r1),sind(r1);-sind(r1),cosd(r1)]; % 变换矩阵
p2 = [cosd(r2),sind(r2);-sind(r2),cosd(r2)];

s1 = p1\[d22(:,2)';d22(:,3)']; % 第一次变换求解
s1 = [d22(:,1)';s1]; % x为原坐标
s2 = p2\[s1(1,:);s1(3,:)];
s2 = [s2(1,:);s1(2,:);s2(2,:)];
% 最终变换坐标
% plot3(s2(1,:),s2(2,:),s2(3,:),'color',[0.8,0.8,0.8]);
end



function [l,x,y,z]=z00_f1(aa,z00,i,j)
% 在小口径求解方程,升级版
load data d5
clear s
syms t
s = solve(d5(j,3).*t==aa*t^2.*(d5(j,1).^2+d5(j,2).^2)+z00(i),t);
s = min(abs(double(s)));
x = d5(j,1).*double(s);
y = d5(j,2).*double(s);
z = d5(j,3).*double(s);
l = sqrt((d5(j,1)-x).^2 +...
   (d5(j,2)-y).^2+...
   (d5(j,3)-z).^2);
if z<d5(j,3)
   l=-l; 
end

end


function [l]=aa_f(a,z0,i,j)
% 求aa参数
clear s
load data d1
syms t
s = solve(d1(j,3).*t==a(i)*t^2.*(d1(j,1).^2+d1(j,2).^2)+z0,t);
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