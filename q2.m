clear,clc
% ����alpha,beta,ԭ����ת��Ϊ������
alpha = 36.795; beta = 78.169; % (��)

%% ����ٶ��������ƶ�
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
title('�ٶ�������Լ��');
xlabel('���������������ڵ����');
ylabel('�ٶ��������ƶ�����');
% %% ����ƽ��ͼ
% aa = 0.00181; zz0 = -294.4367;
% s_plot1(aa,zz0);
    
%% ��ͼ
clear,clc
alpha = 36.795; beta = 78.169; % (��)
load data d5 % �����ڵ����꣨ԭ��);�ھ���������
aa = 0.00181; n = 541;
zz0 = -300.8271;
load data d22
x = -250:5:250; y = -250:5:250;
[X1,Y1] = meshgrid(x,y);
Z1 = (aa*(X1.^2+Y1.^2)+zz0);
Z1(X1<0 & Z1>-300.9 | Z1>-240) = NaN;
% ֻ�����²�
figure,mesh(X1,Y1,Z1);  hold on; % ����Ϊ������
[x,y,z] = deal(d22(:,1),d22(:,2),d22(:,3)); % ԭ��
[X2,Y2] = meshgrid(x(1:5:end),y(1:5:end));
Z2 = griddata(x,y,z,X2,Y2);
Z2(X2<0 & Z2<-100) = NaN;
surf(X2,Y2,Z2);  % ֻ���ư�߻�׼����
shading interp;
xlabel('x');ylabel('y');zlabel('z');
axis([-100 250 -250 250 -301 -160]);
title('��׼�����빤������')

l = zeros(n,1);d = zeros(n,1);
% ���������Ľ���
xs = zeros(n,1);ys = zeros(n,1);zs = zeros(n,1);
for i=1
    parfor j=1:n
        [l(j),xs(j),ys(j),zs(j)]=z00_f1(aa,zz0,i,j);
    end
    k = sum(l>-0.6 & l<0.6);
end
% figure,plot(1:n,l);
% xlabel('�����浥�������ڵ����');ylabel('�ٶ����ƶ�����');
% ����С�ھ�����������任
[xyz_w] = shift_f(alpha,beta,[xs ys zs]);xyz_w = xyz_w';
[Xs,Ys] = meshgrid(xyz_w(1:5:end,1),xyz_w(1:5:end,2));
Zs = griddata(xyz_w(1:5:end,1),xyz_w(1:5:end,2),xyz_w(1:5:end,3),Xs,Ys);
figure,surf(Xs,Ys,Zs);hold on; shading interp;
% С�ھ���Ӧ��׼����任
[xyz_b] = shift_f(alpha,beta,d5);xyz_b = xyz_b';

% ���ƻ�׼����
H = -sqrt(300^2 -(500/2)^2);
plot3(cosd(0:360)*500/2,sind(0:360)*500/2,...
    ones(1,361)*H,'color',[0.7 0.7 0.7])
load data d3 d21 d22
x_b=d22(:,1);y_b=d22(:,2);z_b=d22(:,3); % bais״̬
[~,tri] = ismember(d3,d21(2:end,1));
plot3(x_b(tri(:,[1:3,1]))',y_b(tri(:,[1:3,1]))',z_b(tri(:,[1:3,1]))',...
    'color',[0.8,0.8,0.8]);
plot3(0,0,-300.4,'ro'); % ����ԭ������
hold on
view([0,0,1]); axis equal
xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');

%% ������
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

%% ��������
% function s_plot1(aa,zz0)
% % �����м����ͼ
% load data d22
% x = -250:5:250; 
% z = aa*x.^2+zz0;
% [xs;zs] = 
% plot(x,z,'r');hold on
% y = -sqrt(300^2-x.^2);
% plot(x,y,'b');
% axis([-250 250 -301 -160]);
% title('��׼�����빤������')
% end


function [s2] = shift_f(alpha,beta,d22)
% ����alpha,beta,ԭ����ת��Ϊ������
r1 = atand(sind(alpha)/tand(beta));
r2 = asind(cosd(alpha)*cosd(beta)); % �ռ�Ƕȹ�ϵ
p1 = [cosd(r1),sind(r1);-sind(r1),cosd(r1)]; % �任����
p2 = [cosd(r2),sind(r2);-sind(r2),cosd(r2)];

s1 = p1\[d22(:,2)';d22(:,3)']; % ��һ�α任���
s1 = [d22(:,1)';s1]; % xΪԭ����
s2 = p2\[s1(1,:);s1(3,:)];
s2 = [s2(1,:);s1(2,:);s2(2,:)];
% ���ձ任����
% plot3(s2(1,:),s2(2,:),s2(3,:),'color',[0.8,0.8,0.8]);
end



function [l,x,y,z]=z00_f1(aa,z00,i,j)
% ��С�ھ���ⷽ��,������
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
% ��aa����
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