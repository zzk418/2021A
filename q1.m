% ���Ʒ���ʾ��ͼ
d22 = xlsread('����1.csv'); % ԭ˳���Ӧ����
plot3(d22(:,1),d22(:,2),d22(:,3));

% clear,clc
% a0 = 0.00170; z0 = -300.8; % ��ʼ��
% a = 0.6*a0:0.0001:1.3*a0;
% n = 130;
% l = zeros(n,1); % �ٶ����ƶ�
% x = zeros(n,1);y = zeros(n,1);z = zeros(n,1);
% k = zeros(n,1);
% for i=1:length(a)% ����a
%     parfor j=1:n
%         l(j) = aa_f(a,z0,i,j); % ��aa����
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
% [vmax,imax] = max(k);
% aa = a(imax); % ���¼���l
% 
% 
% %% a�����l
% clear,clc
% aa = 0.00176; n = 130;z0 = -300.8; load data d1
% l = zeros(n,1);
% for i=1:length(aa)% ����a
%     parfor j=1:n
%         l(j) = aa_f(aa,z0,i,j); % ��aa����
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
%    [vmax,imax] = max(k);
%    
% %% ������ʼz0
% clear,clc
% z0 = -300.4;
% aa = 0.00176; n = 130;
% z00 = z0-0.5:0.1:z0+0.5;
% k = zeros(n,1);l = zeros(n,1);
% for i=1:length(z00)% ����z
%     parfor j=1:n
%         l(j)=z00_f(aa,z00,i,j);
%     end
%     k(i) = sum(l>-0.6 & l<0.6);
% end
% [vmax,imax] = max(k);
% zz0 = z00(imax);

% ����Լ������һԼ����q11
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
xlabel('�����浥�������ڵ����');
ylabel('�ٶ����ƶ�����');
title('�ٶ���Լ�����')

%% ���ƽ���ͼ
clear, clc
aa = 0.00181;zz0 = -300.8271;
s_plot1(aa,zz0)% ����ƽ��ͼ

%% �����ȷ���a
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
title('����a�����ȼ���')
xlabel('a����ֵ');ylabel('����Լ�������������ڵ���');
[s_vmax,s_imax] = max(k);

%% z0������
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
title('����z0�����ȼ���')
xlabel('z0����ֵ');ylabel('����Լ�������������ڵ���');
[s_vmax,s_imax] = max(k);

%% �Ӻ�������
function s_plot(aa,zz0)
% ���ƽ���ͼ
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

axis([-100 250 -250 250 -301 -160]);
title('��׼�����빤������')
end

function s_plot1(aa,zz0)
% ���ƽ���ͼ
load data d22
x = -200:5:200; 
y = aa*x.^2+zz0;
plot(x,y,'r');hold on
y = -sqrt(300^2-x.^2);
plot(x,y,'b');
axis([-250 250 -301 -160]);
title('��׼�����빤������')
end

function [l,d]=z00_f(aa,z00,i,j)
% ��ⷽ��
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
% ��aa����
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
% for i=1:length(a)% ����a
%     for j=1:length(z)% ����z
%         % ���������㶼��������������
%         % ��ⷴ�����
%         % �������Լ�������˳�
%     end
% end


% [x,y,z] = deal(d22(:,1),d22(:,2),d22(:,3)); % ԭ��
% [X,Y] = meshgrid(x(1:5:end),y(1:5:end));
% Z=griddata(x,y,z,X,Y);
% surf(X,Y,Z);
% shading interp 













