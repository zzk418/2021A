% load q2.mat
% load data2.mat
% 
% n_index = [];
% for i = 1:4300
%    for j = 1:3
%     if find(i_index==index(i,j))
%        n_index = [n_index;index(i,:)]; 
%        break;
%     end
%    end    
% end

clear 
load q3.mat % d222变换前 d22变换后

[x1,y1,r1,ddet1,degree1] = solve_q3(d222,n_index);
[x2,y2,r2,ddet2,degree2] = solve_q3(d22,n_index);

d1 = length(r1(r1<1))/length(r1);
d2 = length(r2(r2<1))/length(r2);

%% 绘图
alpha = 36.795; beta = 78.169; % (度)
aa = 0.00181; n = 678;
zz0 = -299.4367;
load f3 % 对应标号 1为抛物面，2为球面
load q3% 对应变换后坐标
xyz_w = zeros(length(nn_index1),3);
xyz_b = zeros(length(nn_index2),3);

for i=1:length(nn_index1)
    xyz_w(i,1) = d22(nn_index1(i),1);
    xyz_w(i,2) = d22(nn_index1(i),2);
    xyz_w(i,3) = d22(nn_index1(i),3);
end
figure,plot3(xyz_w(:,1),xyz_w(:,2),xyz_w(:,3),'bd');hold on; 

for i=1:length(nn_index2)
    xyz_b(i,1) = d22(nn_index2(i),1);
    xyz_b(i,2) = d22(nn_index2(i),2);
    xyz_b(i,3) = d22(nn_index2(i),3);
end
plot3(xyz_b(:,1),xyz_b(:,2),xyz_b(:,3),'rd');hold on
% legend('工作抛物面反射板','')
% 绘制基准球面
H = -sqrt(300^2 -(500/2)^2);
plot3(cosd(0:360)*500/2,sind(0:360)*500/2,...
    ones(1,361)*H,'color',[0.7 0.7 0.7])
load data d3 d21 d22
x_b=d22(:,1);y_b=d22(:,2);z_b=d22(:,3); % bais状态
[~,tri] = ismember(d3,d21(2:end,1));
plot3(x_b(tri(:,[1:3,1]))',y_b(tri(:,[1:3,1]))',z_b(tri(:,[1:3,1]))',...
    'color',[0.8,0.8,0.8]);hold on
plot3(0,0,-299.4367,'ro'); % 绘制原来中心
title('第三问反射区域俯视图');

xlabel('x(m)');ylabel('y(m)');zlabel('z(m)');
view([0,0,1]); axis equal

%% 子函数
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


function [x,y,r,delta,degree,nn_index] = solve_q3(d222,n_index)

p1 = zeros(1444,3);
p2 = zeros(1444,3);
v1 = zeros(1444,3);
v2 = zeros(1444,3);
z = zeros(1444,1);
nn_index = [];
    for i = 1:1444
    n1 = n_index(i,1);
    n2 = n_index(i,2);
    n3 = n_index(i,3);
    
    %p1为抛物面中心坐标，篇p2为球面中心坐标
    p1(i,:) = [d222(n1,1)+d222(n2,1)+d222(n3,1),d222(n1,2)+d222(n2,2)+d222(n3,2),...
        d222(n1,3)+d222(n2,3)+d222(n3,3)]/3;

%         
    a1 = (d222(n2,2)-d222(n1,2))*(d222(n3,3)-d222(n1,3))- ...
            (d222(n3,2)-d222(n1,2))*(d222(n2,3)-d222(n1,3));

%         
    b1 = (d222(n2,3)-d222(n1,3))*(d222(n3,1)-d222(n1,1))- ...
            (d222(n3,3)-d222(n1,3))*(d222(n2,1)-d222(n1,1));    

        
    c1 = (d222(n2,1)-d222(n1,1))*(d222(n3,2)-d222(n1,2))- ...
            (d222(n3,1)-d222(n1,1))*(d222(n2,2)-d222(n1,2));           

        
     %法向量
     v1(i,:) = [a1/sqrt(a1^2+b1^2+c1^2),b1/sqrt(a1^2+b1^2+c1^2),c1/sqrt(a1^2+b1^2+c1^2)];
%      %入射光与平面法向夹角
     theta1_1(i) = acosd(-c1/sqrt(a1^2+b1^2+c1^2));
     
     
     %反射光单位方向向量
     cos_theta1x_1 = -sign(p1(i,1))*a1/sqrt(a1^2+b1^2+c1^2);
     cos_theta1y_1 = -sign(p1(i,2))*b1/sqrt(a1^2+b1^2+c1^2);
     cos_theta1s_1 = cosd(theta1_1(i));
      v1(i,:) = [cos_theta1x_1 cos_theta1y_1 cos_theta1s_1];
   
        v_1 = [0 0 -1];
        v_2 = v_1 - 2*dot(v_1,v1(i,:))*v1(i,:);
        ddet(i) = det([v_1;v1(i,:);v_2]);
        degree(i,1) = acosd(dot(v_1,v1(i,:))/(norm(v_1)*norm(v1(i,:))));
         degree(i,2) = acosd(dot(v_2,v1(i,:)/(norm(v1(i,:))*norm(v_2))));
        
        t = (-160.2 - p1(i,3))/v_2(3);
        x(i) = v_2(1)*t + p1(i,1);
        y(i) = v_2(2)*t +p1(i,2);
        r(i) = sqrt(x(i)^2 + y(i)^2)/2;
        
        if r(i)<0.5
            z(i) = 1;
            nn_index =[nn_index n_index(i,:)];
        end
      
    end
        nsum = 0;
        n = 0;
        for i = 1:1440
            n = n + z(i)*abs(sin(theta1_1(i)));
            nsum = nsum + abs(sin(theta1_1(i)));
        end
        
        delta = n/nsum;
        nn_index = unique(nn_index);

end