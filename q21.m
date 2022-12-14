clear,clc;
load q2.mat
load data2.mat

x = d22(:,1);
y = d22(:,2);
z = d22(:,3);
r = zeros(2226,1);

z0 = -300.4;
% aa = 0.0017;
z00 =  -300.8271;
%z0-0.5:0.1:z0+0.5;
n = length(i_index);
l = zeros(n,1); % 促动器移动
k = zeros(n,1);
dmin = 1000;
for i = 1:28
num(i) = [5*(1+i)*i/2+1];
end
num(29) = 2226;

for i=1:length(z00)% 遍历z
    aa = 1/(4*(-160.2 - z00(i)));
    x = d22(:,1);y = d22(:,2);z = d22(:,3);r = zeros(n,1);
    stop = 1;
    for jj=1:n
        j = i_index(jj);
        clear s
        syms t
        s = solve(d22(j,3).*t==aa*t^2.*(d22(j,1).^2+d22(j,2).^2)+z00(i),t);
        s1 = min(abs(double(s)));
        x(j) = d22(j,1).*double(s1);
        y(j) = d22(j,2).*double(s1);
        z(j) = d22(j,3).*double(s1);
        l(j) = sqrt((d22(j,1)-x(j)).^2 +...
           (d22(j,2)-y(j)).^2+...
           (d22(j,3)-z(j)).^2);
        if z(j)<d22(j,3)
            l=-l; 
        end
        r(j) =sqrt( x(j)^2 + y(j)^2); 
        %连续找39个伸长量大于0.6或者半径大于150的（说明不再有点可以伸长至抛物线）
        if l(j)<-0.6 || l(j)>0.6
            x(j) = d22(j,1);
            y(j) = d22(j,2);
            z(j) = d22(j,3);
            l(j) = 0;
        else
        
        if r(j)> 150
            stop = stop + 1;
        else
            stop = 1;
        end
        
        end
        
        if stop >= 39
            stopnum = jj;
            for iii = 1:39
                ii = i_index(jj-iii+1);
            x(ii) = d22(ii,1);
            y(ii) = d22(ii,2);
            z(ii) = d22(ii,3);
            l(ii) = 0;
            end 
            break;
        end  
        
    end
    
    %计算最大变化幅度
    newdistance = zeros(4300,3);
    for ii = 1:4300
    distance1 = sqrt((x(index(ii,1))- x(index(ii,2)))^2+...
                    (y(index(ii,1))- y(index(ii,2)))^2+...
                     (z(index(ii,1))- z(index(ii,2)))^2);
    distance2 = sqrt((x(index(ii,1))- x(index(ii,3)))^2+...
                    (y(index(ii,1))- y(index(ii,3)))^2+...
                     (z(index(ii,1))- z(index(ii,3)))^2);
    distance3 = sqrt((x(index(ii,3))- x(index(ii,2)))^2+...
                    (y(index(ii,3))- y(index(ii,2)))^2+...
                     (z(index(ii,3))- z(index(ii,2)))^2);
                
    newdistance(ii,:) = [distance1,distance2,distance3]; 
    end
    distance_change = newdistance - distance;
    delta = abs(distance_change)./distance;
    [max_delta,max_delta_index]= max(delta); 
    delta(find(delta>0.0007))=1;

    mean1 = mean(r(527:601));
    mean2 = mean(r(602:682));
    dsum = sum(sum(delta==1));
    lsum = sum(abs(l));
    
    if dsum < dmin
        zt = i;
        dmin = dsum;
    end
       
end
    d222 = [x,y,z];
  alpha = 36.795; beta = 78.169; % (度)
% 输入alpha,beta,原坐标转化为新坐标
r1 = atand(sind(alpha)/tand(beta));
r2 = asind(cosd(alpha)*cosd(beta)); % 空间角度关系
p1 = [cosd(r1),sind(r1);-sind(r1),cosd(r1)]; % 变换矩阵
p2 = [cosd(r2),sind(r2);-sind(r2),cosd(r2)];
s1 = p2*[d222(:,1)';d222(:,3)']; % 第一次变换求解
s1 = [s1(1,:);d222(:,2)';s1(2,:)]; % x为原坐标
s2 = p1*[s1(2,:);s1(3,:)];
s2 = [s1(1,:)' s2(1,:)' s2(2,:)'];
