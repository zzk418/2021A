% 相邻主索节点变化幅度
clear,clc
load data2.mat % 节点序号转换
z0 = -300.4;
% aa = 0.0017;
z00 =  -300.85; %z0-0.5:0.1:z0+0.5;
load data.mat d22
n = 2226;
l = zeros(n,1); % 促动器移动
k = zeros(n,1);
num = zeros(28,1);
for i = 1:28
    num(i) = [5*(1+i)*i/2+1];
end
num(29) = 2226;

for i=1:length(z00)
    aa = 1/(4*(-160.2 - z00(i)));
    x = d22(:,1);y = d22(:,2);z = d22(:,3);r = zeros(n,1);
    stop = 1;
    for j=1:n
        clear s
        syms t
        s = solve(d22(j,3).*t==aa*t^2.*(d22(j,1).^2+d22(j,2).^2)+z00(i),t);
        s = min(abs(double(s)));
        x(j) = d22(j,1).*double(s);
        y(j) = d22(j,2).*double(s);
        z(j) = d22(j,3).*double(s);
        l(j) = sqrt((d22(j,1)-x(j)).^2 +...
           (d22(j,2)-y(j)).^2+...
           (d22(j,3)-z(j)).^2);
        if z(j)<d22(j,3)
            l=-l; 
        end
        r(j) =sqrt( x(j)^2 + y(j)^2); 
        %连续找39个伸长量大于0.6或者半径大于150的（说明不再有点可以伸长至抛物线）
        if l(j)<-0.6 || l(j)>0.6||max(r)> 150
            stop = stop + 1;
            x(j) = d22(j,1);
            y(j) = d22(j,2);
            z(j) = d22(j,3);
        else
            stop = 1;
        end
        
        if stop == 39
            stopnum = j;
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
    delta = abs(newdistance - distance)./distance;
    [max_delta_index,max_delta]= max(delta); 
    
    if max(r) >150 & max_delta < 0.0007
        break;
    end
end
rate = sum(delta(:,1)>0.0007)/n;
%% 储存约束满足相对情况
z0 = [-300.80 -300.81 -300.82 -300.83 -300.84 -300.86];
rate = [0.2071 0.1752 0.0444 0.0876 0.0921 0.0988]; % z0的相对满足率
z0i = -300.80:-0.01:-300.86;
ratei = interp1(z0,rate,z0i);
figure,plot(z0i,ratei);
xlabel('z0值');
ylabel('相邻主索节点变化幅度不满足率');
title('相邻主索节点变化幅度约束情况')