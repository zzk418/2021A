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

