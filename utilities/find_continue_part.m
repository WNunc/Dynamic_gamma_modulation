function data= find_continue_part(a)
% 输入一维向量，找到其中连续的部分，按片段长度排序
% 输出cell格式
%   此处显示详细说明
inter=find(diff(a)~=1);
if isempty(inter)
    data{1} = a;
else
for i=1:length(inter)+1
    if i==1
       data{1}=a(1:inter(1));
    elseif i==length(inter)+1
           data{i}=a(inter(i-1)+1:end);
    else
        data{i}=a(inter(i-1)+1:inter(i));
    end
end
for ii = 1:length(data)
    l(ii) = length(data{ii});
end
[~,I] = sort(l,'descend');
data = data(I);
end

