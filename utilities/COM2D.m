function [objCOMx,objCOMy] = COM2D(Mat,flag)
%COM2D 此处显示有关此函数的摘要
%   此处显示详细说明
[N_x,N_y] = size(Mat);
if nargin == 1
    flag = 0;
end
sum_x = 0;
sum_y = 0;
area = 0;
switch flag
    case 1
        for i=1:N_x
            for j=1:N_y
                if fieldcomp(i,j)==iob
                    sum_x=sum_x+i * Mat(i,j);  %计算第Ｋ区域的横坐标总和
                    sum_y=sum_y+j * Mat(i,j);  %计算第Ｋ区域的纵坐标总和
                    area = area + Mat(i,j);
                end
            end
        end
    case 0
        for i=1:N_x
            for j=1:N_y
                sum_x=sum_x+i * Mat(i,j);  %计算第Ｋ区域的横坐标总和
                sum_y=sum_y+j * Mat(i,j);  %计算第Ｋ区域的纵坐标总和
                area = area + Mat(i,j);
            end
        end
end

objCOMx=sum_x/area;  %计算第Ｋ区域的质心横坐标
objCOMy=sum_y/area;%计算第Ｋ区域的质心纵坐标
end

