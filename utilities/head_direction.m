function [v_dir,U,V]= head_direction(x,y,length_num)
% ����˲ʱ�ٶȶ�Ӧ�ķ���
N = length(x);
v_dir = zeros(N,1);%����N��1�е�ȫ�����
U = zeros(N,1);%����N��1�е�ȫ�����
V = zeros(N,1);%����N��1�е�ȫ�����
num = (length_num+1):1:(N-length_num);
for ii = num
    x1 = x(ii+length_num)-x(ii-length_num);
    y1 = y(ii+length_num)-y(ii-length_num);
    U(ii) = x1;
    V(ii) = y1;
    if x1 == 0 %�������к����㲻������ֵ�����
        if y1 == 0
            if ii ==1
                v_dir(ii) = 0;
            else
                v_dir(ii) = v_dir(ii-1); %��λ��Ϊ0����ͷ����������һʱ�̵�ͷ��������ͬ
            end
        elseif y1 > 0
            v_dir(ii) = 90;
        else 
            v_dir(ii) = 270;
        end
    elseif x1 > 0  
        v_dir(ii) = rad2deg(atan(y1/x1));
        if v_dir(ii) < 0 
            v_dir(ii) = v_dir(ii) + 360;
        end
    else 
        v_dir(ii) = rad2deg(atan(y1/x1)) + 180;
    end
end
for i = 1:length_num
    v_dir(i) = v_dir(1+length_num);
    v_dir(end-i+1) = v_dir(end-length_num);
end