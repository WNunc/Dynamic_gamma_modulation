% 计算头方向
% 统计
% (1)running时的头方向
% (2)sequence中的头方向

clear
close all
directories_allData_v0_allgood


length_num = 1;
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    
    path_ns = path{ns};
    disp(path_ns);
    cd(path_ns);
    load('Data_video.mat', 'data_video')
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Ang_StartZone_center')
    x = data_video(:,2);
    y = data_video(:,3);
    
    [x_rotated,y_rotated] = rotcoord(x,y,-Ang_StartZone_center);
    pos_xy = [x_rotated,y_rotated];
    
    [v_dir,U,V]= head_direction(x_rotated,y_rotated,length_num);
    u = cos(deg2rad(v_dir));%x 分量
    v = sin(deg2rad(v_dir));%y 分量
    m = [v_dir,U,V,u,v];
    figure
    
    quiver(x,y,u,v)
    hold on
    quiver(x,y,U,V,2)
    
    quiver(x_rotated,y_rotated,u,v)
    quiver(x_rotated,y_rotated,U,V,2)
    title([{'head_direction'},path(ns)])
    legend({'raw norm','raw','rotate norm','rotate',},'Location','northeastoutside')
    axis equal
    axis square
    xlim([-80,80])
    ylim([-80,80])
    saveas(gcf,['C:\Users\tju\OneDrive\gamma整理\投稿2\comment\HD-',num2str(ns)])
    save('Data_HeadDirection.mat','v_dir','pos_xy','U','V','length_num','Ang_StartZone_center')
end




function [x_rotated,y_rotated] = rotcoord(x,y,angle_rad)
% 将坐标点转换为复数数组，作为极坐标
z = complex(x, y);

% 将极坐标按逆时针旋转90度
z_rotated = abs(z) .* exp(1i * (angle(z) - angle_rad));

% 将旋转后的极坐标转换为实数数组
x_rotated = real(z_rotated);
y_rotated = imag(z_rotated);
end