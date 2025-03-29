function lg = movdirection(lap_pos)
% 输入路径 一维角度bin
% 输出正向和负向路径信息lg，1=正向，-1=负向
%%
diff_pos = diff(lap_pos);
posdif = [0, diff_pos];
lg = 0;
for T = 1:length(posdif)-1
    if posdif(T+1)~=lg(T)&&posdif(T+1)~=0
        lg(T+1) = posdif(T+1);
    else
        lg(T+1) = lg(T);
    end
end
lg = lg(1:length(posdif));
lg(lg>0)=1;
lg(lg<0)=-1;
% figure;plot(lg)
% hold on
% plot(lap_pos)
% hold off