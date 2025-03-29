function [ts_v,cycle_lvel,cycle_avel] = get_thetacycle_vel(ts,linvel,angvel,cyclets)
% input
% ts--timestamp that corresponds of vel 
% linvel, angvel--speed and angle speed
% cyclets--timestamp that cycle start and end
% output
% ts_v--这一个cycle中找到的速度的timestamp
% cycle_lvel--cycle中的平均速度
% cycle_avel--cycle中的平均角速度
%   此处显示详细说明
ind = find(ts>=cyclets(1) & ts<=cyclets(2));
ts_v = ts(ind);
cycle_lvel = mean(linvel(ind));
cycle_avel = mean(angvel(ind));
end

