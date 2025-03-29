 
%============================确定动物最终选择的位置编号，以及具体角度===========================%
ang_vel_limit = 5;
ind_rewardloc_offline=[];
ang_sample_test=[];
ts_restart=[];

for nl = 1:n_test
    %sample
    if ~isnan(ts_sample_reward(nl))
    ind=find(data_angle{1,2}{nl}(:,1)==ts_sample_reward(nl)/10^6);
    ang_sample_test(nl,1) = data_angle{1,2}{nl}(ind,2);
    a = abs(ang_sample_test(nl,1)-Ang_RewardLoc_ontrack(:,1));
    ind_rewardloc_offline(nl,1)=find(a==min(a));
    else
        ang_sample_test(nl,1) = nan;
        ind_rewardloc_offline(nl,1) = nan;
    end
    
    %test
    if ~isnan(ts_test_reward(nl))
    ind=find(data_angle{1,3}{nl}(:,1)==ts_test_reward(nl)/10^6);
    ang_sample_test(nl,2) = data_angle{1,3}{nl}(ind,2);
    a = abs(ang_sample_test(nl,2)-Ang_RewardLoc_ontrack(:,1));
    ind_rewardloc_offline(nl,2)=find(a==min(a));
    else
        ang_sample_test(nl,2) = nan;
        ind_rewardloc_offline(nl,2) = nan;
    end
end
%  save(trackdata_ns,'ang_sample_test','ind_rewardloc_offline','-append');
%ind_rewardloc_offline和sign_correct_test仍然对不上！！！！

%===========================Re-start time: when the rat restarts after getting rewards===============================%
            % Starting from reward location +1
            for nl = 1:n_test
            % Sample trials
            if ~isnan(sign_correct_sample(nl))
                ind_rewardloc = ind_rewardloc_offline(nl,1);
            else
                ind_rewardloc = Ind_rewardloc;
            end
            ang_reward_sample_border = Ang_RewardLoc_ontrack(ind_rewardloc+1,1);
            ind = find(data_angle{1,2}{nl}(:,1) > ts_sample_reward(nl)/10^6 &...
                data_angle{1,2}{nl}(:,2) >= ang_reward_sample_border &...
                data_angle{1,2}{nl}(:,3) >= ang_vel_limit);
            if ~isempty(ind)
            ind = ind(1);
            ts_restart(nl,1) = data_angle{1,2}{nl}(ind,1);%第一列sample
            else
                ts_restart(nl,1) = nan;
            end
            end
            
           for nl = 1:n_test
            % Test trials
            if ~isnan(sign_correct_test(nl))
                ind_rewardloc = ind_rewardloc_offline(nl,2);
            else
                ind_rewardloc = Ind_rewardloc;
            end
            ang_reward_test_border = Ang_RewardLoc_ontrack(ind_rewardloc-1,1);
            ind = find(data_angle{1,3}{nl}(:,1) > ts_test_reward(nl)/10^6 &...
                data_angle{1,3}{nl}(:,2) <= ang_reward_test_border &...
                data_angle{1,3}{nl}(:,3) >= ang_vel_limit);
             if ~isempty(ind)
            ind = ind(1);
            ts_restart(nl,2) = data_angle{1,3}{nl}(ind,1);%第二列test
            else
                ts_restart(nl,2) = nan;
            end
           end
          