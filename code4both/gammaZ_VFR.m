
vel_exp = Vel_win_part{nd}{nl};
for ia = 1:narea
%       if ~isempty(Sz0{ia,1}{nd}{nl})
        Sz0_lap = Sz_vel_part{ia,1}{nd}{nl}';
        
        % remove the nan
        ind = find(~isnan(mean(Sz0_lap)));
        Sz0_lap1 = Sz0_lap(:,ind);
        vel_lap1 = vel_exp(ind);
        
        % remove the outliers
        ind = find(mean(Sz0_lap1)<=zscoremax);
        Sz0_lap1 = Sz0_lap1(:,ind);
        vel_lap1 = vel_lap1(ind);

        VFR0_lap = [];
        
        for nb=1:nbin
            VFR0_lap(:,nb) = P_estimator(Sz0_lap1,vel_lap1,vr_center(nb),invh);
        end
        VFR_exp{ia}{nd,nl} = VFR0_lap;
        VFR{ia}{nd,ns}(:,:,nl)=VFR0_lap;
%      end
end