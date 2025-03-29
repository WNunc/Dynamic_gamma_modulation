%% Figure look like:
%             Sample #                             Test #
% fast gamma power during sample       fast gamma power during test
% slow gamma power during sample       slow gamma power during test
%   running speed during sample          running speed during test
%   theta cycles during sample            theta cycles during test
% bayesian decoding during sample       bayesian decoding during test

% Edited on 06/18/2020
% Modified on plot_decoding_sampletest_v1
% plot the stop zone, but not the correct reward location

%%
disp(strcat('ntrial = ',num2str(nl)))
ang_reward = Ang_RewardLoc_ontrack(Ind_rewardloc);
if ~isnan(Ind_rewardloc_sample{1}(nl))
    ind0 = Ind_rewardloc_sample{1}(nl);
    ang_stop_sample = Ang_RewardLoc_ontrack(ind0-1:ind0);
    ang_stop_sample5 = Ang_RewardLoc_ontrack(ind0-5);
else
    ang_stop_sample = Ang_RewardLoc_ontrack(Ind_rewardloc-1:Ind_rewardloc);
    ang_stop_sample5 = Ang_RewardLoc_ontrack(Ind_rewardloc-5);
end
if ~isnan(Ind_rewardloc_test{1}(nl))
    ind0 = Ind_rewardloc_test{1}(nl);
    ang_stop_test = Ang_RewardLoc_ontrack(ind0-1:ind0);
    ang_stop_test5 = Ang_RewardLoc_ontrack(max(ind0-5,1));
else
    ang_stop_test = Ang_RewardLoc_ontrack(Ind_rewardloc-1:Ind_rewardloc);
    ang_stop_test5 = Ang_RewardLoc_ontrack(max(Ind_rewardloc-5,1));
end
%% Plot the gamma power over time before getting reward
t_lap_sample = scores_sample{nl,6};
ind_t_lap_sample = find(t_lap_sample >= ts_start_stop_sample_nl(1)...
    & t_lap_sample <= ts_start_stop_sample_nl(2));
t_lap_sample = t_lap_sample(ind_t_lap_sample);
t_lap0_sample = t_lap_sample - t_lap_sample(1);
t_lap_test = scores_test{nl,6};
ind_t_lap_test = find(t_lap_test >= ts_start_stop_test_nl(1)...
    & t_lap_test <= ts_start_stop_test_nl(2));
t_lap_test = t_lap_test(ind_t_lap_test);
t_lap0_test = t_lap_test - t_lap_test(1);

% Fast gamma power--sample
fgpow = scores_sample{nl,14};
fgpow_lap = nan(size(fgpow,1),length(t_lap_sample));
for nbin = 1:length(t_lap_sample)
    ind_bins_eeg = find(scores_sample{nl,10} >= t_lap_sample(nbin)-dt/2 & scores_sample{nl,10} <= t_lap_sample(nbin)+dt/2);
    fgpow_lap(:,nbin) = mean(fgpow(:,ind_bins_eeg),2);
end
fgpow_lap_sample_mean = mean(fgpow_lap,1);
% Fast gamma power--test
fgpow = scores_test{nl,14};
fgpow_lap = nan(size(fgpow,1),length(t_lap_test));
for nbin = 1:length(t_lap_test)
    ind_bins_eeg = find(scores_test{nl,10} >= t_lap_test(nbin)-dt/2 & scores_test{nl,10} <= t_lap_test(nbin)+dt/2);
    fgpow_lap(:,nbin) = mean(fgpow(:,ind_bins_eeg),2);
end
fgpow_lap_test_mean = mean(fgpow_lap,1);
fg_min = floor(min(min(fgpow_lap_sample_mean),min(fgpow_lap_test_mean)));
fg_max = ceil(max(max(fgpow_lap_sample_mean),max(fgpow_lap_test_mean)));

subplot(7,2,1)
% plot fg gamma power in sample lap
plot(t_lap0_sample,fgpow_lap_sample_mean, 'color', fg_color);
xlim([t_lap0_sample(1),t_lap0_sample(end)]);
ylim([fg_min,fg_max]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Z fg';'power'})
title(strcat('Sample ',num2str(nl)))
subplot(7,2,2)
% plot fg gamma power in test lap
plot(t_lap0_test,fgpow_lap_test_mean, 'color', fg_color);
xlim([t_lap0_test(1),t_lap0_test(end)]);
ylim([fg_min,fg_max]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Z fg';'power'})
if Sign_correct_test{1}(nl) == 1
    title(strcat('Test ',num2str(nl),' (correct)'))
elseif Sign_correct_test{1}(nl) == 0
    if Ind_rewardloc_test{1}(nl) > Ind_rewardloc
        ind_error = abs(Ind_rewardloc_test{1}(nl)-Ind_rewardloc);
        title(strcat('Test ',num2str(nl),' (error +',num2str(ind_error),')'))
    else
        ind_error = abs(Ind_rewardloc_test{1}(nl)-Ind_rewardloc);
        title(strcat('Test ',num2str(nl),' (error -',num2str(ind_error),')'));
    end
elseif isnan(Sign_correct_test{1}(nl))
    title(strcat('Test ',num2str(nl),' (non-stop)'))
end

% Slow gamma power
sgpow = scores_sample{nl,15};
sgpow_lap = nan(size(sgpow,1),length(t_lap_sample));
for nbin = 1:length(t_lap_sample)
    ind_bins_eeg = find(scores_sample{nl,10} >= t_lap_sample(nbin)-dt/2 & scores_sample{nl,10} <= t_lap_sample(nbin)+dt/2);
    sgpow_lap(:,nbin) = mean(sgpow(:,ind_bins_eeg),2);
end
sgpow_lap_sample_mean = mean(sgpow_lap,1);
sgpow = scores_test{nl,15};
sgpow_lap = nan(size(sgpow,1),length(t_lap_test));
for nbin = 1:length(t_lap_test)
    ind_bins_eeg = find(scores_test{nl,10} >= t_lap_test(nbin)-dt/2 & scores_test{nl,10} <= t_lap_test(nbin)+dt/2);
    sgpow_lap(:,nbin) = mean(sgpow(:,ind_bins_eeg),2);
end
sgpow_lap_test_mean = mean(sgpow_lap,1);
sg_min = floor(min(min(sgpow_lap_sample_mean),min(sgpow_lap_test_mean)));
sg_max = ceil(max(max(sgpow_lap_sample_mean),max(sgpow_lap_test_mean)));

subplot(7,2,3)
% plot sg gamma power in sample lap
plot(t_lap0_sample,sgpow_lap_sample_mean, 'color', sg_color);
xlim([t_lap0_sample(1),t_lap0_sample(end)]);
ylim([sg_min,sg_max]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Z sg';'power'})

subplot(7,2,4)
% plot sg gamma power in test lap
plot(t_lap0_test,sgpow_lap_test_mean, 'color', sg_color);
xlim([t_lap0_test(1),t_lap0_test(end)]);
ylim([sg_min,sg_max]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Z sg';'power'})


%% Running speed (rad/s)
vel_sample = scores_sample{nl,8}(4,ind_t_lap_sample);
vel_test = scores_test{nl,8}(4,ind_t_lap_test);
min_vel = floor(min(min(vel_sample),min(vel_test)));
max_vel = ceil(max(max(vel_sample),max(vel_test)));

% plot running speed in sample lap
subplot(7,2,5)
plot(t_lap0_sample,vel_sample, 'color', 'k');
xlim([t_lap0_sample(1),t_lap0_sample(end)]);
ylim([min_vel,max_vel]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Speed';'(rad/s)'})

% plot running speed in test lap
subplot(7,2,6)
plot(t_lap0_test,vel_test, 'color', 'k');
xlim([t_lap0_test(1),t_lap0_test(end)]);
ylim([min_vel,max_vel]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Speed';'(rad/s)'})

%% Plot theta cycles
% get theta EEG in sample lap
t_EEG_sample = scores_sample{nl,10};
ind_t_EEG_sample = find(t_EEG_sample >= ts_start_stop_sample_nl(1)...
    & t_EEG_sample <= ts_start_stop_sample_nl(2));
t_EEG_sample = t_EEG_sample(ind_t_EEG_sample);
t_EEG_sample = t_EEG_sample - t_lap_sample(1);
EEG_sample = scores_sample{nl,11}(:,ind_t_EEG_sample);
smo = 1000; %in index points
thetadelta = nan(size(EEG_sample,1),1);
for ntt = 1:size(EEG_sample,1)
    TFRt = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,6,10); %theta
    TFRd = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,2,4);%delta
    thetadelta(ntt,1) = mean(smooth(TFRt./TFRd,smo));
end
[~,ind_EEG] = max(thetadelta);
theta_EEG_sample = scores_sample{nl,18}(ind_EEG,ind_t_EEG_sample);

subplot(7,2,7)
theta_EEG_sample = theta_EEG_sample-mean(theta_EEG_sample);
plot(t_EEG_sample,theta_EEG_sample,'k')
xlim([t_lap0_sample(1),t_lap0_sample(end)]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Theta';'cycles'})

% get theta EEG in test lap
t_EEG_test = scores_test{nl,10};
ind_t_EEG_test = find(t_EEG_test >= ts_start_stop_test_nl(1)...
    & t_EEG_test <= ts_start_stop_test_nl(2));
t_EEG_test = t_EEG_test(ind_t_EEG_test);
t_EEG_test = t_EEG_test - t_lap_test(1);
EEG_test = scores_test{nl,11}(:,ind_t_EEG_test);
smo = 1000; %in index points
thetadelta = nan(size(EEG_test,1),1);
for ntt = 1:size(EEG_test,1)
    TFRt = TFR_frequency_band(EEG_test(ntt,:)',2000,5,6,10); %theta
    TFRd = TFR_frequency_band(EEG_test(ntt,:)',2000,5,2,4);%delta
    thetadelta(ntt,1) = mean(smooth(TFRt./TFRd,smo));
end
[~,ind_EEG] = max(thetadelta);
theta_EEG_test = scores_test{nl,18}(ind_EEG,ind_t_EEG_test);

subplot(7,2,8)
theta_EEG_test = theta_EEG_test-mean(theta_EEG_test);
plot(t_EEG_test,theta_EEG_test,'k')
xlim([t_lap0_test(1),t_lap0_test(end)]);
set(gca,'XTick',[]);
set(gca,'fontsize',24);
ylabel({'Theta';'cycles'})

%% Plot Bayesian decoding and theta cycles
subplot(7,2,9:2:13)
% plot decoding in sample lap
imagesc(t_lap0_sample,scores_sample{nl,5},scores_sample{nl,3}(:,ind_t_lap_sample))
%         y = colorbar;
%         y0 = ylabel(y,'Probability','FontSize',24);
%         set(y0, 'Units', 'Normalized', 'Position', [1.2, 0.5, 0]);
%         min0 = nanmin(nanmin(scores_sample{nsample,3}));
%         max0 = nanmax(nanmax(scores_sample{nsample,3}));
%         set(y,'YTick',[min0,max0]);
%         set(y,'YTickLabel',{'min','max'});
hold on
plot(t_lap0_sample,scores_sample{nl,5}(scores_sample{nl,4}(ind_t_lap_sample)),'w')
if ~isnan(Ind_rewardloc_sample{1}(nl))
    plot(t_lap0_sample,ones(size(t_lap0_sample))*ang_stop_sample(1),'--r')
    plot(t_lap0_sample,ones(size(t_lap0_sample))*ang_stop_sample(2),'--r')
else
    plot(t_lap0_sample,ones(size(t_lap0_sample))*ang_stop_sample(1),'--b')
    plot(t_lap0_sample,ones(size(t_lap0_sample))*ang_stop_sample(2),'--b')
end
% plot(t_lap0_sample,ones(size(t_lap0_sample))*ang_stop_sample5,'--g')
hold off
axis xy
ylim([0,2*pi])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','1/2 pi','pi','2/3 pi','2pi'});
set(gca,'fontsize',24);
xlabel('Time (s)')
ylabel('Location on track (rad)')


subplot(7,2,10:2:14)
% plot decoding in test lap
imagesc(t_lap0_test,scores_test{nl,5},scores_test{nl,3}(:,ind_t_lap_test))
%         y = colorbar;
%         y0 = ylabel(y,'Probability','FontSize',24);
%         set(y0, 'Units', 'Normalized', 'Position', [1.2, 0.5, 0]);
%         min0 = nanmin(nanmin(scores_test{ntest,3}));
%         max0 = nanmax(nanmax(scores_test{ntest,3}));
%         set(y,'YTick',[min0,max0]);
%         set(y,'YTickLabel',{'min','max'});
hold on
plot(t_lap0_test,scores_test{nl,5}(scores_test{nl,4}(ind_t_lap_test)),'w')
if ~isnan(Ind_rewardloc_test{1}(nl))
    plot(t_lap0_test,ones(size(t_lap0_test))*ang_stop_test(1),'--r')
    plot(t_lap0_test,ones(size(t_lap0_test))*ang_stop_test(2),'--r')
else
    plot(t_lap0_test,ones(size(t_lap0_test))*ang_stop_test(1),'--b')
    plot(t_lap0_test,ones(size(t_lap0_test))*ang_stop_test(2),'--b')    
end
% plot(t_lap0_test,ones(size(t_lap0_test))*ang_stop_test5,'--g')
hold off
axis xy
ylim([0,2*pi])
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',{'0','1/2 pi','pi','2/3 pi','2pi'});
set(gca,'fontsize',24);
xlabel('Time (s)')
ylabel('Location on track (rad)')