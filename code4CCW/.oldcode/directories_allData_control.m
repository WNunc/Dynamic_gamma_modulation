 Dir = 'H:\neuralynx\';

isession = 0;
Ind_Rat = [];
path = {};
date = {};
trackdata = {};
CSClist_CA1 = {};
CSClist_PFC = {};
CSClist_PRL = {};

CSClist_HPC_R = {};
CSClist_PFC_R = {};
rewardlast = [];
rewardcurrent = [];

%% Rat 002 only hippocampus    
%  TT2、6、7、8、9、10、11、12、13、14 at CA1, TT5 at CA2
isession_Rat002 = 0;

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-27-CT\');
date_Rat002{isession_Rat002,1}='2020-03-27';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-27.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 8;%没有记录3-27之前的奖励点
rewardcurrent_Rat002(isession_Rat002,1) = 12;
% sign_correct_test = [0;0;1;0;0;0;0;NaN]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-30-CT\');
date_Rat002{isession_Rat002,1}='2020-03-30';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-30.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 9;
rewardcurrent_Rat002(isession_Rat002,1) = 7;
% sign_correct_test = [1;0;1;0;1;0;NaN;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-31-CT\');
date_Rat002{isession_Rat002,1}='2020-03-31';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-31.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 7;
rewardcurrent_Rat002(isession_Rat002,1) = 11;
% sign_correct_test = [NaN;1;0;1;1;0;1;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-01-CT\');
date_Rat002{isession_Rat002,1}='2020-04-01';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-01.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 11;
rewardcurrent_Rat002(isession_Rat002,1) = 8;
% sign_correct_test = [0;0;1;1;0;0;NaN;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-02-CT\');
date_Rat002{isession_Rat002,1}='2020-04-02';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-02.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 8;
rewardcurrent_Rat002(isession_Rat002,1) = 12;
% sign_correct_test = [0;1;0;0;1;0;1;0]

isession_Rat002 = isession_Rat002+1; 
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-05-CT\');
date_Rat002{isession_Rat002,1}='2020-04-05';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-05.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 12;
rewardcurrent_Rat002(isession_Rat002,1) = 9;
% sign_correct_test = [1;1;1;1;0;1;NaN;0]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-06-CT\');
date_Rat002{isession_Rat002,1}='2020-04-06';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-06.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 9;
rewardcurrent_Rat002(isession_Rat002,1) = 11;
% sign_correct_test = [1;1;0;0;1;0;0;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-08-CT\');
date_Rat002{isession_Rat002,1}='2020-04-08';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-08.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 11;
rewardcurrent_Rat002(isession_Rat002,1) = 7;
% sign_correct_test = [0;0;0;0;1;1;NaN;0]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-09-CT\');
date_Rat002{isession_Rat002,1}='2020-04-09';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-09.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 7;
rewardcurrent_Rat002(isession_Rat002,1) = 10;
% sign_correct_test = [0;1;1;1;0;1;0;0]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-13-CT\');
date_Rat002{isession_Rat002,1}='2020-04-13';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-13.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 10;
rewardcurrent_Rat002(isession_Rat002,1) = 12;
% sign_correct_test = [0;1;1;1;1;1;1;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-14-CT\');
date_Rat002{isession_Rat002,1}='2020-04-14';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-14.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 12;
rewardcurrent_Rat002(isession_Rat002,1) = 8;
% sign_correct_test = [1;0;0;0;1;1;1;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-05-05-CT\');
date_Rat002{isession_Rat002,1}='2020-05-05';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-05-05.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 8;
rewardcurrent_Rat002(isession_Rat002,1) = 10;
% sign_correct_test = [1;0;0;0;1;1;1;0]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-05-07-CT\');
date_Rat002{isession_Rat002,1}='2020-05-07';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-05-07.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA2{isession_Rat002,1}=[5];
csclist_Rat002_PFC{isession_Rat002,1}=[];
csclist_Rat002_PRL{isession_Rat002,1}=[];
csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
rewardlast_Rat002(isession_Rat002,1) = 10;
rewardcurrent_Rat002(isession_Rat002,1) = 7;
% sign_correct_test = [1;0;0;0;0;0;1;1]

isession = isession+isession_Rat002;
Ind_Rat002 = ones(isession_Rat002,1)*002;
Ind_Rat = [Ind_Rat;Ind_Rat002];
path = [path;path_Rat002];
trackdata = [trackdata;trackdata_Rat002];
CSClist_CA1 = [CSClist_CA1;csclist_Rat002_CA1];
CSClist_PFC = [CSClist_PFC;csclist_Rat002_PFC];
CSClist_PRL = [CSClist_PRL;csclist_Rat002_PRL];
CSClist_HPC_R = [CSClist_HPC_R;csclist_Rat002_R2];
rewardlast = [rewardlast;rewardlast_Rat002];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat002];
date = [date;date_Rat002];

%% Rat 003；右右
%  not easy to find trace, we make the default that all tetrodes have cells
%  are at CA1
isession_Rat003 = 0;

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-04-07-CT\');
date_Rat003{isession_Rat003,1}='2020-04-07';
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-04-07.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[2,4,5,6,7,8,9,10,12,13,14,15];
csclist_Rat003_PFC{isession_Rat003,1}=[];
csclist_Rat003_PRL{isession_Rat003,1}=[];
csclist_Rat003_R2{isession_Rat003,1}='HS3R1.ncs';
rewardlast_Rat003(isession_Rat003,1) = 12;%昨天的奖励点
rewardcurrent_Rat003(isession_Rat003,1) = 11;%今天的奖励点
% 62 CA1 place cells, 6 CA1 interneurons
% sign_correct_test = [0,NaN,1,1,1,1,1,1]

% isession_Rat003 = isession_Rat003+1;%test 8
% matlab报错，test8没有时间戳,没有存上Ts_sample和Ts_test,因此手动存了,手动存了还是报错，先舍弃再说
% path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-04-30-CT\');
% date_Rat003{isession_Rat003,1}='2020-04-30';
% trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-04-30_1.mat');%手动存了Ts_sample和Ts_test
% csclist_Rat003_CA1{isession_Rat003,1}=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17];
% csclist_Rat003_PRL{isession_Rat003,1}=[];
% csclist_Rat003_R2{isession_Rat003,1}='HS3R1.ncs';%海马的reference
% rewardlast_Rat003(isession_Rat003,1) = 13;
% rewardcurrent_Rat003(isession_Rat003,1) = 10;
% % 76 CA1 place cells, 8 CA1 interneurons
% % sign_correct_test = [1,0,0,0,1,0,0,NaN];

isession_Rat003 = isession_Rat003+1; 
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-05-03-CT\');
date_Rat003{isession_Rat003,1}='2020-05-03';
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-05-03.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[1,2,3,4,5,6,7,9,10,12,13,14,15,16];
csclist_Rat003_PFC{isession_Rat003,1}=[];
csclist_Rat003_PRL{isession_Rat003,1}=[];
csclist_Rat003_R2{isession_Rat003,1} ='HS3R1.ncs';
rewardlast_Rat003(isession_Rat003,1) = 10;
rewardcurrent_Rat003(isession_Rat003,1) = 12;
% 75 CA1 place cells, 6 CA1 interneurons
% % sign_correct_test = [0,1,0,0,0,NaN,NaN,0];

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-06-15-CT\');
date_Rat003{isession_Rat003,1}='2020-06-15';
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-06-15.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[2,4,6,9,10,12,13,14,17];
csclist_Rat003_PRL{isession_Rat003,1}=[];
csclist_Rat003_PFC{isession_Rat003,1}=[];
csclist_Rat003_R2{isession_Rat003,1} ='HS3R1.ncs';
rewardlast_Rat003(isession_Rat003,1) = 10;
rewardcurrent_Rat003(isession_Rat003,1) = 11;
% 38 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,0,1,1,0,NaN,0,0];

isession_Rat003 = isession_Rat003+1;% posttest3时结束实验
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-06-17-CT\');%某一天的数据
date_Rat003{isession_Rat003,1}='2020-06-17';
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-06-17.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[2,6,7,9,10,13,14,17];
csclist_Rat003_PFC{isession_Rat003,1}=[];
csclist_Rat003_PRL{isession_Rat003,1}=[];
csclist_Rat003_R2{isession_Rat003,1}='HS3R1.ncs';
rewardlast_Rat003(isession_Rat003,1) = 11;
rewardcurrent_Rat003(isession_Rat003,1) = 7;
% 38 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [NaN,0,NaN,0,1,0,1,1]

isession = isession+isession_Rat003;
Ind_Rat003 = ones(isession_Rat003,1)*3;
Ind_Rat = [Ind_Rat;Ind_Rat003];
path = [path;path_Rat003];
trackdata = [trackdata;trackdata_Rat003];
CSClist_CA1 = [CSClist_CA1;csclist_Rat003_CA1];
CSClist_PFC = [CSClist_PFC;csclist_Rat003_PFC];
CSClist_PRL = [CSClist_PRL;csclist_Rat003_PRL];
CSClist_HPC_R = [CSClist_HPC_R;csclist_Rat003_R2];
rewardlast = [rewardlast;rewardlast_Rat003];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat003];
date = [date;date_Rat003];

%% Rat004 后后，海马和前额叶植入，前额叶大血管破裂所以前额叶tetrode被血块挡住
%  he was sac at age 15-months, so it is not easy to determie whether he is
%  a normal rats or not, abandon
for i =1
isession_Rat004 = 0;

isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-14-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-14';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201014_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];%2,12
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 7; 
rewardcurrent_Rat004(isession_Rat004,1) = 12;
% sign_correct_test = [NaN;1;NaN;0;0;NaN;0;0]
% 21 CA1 place cells, 2 CA1 interneurons；
% 9 PFC neurons；

isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-15-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-15';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201015_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 12; % 正式纪录前一天，并没有记录
rewardcurrent_Rat004(isession_Rat004,1) = 9;
% sign_correct_test = [NaN;0;1;0;1;1;1;0]
% 29 CA1 place cells, 3 CA1 interneurons；
% 15 PFC neurons；
% 
isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-16-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-16';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201016_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 9; % 正式纪录前一天，并没有记录
rewardcurrent_Rat004(isession_Rat004,1) = 7;
% sign_correct_test = [1;0;1;0;1;0;1;0]
% 36 CA1 place cells, 2 CA1 interneurons；
% 15 PFC neurons；

isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-18-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-18';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201018_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 10; % 正式纪录前一天，并没有记录
rewardcurrent_Rat004(isession_Rat004,1) = 13;
% sign_correct_test = [1;0;0;0;NaN;0;0;0]
% 42 CA1 place cells, 1 CA1 interneurons；
% 17 PFC neurons；

isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-19-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-19';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201019_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 13; % 正式纪录前一天，并没有记录
rewardcurrent_Rat004(isession_Rat004,1) = 8;
% sign_correct_test = [0;0;1;1;0;NaN;NaN;0]
% 50 CA1 place cells, 1 CA1 interneurons；
% 16 PFC neurons；

isession_Rat004 = isession_Rat004+1;
path_Rat004{isession_Rat004,1} = strcat(Dir,'Rat004\2020-10-20-CT\');
date_Rat004{isession_Rat004,1} = '2020-10-20';
trackdata_Rat004{isession_Rat004,1} = strcat(path_Rat004{isession_Rat004,1},'20201020_CT_tracking.mat');
csclist_Rat004_CA1{isession_Rat004,1}=[12];
csclist_Rat004_CA2{isession_Rat004,1}=[13,16];
csclist_Rat004_CA3{isession_Rat004,1}=[14,15];
csclist_Rat004_PFC{isession_Rat004,1}=[8];%8,9,10
csclist_Rat004_R1_CA1{isession_Rat004,1} = 'HS1R1.ncs';
csclist_Rat004_R3_PFC{isession_Rat004,1} = 'HS3R1.ncs';
rewardlast_Rat004(isession_Rat004,1) = 8; % 正式纪录前一天，并没有记录
rewardcurrent_Rat004(isession_Rat004,1) = 12;
% sign_correct_test = [1;0;1;0;0;0;0;0]
% 47 CA1 place cells, 2 CA1 interneurons；
% 18 PFC neurons；

% isession_1 = isession_1+isession_Rat004;
% Ind_Rat004 = ones(isession_Rat004,1)*004;
% Ind_Rat_1 = [Ind_Rat_1;Ind_Rat004];
% path_1 = [path_1;path_Rat004];
% date_1 = [date_1;date_Rat004];
% trackdata_1 = [trackdata_1;trackdata_Rat004];
% CSClist_CA1_1 = [CSClist_CA1_1;csclist_Rat004_CA1];
% CSClist_PFC = [CSClist_PFC;csclist_Rat004_PFC];
end
%% Rat009 小胖，海马信号异常，只能用前额叶
%  TT14、15、16、17、18 at Prl,TT13 at Cg
for i = 1
isession_Rat009 = 0;

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-03-26_CT\');
date_Rat009{isession_Rat009,1} = '2021-03-26';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210326_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18]; 
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 8; 
rewardcurrent_Rat009(isession_Rat009,1) = 11;
% sign_correct_test = [1;0;1;0;0;1;1;1;1;0;0;1;NaN;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-03-27_CT\');
date_Rat009{isession_Rat009,1} = '2021-03-27';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210327_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18]; 
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 11; 
rewardcurrent_Rat009(isession_Rat009,1) = 9;
% sign_correct_test = [0;0;0;0;0;0;0;0;1;0;0;0;0;0;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-03-30_CT\');
date_Rat009{isession_Rat009,1} = '2021-03-30';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210330_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18]; 
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 9; 
rewardcurrent_Rat009(isession_Rat009,1) = 13;
% sign_correct_test = [1;0;0;1;1;1;1;0;0;1;0;0;1;1;0;1;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-01_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-01';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210401_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18]; 
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 7; 
rewardcurrent_Rat009(isession_Rat009,1) = 12;
% sign_correct_test = [0;0;0;0;1;0;1;0;0;0;1;1;1;1;0;1;0;1;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-03_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-03';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210403_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18]; 
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 12; 
rewardcurrent_Rat009(isession_Rat009,1) = 10;
% sign_correct_test = [0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-05_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-05';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210405_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 10; 
rewardcurrent_Rat009(isession_Rat009,1) = 7;
% sign_correct_test = [1;0;1;1;0;1;1;0;1;1;1;1;1;0;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-06_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-06';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210406_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 7; 
rewardcurrent_Rat009(isession_Rat009,1) = 13;
% sign_correct_test = [1;1;0;1;1;0;0;0;0;1;0;0;1;0;1;0;0;1;0;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-07_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-07';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210407_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 13; 
rewardcurrent_Rat009(isession_Rat009,1) = 8;
% sign_correct_test = [1;0;1;1;0;0;0;1;0;1;1;0;0;1;0;0;1;1;1;0]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-08_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-08';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210408_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[13,14,15,16,17];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 8; 
rewardcurrent_Rat009(isession_Rat009,1) = 12;
% sign_correct_test = [1;0;1;0;0;0;1;0;0;0;1;0;0;0;0;1;1;1;0;0]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-10_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-10';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210410_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 12; 
rewardcurrent_Rat009(isession_Rat009,1) = 8;
% sign_correct_test = [1;0;0;1;1;1;1;1;0;0;1;1;0;0;0;1;1;1;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-13_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-13';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210413_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[13,14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 13; 
rewardcurrent_Rat009(isession_Rat009,1) = 9;
% sign_correct_test = [NaN;0;0;0;0;0;0;0;1;1;0;0;1;0;1;1;0;0;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons；

% isession_Rat009 = isession_Rat009+1;
% path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-14_CT\');
% date_Rat009{isession_Rat009,1} = '2021-04-14';
% trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210414_CT_tracking.mat');
% csclist_Rat009_CA1{isession_Rat009,1}=[];
% csclist_Rat009_PRL{isession_Rat009,1}=[13,14,15,16,17,18];
% csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
% rewardlast_Rat009(isession_Rat009,1) = 9; 
% rewardcurrent_Rat009(isession_Rat009,1) = 12;
% % sign_correct_test = [1;0;1;0;0;0;0;1;1;0;0;1;0;0;1;0;0;0;0;0]
% % CA1 place cells,CA1 interneurons；
% % PFC neurons；

isession_Rat009 = isession_Rat009+1;
path_Rat009{isession_Rat009,1} = strcat(Dir,'Rat009\2021-04-19_CT\');
date_Rat009{isession_Rat009,1} = '2021-04-19';
trackdata_Rat009{isession_Rat009,1} = strcat(path_Rat009{isession_Rat009,1},'20210419_CT_tracking.mat');
csclist_Rat009_CA1{isession_Rat009,1}=[];
csclist_Rat009_PRL{isession_Rat009,1}=[13,14,15,16,17,18];
csclist_Rat009_R1{isession_Rat009,1} = 'HS3R1.ncs';
rewardlast_Rat009(isession_Rat009,1) = 11; 
rewardcurrent_Rat009(isession_Rat009,1) = 9;
% sign_correct_test = [0;0;1;1;0;1;1;1;1;1;1;1;0;1;1;0;1;1;1;1]
% CA1 place cells,CA1 interneurons；
% PFC neurons

% isession = isession+isession_Rat009;
% Ind_Rat009 = ones(isession_Rat009,1)*009;
% Ind_Rat = [Ind_Rat;Ind_Rat009];
% path = [path;path_Rat009];
% date = [date;date_Rat009];
% trackdata = [trackdata;trackdata_Rat009];
% CSClist_CA1 = [CSClist_CA1;csclist_Rat009_CA1];
% CSClist_PRL = [CSClist_PRL;csclist_Rat009_PRL];
% CSClist_PFC_R = [CSClist_PFC_R;csclist_Rat009_R1];
end
%% Rat 010 前额叶额bundle没有target上前额叶脑区，只有海马可以用
%  TT10、11、12、13、14、15 at CA1, TT16、17、18 at CA2
isession_Rat010 = 0;

isession_Rat010 = isession_Rat010+1; 
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-06-30_CT\');
date_Rat010{isession_Rat010,1} = '2021-06-30';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210630_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; 
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-01_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-01';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210701_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-04_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-04';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210704_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-05_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-05';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210705_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-06_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-06';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210706_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-12_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-12';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210712_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-12_CT_2\');
date_Rat010{isession_Rat010,1} = '2021-07-12-2';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210712_CT_2_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-13_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-13';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210713_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-13_CT_2\');
date_Rat010{isession_Rat010,1} = '2021-07-13-2';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210713_CT_2_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession_Rat010 = isession_Rat010+1; %Day
path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-15_CT\');
date_Rat010{isession_Rat010,1} = '2021-07-15';
trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210715_CT_tracking.mat');
csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
csclist_Rat010_CA3{isession_Rat010,1}=[];
csclist_Rat010_PFC{isession_Rat010,1}=[];
csclist_Rat010_PRL{isession_Rat010,1}=[];
csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

isession = isession+isession_Rat010;
Ind_Rat010 = ones(isession_Rat010,1)*010;
Ind_Rat = [Ind_Rat;Ind_Rat010];
path = [path;path_Rat010];
date = [date;date_Rat010];
trackdata = [trackdata;trackdata_Rat010];
CSClist_CA1 = [CSClist_CA1;csclist_Rat010_CA1];
CSClist_PFC = [CSClist_PFC;csclist_Rat010_PFC];
CSClist_PRL = [CSClist_PRL;csclist_Rat010_PRL];
CSClist_HPC_R = [CSClist_HPC_R;csclist_Rat010_R2];

%% Rat011  海马和前额叶已经确认可以用
%  all hippocampus tetrode at CA1, all prefrontal except TT9 at Prl and Cg2
%  use tetrodes that have cells for every days' analyse
isession_Rat011 = 0;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-16_CT_6\');
date_Rat011{isession_Rat011,1}='2021-07-16';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-16-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,6,7,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=11;
rewardcurrent_Rat011(isession_Rat011,1)=8;

isession_Rat011 = isession_Rat011+1;% 轨道左侧有50Hz
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-20_CT_9\');
date_Rat011{isession_Rat011,1}='2021-07-20';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-20-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,6,7,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=9;
rewardcurrent_Rat011(isession_Rat011,1)=12;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-21_CT_11\');
date_Rat011{isession_Rat011,1}='2021-07-21';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-21-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,6,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=12;
rewardcurrent_Rat011(isession_Rat011,1)=7;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-22_CT_11\');
date_Rat011{isession_Rat011,1}='2021-07-22';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-22-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,6];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1) = 11;
rewardcurrent_Rat011(isession_Rat011,1) = 8;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-23_CT_11\');
date_Rat011{isession_Rat011,1}='2021-07-23';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-23-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=11;
rewardcurrent_Rat011(isession_Rat011,1)=8;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-24_CT_8\');
date_Rat011{isession_Rat011,1}='2021-07-24';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-24-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,13,15,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=8;
rewardcurrent_Rat011(isession_Rat011,1)=10;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-25_CT_9\');
date_Rat011{isession_Rat011,1}='2021-07-25';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-25-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=10;
rewardcurrent_Rat011(isession_Rat011,1)=13;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-26_CT_8\');
date_Rat011{isession_Rat011,1}='2021-07-26';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-26-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=13;
rewardcurrent_Rat011(isession_Rat011,1)=9;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-27_CT_6\');
date_Rat011{isession_Rat011,1}='2021-07-27';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-27-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=9;
rewardcurrent_Rat011(isession_Rat011,1)=12;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-28_CT_PM_7\');
date_Rat011{isession_Rat011,1}='2021-07-28-PM';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-28-CT-PM.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,6];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=7;
rewardcurrent_Rat011(isession_Rat011,1)=10;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-29_CT_AM_15\');
date_Rat011{isession_Rat011,1}='2021-07-29-AM';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-29-CT-AM.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[2,3,5,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=10;
rewardcurrent_Rat011(isession_Rat011,1) =13;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-29_CT_PM_16\');
date_Rat011{isession_Rat011,1}='2021-07-29-PM';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-29-CT-PM.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[3,5,7,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=13;
rewardcurrent_Rat011(isession_Rat011,1)=9;

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-30_CT_PM_15\');
date_Rat011{isession_Rat011,1}='2021-07-30-PM';
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-30-CT-PM.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat011_PFC{isession_Rat011,1}=[3,5,6,7,8];
csclist_Rat011_PRL{isession_Rat011,1}=[1,2,3,4,5,6,7,8];
csclist_Rat011_R2{isession_Rat011,1}='HS3R1.ncs';
csclist_Rat011_R1{isession_Rat011,1}='HS1R1.ncs';
rewardlast_Rat011(isession_Rat011,1)=11;
rewardcurrent_Rat011(isession_Rat011,1)=8;

Ind_Rat011 = ones(isession_Rat011,1)*011;

isession = isession+isession_Rat011;
Ind_Rat = [Ind_Rat;Ind_Rat011];
path = [path;path_Rat011];
date = [date;date_Rat011];
trackdata = [trackdata;trackdata_Rat011];
CSClist_CA1 = [CSClist_CA1;csclist_Rat011_CA1];
CSClist_PFC = [CSClist_PFC;csclist_Rat011_PFC];
CSClist_PRL = [CSClist_PRL;csclist_Rat011_PRL];

rewardlast = [rewardlast;rewardlast_Rat011];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat011];

%% Rat016  海马和前额叶已经确认可以用
%  all hippocampus tetrode at CA1, all prefrontal except TT1,2 at Prl and Cg2
%  use tetrodes that have cells for every days' analyse
isession_Rat016 = 0;

% training
% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-13_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-13';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211213_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[14,16];
% csclist_Rat016_PFC{isession_Rat016,1}=[4,5,6,7,8];
% csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% rewardlast_Rat016(isession_Rat016,1)=12;
% rewardcurrent_Rat016(isession_Rat016,1)=9;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-14_CT\');
date_Rat016{isession_Rat016,1}='2021-12-14';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211214_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[12,14,16,17];
csclist_Rat016_PFC{isession_Rat016,1}=[6,7];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=9;
rewardcurrent_Rat016(isession_Rat016,1)=12;

%sample-test第一圈angle是空的
% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-15_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-15';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211215_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[14,16];
% csclist_Rat016_PFC{isession_Rat016,1}=[4,5,6,7,8];
% csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% rewardlast_Rat016(isession_Rat016,1)=12;
% rewardcurrent_Rat016(isession_Rat016,1)=8;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-16_CT\');
date_Rat016{isession_Rat016,1}='2021-12-16';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211216_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[12,14,15,16];
csclist_Rat016_PFC{isession_Rat016,1}=[5,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=8;
rewardcurrent_Rat016(isession_Rat016,1)=11;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-17_CT\');
date_Rat016{isession_Rat016,1}='2021-12-17';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211217_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[11,12,14,15,16];
csclist_Rat016_PFC{isession_Rat016,1}=[2,4,5,7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=11;
rewardcurrent_Rat016(isession_Rat016,1)=7;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-18_CT\');
date_Rat016{isession_Rat016,1}='2021-12-18';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211218_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[12,13,14,15,16];
csclist_Rat016_PFC{isession_Rat016,1}=[2,5];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=7;
rewardcurrent_Rat016(isession_Rat016,1)=10;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-20_CT\');
date_Rat016{isession_Rat016,1}='2021-12-20';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211220_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,13,14,15,16];
csclist_Rat016_PFC{isession_Rat016,1}=[2,5,7];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=10;
rewardcurrent_Rat016(isession_Rat016,1)=12;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-21_CT\');
date_Rat016{isession_Rat016,1}='2021-12-21';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211221_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[11,13,15,16,17];
csclist_Rat016_PFC{isession_Rat016,1}=[2,5];
csclist_Rat016_PRL{isession_Rat016,1}=[2,3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=12;
rewardcurrent_Rat016(isession_Rat016,1)=9;

%从今天开始cell变多了
isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-27_CT\');
date_Rat016{isession_Rat016,1}='2021-12-27';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211227_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[5,7];
csclist_Rat016_PRL{isession_Rat016,1}=[2,3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=9;
rewardcurrent_Rat016(isession_Rat016,1)=9;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-28_CT\');
date_Rat016{isession_Rat016,1}='2021-12-28';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211228_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[7];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=9;
rewardcurrent_Rat016(isession_Rat016,1)=7;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-29_CT\');
date_Rat016{isession_Rat016,1}='2021-12-29';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211229_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[2,7];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=7;
rewardcurrent_Rat016(isession_Rat016,1)=10;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-30_CT\');
date_Rat016{isession_Rat016,1}='2021-12-30';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211230_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[2,7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=10;
rewardcurrent_Rat016(isession_Rat016,1)=12;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-31_CT\');
date_Rat016{isession_Rat016,1}='2021-12-31';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211231_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[2,4,7];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=12;
rewardcurrent_Rat016(isession_Rat016,1)=8;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-02_CT\');
date_Rat016{isession_Rat016,1}='2022-01-02';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220102_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[2,3,4,7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=8;
rewardcurrent_Rat016(isession_Rat016,1)=11;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-04_CT\');
date_Rat016{isession_Rat016,1}='2022-01-04';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220104_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,17];
csclist_Rat016_PFC{isession_Rat016,1}=[7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=11;
rewardcurrent_Rat016(isession_Rat016,1)=13;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-05_CT\');
date_Rat016{isession_Rat016,1}='2022-01-05';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220105_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=13;
rewardcurrent_Rat016(isession_Rat016,1)=9;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-06_CT\');
date_Rat016{isession_Rat016,1}='2022-01-06';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220106_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,15,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[6,7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=9;
rewardcurrent_Rat016(isession_Rat016,1)=7;

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-07_CT\');
date_Rat016{isession_Rat016,1}='2022-01-07';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220107_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,16,17,18];
csclist_Rat016_PFC{isession_Rat016,1}=[2,4,7,8];
csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
rewardlast_Rat016(isession_Rat016,1)=7;
rewardcurrent_Rat016(isession_Rat016,1)=10;

Ind_Rat016 = ones(isession_Rat016,1)*016;

isession = isession+isession_Rat016;
Ind_Rat = [Ind_Rat;Ind_Rat016];
path = [path;path_Rat016];
date = [date;date_Rat016];
trackdata = [trackdata;trackdata_Rat016];
CSClist_CA1 = [CSClist_CA1;csclist_Rat016_CA1];
CSClist_PFC = [CSClist_PFC;csclist_Rat016_PFC];
CSClist_PRL = [CSClist_PRL;csclist_Rat016_PRL];
