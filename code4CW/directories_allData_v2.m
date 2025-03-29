% for theta seqence with gamma CW

Dir = 'H:\neuralynx\';

isession = 0;
Ind_Rat = [];
path = {};
trackdata = {};
CSClist_CA1 = {};
csclist_CA1_mid = {};
Seq_cutPhase = {};% 第一列CW，第二列CCW
% CSClist_CA1_sort = {};
% CSClist_CA1_lowfir = {};
% CSClist_R2 = {};
% rewardlast = [];
% rewardcurrent = [];

%% Rat 002
isession_Rat002 = 0;
cutPhase_CW = 200;
cutPhase_CCW = 200;
% isession_Rat002 = isession_Rat002+1;
% path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-27-CT\');
% date_Rat002{isession_Rat002,1}='2020-03-27';
% trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-27.mat');
% csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
% % csclist_Rat002_CA2{isession_Rat002,1}=[5];
% % csclist_Rat002_PFC{isession_Rat002,1}=[];
% % csclist_Rat002_PRL{isession_Rat002,1}=[];
% % csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% % rewardlast_Rat002(isession_Rat002,1) = 8;%没有记录3-27之前的奖励点
% % rewardcurrent_Rat002(isession_Rat002,1) = 12;
% % 39 CA1 place cells
% % sign_correct_test = [0;0;1;0;0;0;0;NaN]

% isession_Rat002 = isession_Rat002+1;
% path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-30-CT\');
% date_Rat002{isession_Rat002,1}='2020-03-30';
% trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-30.mat');
% csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
% csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
% Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
% Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% % csclist_Rat002_CA2{isession_Rat002,1}=[5];
% % csclist_Rat002_PFC{isession_Rat002,1}=[];
% % csclist_Rat002_PRL{isession_Rat002,1}=[];
% % csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% % rewardlast_Rat002(isession_Rat002,1) = 9;
% % rewardcurrent_Rat002(isession_Rat002,1) = 7;
% % 54 CA1 place cells
% % % sign_correct_test = [1;0;1;0;1;0;NaN;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-03-31-CT\');
date_Rat002{isession_Rat002,1}='2020-03-31';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-03-31.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_PFC{isession_Rat002,1}=[];
% csclist_Rat002_PRL{isession_Rat002,1}=[];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 7;
% rewardcurrent_Rat002(isession_Rat002,1) = 11;
% 40 CA1 place cells
% % sign_correct_test = [NaN;1;0;1;1;0;1;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-01-CT\');
date_Rat002{isession_Rat002,1}='2020-04-01';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-01.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_PFC{isession_Rat002,1}=[];
% csclist_Rat002_PRL{isession_Rat002,1}=[];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 11;
% rewardcurrent_Rat002(isession_Rat002,1) = 8;
% 65 CA1 place cells
% % sign_correct_test = [0;0;1;1;0;0;NaN;1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-02-CT\');
date_Rat002{isession_Rat002,1}='2020-04-02';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-02.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_PFC{isession_Rat002,1}=[];
% csclist_Rat002_PRL{isession_Rat002,1}=[];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 8;
% rewardcurrent_Rat002(isession_Rat002,1) = 12;
% 52 CA1 place cells
% % sign_correct_test = [0;1;0;0;1;0;1;0]

isession_Rat002 = isession_Rat002+1; 
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-05-CT\');
date_Rat002{isession_Rat002,1}='2020-04-05';
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-05.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2,7,8,9,10,11,17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_PFC{isession_Rat002,1}=[];
% csclist_Rat002_PRL{isession_Rat002,1}=[];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS3R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 12;
% rewardcurrent_Rat002(isession_Rat002,1) = 9;
% 65 CA1 place cells
% % sign_correct_test = [1;1;1;1;0;1;NaN;0]

% isession_Rat002 = isession_Rat002+1;
% path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-06-CT\');
% trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-06.mat');
% csclist_Rat002_CA1{isession_Rat002,1}=[2 7 8 9 10 11 17];
% % csclist_Rat002_CA2{isession_Rat002,1}=[5];
% % csclist_Rat002_CA1_sort{isession_Rat002,1}=[6,17,12,16,1,14,9,8];
% % csclist_Rat002_CA1_lowfir{isession_Rat002,1}=[6,12,17];
% % csclist_Rat002_R2{isession_Rat002,1} = 'HS2R1.ncs';
% % rewardlast_Rat002(isession_Rat002,1) = 12;
% % rewardcurrent_Rat002(isession_Rat002,1) = 10;
% % 43 CA1 place cells
% % sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-04-08-CT\');
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-04-08.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2 7 8 9 10 11 17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_CA1_sort{isession_Rat002,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat002_CA1_lowfir{isession_Rat002,1}=[6,12,17];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS2R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 12;
% rewardcurrent_Rat002(isession_Rat002,1) = 10;
% 58 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat002 = isession_Rat002+1;
path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-05-05-CT\');
trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-05-05.mat');
csclist_Rat002_CA1{isession_Rat002,1}=[2 7 8 9 10 11 17];
csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% csclist_Rat002_CA2{isession_Rat002,1}=[5];
% csclist_Rat002_CA1_sort{isession_Rat002,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat002_CA1_lowfir{isession_Rat002,1}=[6,12,17];
% csclist_Rat002_R2{isession_Rat002,1} = 'HS2R1.ncs';
% rewardlast_Rat002(isession_Rat002,1) = 12;
% rewardcurrent_Rat002(isession_Rat002,1) = 10;
% 43 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

% isession_Rat002 = isession_Rat002+1;
% path_Rat002{isession_Rat002,1} = strcat(Dir,'Rat002\2020-05-07-CT\');
% trackdata_Rat002{isession_Rat002,1} = strcat(path_Rat002{isession_Rat002,1},'2020-05-07.mat');
% csclist_Rat002_CA1{isession_Rat002,1}=[2 7 8 9 10 11 17];
% csclist_Rat002_CA1_mid{isession_Rat002,1}=[8];
% Seq_cutPhase_Rat002{isession_Rat002,1} = cutPhase_CW;
% Seq_cutPhase_Rat002{isession_Rat002,2} = cutPhase_CCW;
% % csclist_Rat002_CA2{isession_Rat002,1}=[5];
% % csclist_Rat002_CA1_sort{isession_Rat002,1}=[6,17,12,16,1,14,9,8];
% % csclist_Rat002_CA1_lowfir{isession_Rat002,1}=[6,12,17];
% % csclist_Rat002_R2{isession_Rat002,1} = 'HS2R1.ncs';
% % rewardlast_Rat002(isession_Rat002,1) = 12;
% % rewardcurrent_Rat002(isession_Rat002,1) = 10;
% % 60 CA1 place cells
% % sign_correct_test = [1,0,0,1,0,0,0,1]

isession = isession+isession_Rat002;
Ind_Rat002 = ones(isession_Rat002,1)*2;
Ind_Rat = [Ind_Rat;Ind_Rat002];
path = [path;path_Rat002];
trackdata = [trackdata;trackdata_Rat002];
CSClist_CA1 = [CSClist_CA1;csclist_Rat002_CA1];
csclist_CA1_mid = [csclist_CA1_mid;csclist_Rat002_CA1_mid];
Seq_cutPhase = [Seq_cutPhase;Seq_cutPhase_Rat002];
% CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat002_CA1_sort];
% CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat002_CA1_lowfir];
% CSClist_R2 = [CSClist_R2;csclist_Rat002_R2];
% rewardlast = [rewardlast;rewardlast_Rat002];
% rewardcurrent = [rewardcurrent;rewardcurrent_Rat002];

%% Rat 003
isession_Rat003 = 0;

cutPhase_CW = 150;
cutPhase_CCW = 170;

% isession_Rat003 = isession_Rat003+1;
% path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-04-07-CT\');
% trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-04-07.mat');
% csclist_Rat003_CA1{isession_Rat003,1}=[2 4 5 6 7 8 9 10 12 13 14 15];
% csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
% Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
% Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% % csclist_Rat003_CA1_sort{isession_Rat003,1}=[10,12,7,9];
% % csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[10];
% % csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% % rewardlast_Rat003(isession_Rat003,1) = 9;
% % rewardcurrent_Rat003(isession_Rat003,1) = 14;
% % 62 CA1 place cells
% % sign_correct_test = [0,0,1,1,1,1,1,1]

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-04-30-CT\');
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-04-30.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17];
csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% csclist_Rat003_CA1_sort{isession_Rat003,1}=[10,12,7,9,8];
% csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[10];
% csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% rewardlast_Rat003(isession_Rat003,1) = 10;
% rewardcurrent_Rat003(isession_Rat003,1) = 15;
% 73 CA1 place cells
% sign_correct_test = [1,1,1,1,1,1,1,1]

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-05-03-CT\');
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-05-03.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[1 2 3 4 5 6 7 9 10 12 13 14 15 16];
csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% csclist_Rat003_CA1_sort{isession_Rat003,1}=[10,12,7,9];
% csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[7];
% csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% rewardlast_Rat003(isession_Rat003,1) = 8;
% rewardcurrent_Rat003(isession_Rat003,1) = 13;
% 75 CA1 place cells
% sign_correct_test = [0,0,1,1,1,1,1,1]

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-05-04-CT\');
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-05-04.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[1 2 3 4 5 6 7 9 10 13 14 15 16 17];
csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% csclist_Rat003_CA1_sort{isession_Rat003,1}=[10,12,7,9];
% csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[7];
% csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% rewardlast_Rat003(isession_Rat003,1) = 8;
% rewardcurrent_Rat003(isession_Rat003,1) = 13;
% 77CA1 place cells
% sign_correct_test = [0,0,1,1,1,1,1,1]

isession_Rat003 = isession_Rat003+1;
path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-06-15-CT\');
trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-06-15.mat');
csclist_Rat003_CA1{isession_Rat003,1}=[2 4 6 9 10 12 13 14 17];
csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% csclist_Rat003_CA1_sort{isession_Rat003,1}=[10,12,7,9];
% csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[7];
% csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% rewardlast_Rat003(isession_Rat003,1) = 8;
% rewardcurrent_Rat003(isession_Rat003,1) = 13;
% 37 CA1 place cells, 1 CA1 interneurons
% sign_correct_test = [0,0,1,1,1,1,1,1]

% isession_Rat003 = isession_Rat003+1;
% path_Rat003{isession_Rat003,1} = strcat(Dir,'Rat003\2020-06-17-CT\');
% trackdata_Rat003{isession_Rat003,1} = strcat(path_Rat003{isession_Rat003,1},'2020-06-17.mat');
% csclist_Rat003_CA1{isession_Rat003,1}=[2 6 7 9 10 13 14 17];
% csclist_Rat003_CA1_mid{isession_Rat003,1}=[6];
% Seq_cutPhase_Rat003{isession_Rat003,1} = cutPhase_CW;
% Seq_cutPhase_Rat003{isession_Rat003,2} = cutPhase_CCW;
% % csclist_Rat003_CA1_lowfir{isession_Rat003,1}=[10];
% % csclist_Rat003_R2{isession_Rat003,1} = 'HS3R1.ncs';
% % rewardlast_Rat003(isession_Rat003,1) = 13;
% % rewardcurrent_Rat003(isession_Rat003,1) = 9;
% % 36 CA1 place cells
% % sign_correct_test = [0,1,1,0,1,1,1,1]
% % NOTE: delete this data because of rare behavior data in error trials

isession = isession+isession_Rat003;
Ind_Rat003 = ones(isession_Rat003,1)*3;
Ind_Rat = [Ind_Rat;Ind_Rat003];
path = [path;path_Rat003];
trackdata = [trackdata;trackdata_Rat003];
CSClist_CA1 = [CSClist_CA1;csclist_Rat003_CA1];
csclist_CA1_mid = [csclist_CA1_mid;csclist_Rat003_CA1_mid];
Seq_cutPhase = [Seq_cutPhase;Seq_cutPhase_Rat003];
% CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat003_CA1_sort];
% CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat003_CA1_lowfir];
% CSClist_R2 = [CSClist_R2;csclist_Rat003_R2];
% rewardlast = [rewardlast;rewardlast_Rat003];
% rewardcurrent = [rewardcurrent;rewardcurrent_Rat003];

%% Rat 010
% isession_Rat010 = 0;
% 
% cutPhase_CW = 130;
% cutPhase_CCW = 130;

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-06_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-06';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210706_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';
% % 38 CA1 place cells

% isession_Rat010 = isession_Rat010+1; 
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-06-30_CT\');
% date_Rat010{isession_Rat010,1} = '2021-06-30';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210630_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; 
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-01_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-01';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210701_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-04_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-04';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210704_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-05_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-05';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210705_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-09_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-09';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210709_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-12_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-12';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210712_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-12_CT_2\');
% date_Rat010{isession_Rat010,1} = '2021-07-12-2';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210712_CT_2_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-13_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-13';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210713_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-13_CT_2\');
% date_Rat010{isession_Rat010,1} = '2021-07-13-2';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210713_CT_2_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';

% isession_Rat010 = isession_Rat010+1; %Day
% path_Rat010{isession_Rat010,1} = strcat(Dir,'Rat010\2021-07-15_CT\');
% date_Rat010{isession_Rat010,1} = '2021-07-15';
% trackdata_Rat010{isession_Rat010,1} = strcat(path_Rat010{isession_Rat010,1},'20210715_CT_tracking.mat');
% csclist_Rat010_CA1{isession_Rat010,1}=[10,11,13,14,15];
% csclist_Rat010_CA1_mid{isession_Rat010,1}=[13];
% Seq_cutPhase_Rat010{isession_Rat010,1} = cutPhase_CW;
% Seq_cutPhase_Rat010{isession_Rat010,2} = cutPhase_CCW;
% % csclist_Rat010_CA2{isession_Rat010,1}=[17,18];
% % csclist_Rat010_CA3{isession_Rat010,1}=[];
% % csclist_Rat010_PFC_mid{isession_Rat010,1}=[];
% % csclist_Rat010_PFC{isession_Rat010,1}=[];
% % csclist_Rat010_PRL{isession_Rat010,1}=[];
% % csclist_Rat010_R2{isession_Rat010,1} = 'HS2R1.ncs';


% isession = isession+isession_Rat010;
% Ind_Rat010 = ones(isession_Rat010,1)*010;
% Ind_Rat = [Ind_Rat;Ind_Rat010];
% path = [path;path_Rat010];
% date = [date;date_Rat010];
% trackdata = [trackdata;trackdata_Rat010];
% CSClist_CA1 = [CSClist_CA1;csclist_Rat010_CA1];
% csclist_CA1_mid = [csclist_CA1_mid;csclist_Rat010_CA1_mid];
% Seq_cutPhase = [Seq_cutPhase;Seq_cutPhase_Rat010];
% % CSClist_PFC = [CSClist_PFC;csclist_Rat010_PFC];
% % CSClist_PRL = [CSClist_PRL;csclist_Rat010_PRL];
% % CSClist_HPC_R = [CSClist_HPC_R;csclist_Rat010_R2];
%% Rat 011
isession_Rat011 = 0;

cutPhase_CW = 130;
cutPhase_CCW = 170;

% isession_Rat011 = isession_Rat011+1;
% path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-20_CT_9\');
% trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-20-CT.mat');
% csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
% csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
% Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
% Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% % csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% % csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% % csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% % rewardlast_Rat011(isession_Rat011,1) = 12;
% % rewardcurrent_Rat011(isession_Rat011,1) = 10;
% % 71 CA1 place cells
% % sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-21_CT_11\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-21-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 62 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-22_CT_11\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-22-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 56 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]


isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-23_CT_11\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-23-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 72 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-25_CT_9\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-25-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 66 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-26_CT_8\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-26-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 74 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat011 = isession_Rat011+1;
path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-27_CT_6\');
trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-27-CT.mat');
csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% rewardlast_Rat011(isession_Rat011,1) = 12;
% rewardcurrent_Rat011(isession_Rat011,1) = 10;
% 76 CA1 place cells
% sign_correct_test = [1,0,0,1,0,0,0,1]

% isession_Rat011 = isession_Rat011+1;
% path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-28_CT_PM_7\');
% trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-28-CT-PM.mat');
% csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
% csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
% Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
% Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% % csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% % csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% % csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% % rewardlast_Rat011(isession_Rat011,1) = 12;
% % rewardcurrent_Rat011(isession_Rat011,1) = 10;
% % 73 CA1 place cells
% % sign_correct_test = [1,0,0,1,0,0,0,1]

% isession_Rat011 = isession_Rat011+1;
% path_Rat011{isession_Rat011,1} = strcat(Dir,'Rat011\2021-07-29_CT_AM_15\');
% trackdata_Rat011{isession_Rat011,1} = strcat(path_Rat011{isession_Rat011,1},'2021-07-29-CT-AM.mat');
% csclist_Rat011_CA1{isession_Rat011,1}=[10 11 12 13 14 15 16 17 18];
% csclist_Rat011_CA1_mid{isession_Rat011,1}=[17];
% Seq_cutPhase_Rat011{isession_Rat011,1} = cutPhase_CW;
% Seq_cutPhase_Rat011{isession_Rat011,2} = cutPhase_CCW;
% % csclist_Rat011_CA1_sort{isession_Rat011,1}=[6,17,12,16,1,14,9,8];
% % csclist_Rat011_CA1_lowfir{isession_Rat011,1}=[6,12,17];
% % csclist_Rat011_R2{isession_Rat011,1} = 'HS2R1.ncs';
% % rewardlast_Rat011(isession_Rat011,1) = 12;
% % rewardcurrent_Rat011(isession_Rat011,1) = 10;
% % 64 CA1 place cells
% % sign_correct_test = [1,0,0,1,0,0,0,1]

isession = isession+isession_Rat011;
Ind_Rat011 = ones(isession_Rat011,1)*11;
Ind_Rat = [Ind_Rat;Ind_Rat011];
path = [path;path_Rat011];
trackdata = [trackdata;trackdata_Rat011];
CSClist_CA1 = [CSClist_CA1;csclist_Rat011_CA1];
csclist_CA1_mid = [csclist_CA1_mid;csclist_Rat011_CA1_mid];
Seq_cutPhase = [Seq_cutPhase;Seq_cutPhase_Rat011];
% CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat011_CA1_sort];
% CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat011_CA1_lowfir];
% CSClist_R2 = [CSClist_R2;csclist_Rat011_R2];
% rewardlast = [rewardlast;rewardlast_Rat011];
% rewardcurrent = [rewardcurrent;rewardcurrent_Rat011];

%% Rat016
isession_Rat016 = 0;

cutPhase_CW = 140;
cutPhase_CCW = 140;

%从今天开始cell变多了
% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-27_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-27';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211227_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,16,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[5,7];
% % csclist_Rat016_PRL{isession_Rat016,1}=[2,3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=9;
% % rewardcurrent_Rat016(isession_Rat016,1)=9;
% % 49 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-28_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-28';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211228_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[7];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=9;
% % rewardcurrent_Rat016(isession_Rat016,1)=7;
% % 39 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-29_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-29';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211229_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,17,18];
% % csclist_Rat016_PFC{isession_Rat016,1}=[2,7];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=7;
% % rewardcurrent_Rat016(isession_Rat016,1)=10;
% % 36 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-30_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-30';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211230_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,16,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[2,7,8];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=10;
% % rewardcurrent_Rat016(isession_Rat016,1)=12;
% % % 50 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2021-12-31_CT\');
% date_Rat016{isession_Rat016,1}='2021-12-31';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20211231_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,16,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[2,4,7];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=12;
% % rewardcurrent_Rat016(isession_Rat016,1)=8;
% % 45 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-02_CT\');
% date_Rat016{isession_Rat016,1}='2022-01-02';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220102_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[2,3,4,7,8];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=8;
% % rewardcurrent_Rat016(isession_Rat016,1)=11;
% % 56 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-04_CT\');
% date_Rat016{isession_Rat016,1}='2022-01-04';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220104_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,15,17];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[7,8];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=11;
% % rewardcurrent_Rat016(isession_Rat016,1)=13;
% % 43 CA1 place cells

isession_Rat016 = isession_Rat016+1;
path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-05_CT\');
date_Rat016{isession_Rat016,1}='2022-01-05';
trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220105_CT_tracking.mat');
csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,13,14,16,17,18];
csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% csclist_Rat016_PFC{isession_Rat016,1}=[7,8];
% csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% rewardlast_Rat016(isession_Rat016,1)=13;
% rewardcurrent_Rat016(isession_Rat016,1)=9;
% 49 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-06_CT\');
% date_Rat016{isession_Rat016,1}='2022-01-06';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220106_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,15,16,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[6,7,8];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=9;
% % rewardcurrent_Rat016(isession_Rat016,1)=7;
% % 62 CA1 place cells

% isession_Rat016 = isession_Rat016+1;
% path_Rat016{isession_Rat016,1} = strcat(Dir,'Rat016\2022-01-07_CT\');
% date_Rat016{isession_Rat016,1}='2022-01-07';
% trackdata_Rat016{isession_Rat016,1} = strcat(path_Rat016{isession_Rat016,1},'20220107_CT_tracking.mat');
% csclist_Rat016_CA1{isession_Rat016,1}=[10,11,12,14,16,17,18];
% csclist_Rat016_CA1_mid{isession_Rat016,1}=[14];
% Seq_cutPhase_Rat016{isession_Rat016,1} = cutPhase_CW;
% Seq_cutPhase_Rat016{isession_Rat016,2} = cutPhase_CCW;
% % csclist_Rat016_PFC{isession_Rat016,1}=[2,4,7,8];
% % csclist_Rat016_PRL{isession_Rat016,1}=[3,4,5,6,7,8,9];
% % csclist_Rat016_R2{isession_Rat016,1}='HS3R1.ncs';
% % csclist_Rat016_R1{isession_Rat016,1}='HS1R1.ncs';
% % rewardlast_Rat016(isession_Rat016,1)=7;
% % rewardcurrent_Rat016(isession_Rat016,1)=10;
% % 55 CA1 place cells

isession = isession+isession_Rat016;
Ind_Rat016 = ones(isession_Rat016,1)*16;
Ind_Rat = [Ind_Rat;Ind_Rat016];
path = [path;path_Rat016];
trackdata = [trackdata;trackdata_Rat016];
CSClist_CA1 = [CSClist_CA1;csclist_Rat016_CA1];
csclist_CA1_mid = [csclist_CA1_mid;csclist_Rat016_CA1_mid];
Seq_cutPhase = [Seq_cutPhase;Seq_cutPhase_Rat016];
% CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat016_CA1_sort];
% CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat016_CA1_lowfir];
% CSClist_R2 = [CSClist_R2;csclist_Rat016_R2];
% rewardlast = [rewardlast;rewardlast_Rat016];
% rewardcurrent = [rewardcurrent;rewardcurrent_Rat016];