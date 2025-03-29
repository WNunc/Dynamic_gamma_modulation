function [ts,chNum,sampFreq,numSamp,header] = load_eeg_retrieve_header_cz(file)
%file = 'E:\Chenguang Zheng\DATA\Rat28\2013-09-11_09-18-22-LT\CSC9.ncs';

% Remember to load unsplit CSC files

% Set the field selection for reading CSC files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Channel number
fieldSelection(3) = 1; % Sample Frequency
fieldSelection(4) = 1; % Number of valid samples
fieldSelection(5) = 0; % Samples (EEG data)
% Do we return header, 1 = Yes, 0 = No.
extractHeader = 1;
% 5 different extraction modes, see help file for Nlx2MatCSC_v3

extractMode = 1; % Extract all data

fid_file = fopen(file,'r');
if fid_file < 0
    % fprintf('Header file not found, using default header;');
%     header{15}='-ADBitVolts 6.1037e-008';
    header = [];
    ts=0;
    chNum=0;
    sampFreq=0;
    numSamp=0;
else
    [ts,chNum,sampFreq,numSamp,header] =...
        Nlx2MatCSC(file,fieldSelection,extractHeader,extractMode);
end

% Transform the 2-D samples array to an 1-D array
%M = size(samp,2);
%samples = zeros(512*M,1);

%for jj = 1:M
    %samples(((jj-1)*512)+1:512*jj) = samp(1:512,jj);
%end

samples = 1;
samp = 1;