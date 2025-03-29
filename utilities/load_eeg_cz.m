function [ts,chNum,sampFreq,numSamp,samp,samples] = load_eeg_cz(file)
%file = 'E:\Chenguang Zheng\DATA\Rat28\2013-09-11_09-18-22-LT\begin1\CSC9.ncs';

% Set the field selection for reading CSC files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Channel number
fieldSelection(3) = 1; % Sample Frequency
fieldSelection(4) = 1; % Number of valid samples
fieldSelection(5) = 1; % Samples (EEG data)
% Do we return header, 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatCSC_v3

extractMode = 1; % Extract all data

[ts,chNum,sampFreq,numSamp,samp] =...
   Nlx2MatCSC(file,fieldSelection,extractHeader,extractMode);

% Transform the 2-D samples array to an 1-D array
M = size(samp,2);
samples = zeros(512*M,1);

for jj = 1:M
    samples(((jj-1)*512)+1:512*jj) = samp(1:512,jj);
end
end