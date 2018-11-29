%% Shakthi Visagan 804622954
% Professor Liu, M260 Neuroengineering
% Spike Sorting/Decomposition
% 30 November, 2018 

%% Administrative Code
clc; 
clear all;
close all;
format compact;

%% Input File
M = csvread('EMG_example_2_fs_2k.csv'); %read in csv file
time = M(:,1); % first column is the time series
freq_samp = (time(2)-time(1))^(-1); % calculate the sample frequecy
channel_number = size(M,2)-1; % num of channels in the database

% Plotting 
for i=1:channel_number
    figure('Name', num2str(i), 'Color',[1,1,1]);
    plot(time, M(:,i+1)); %plot each channel
    str= sprintf('Channel %d',i);
    xlabel('seconds');
    title(str);
    xlim([time(1) time(size(time,1))]); % label and title each plots
end

channel_select = 1; % select channel for testing. channel_select <= channel_number
test_input = M(:,channel_select+1); % test_input will go through all the individual sections

%% Filter the Signal by Removing Frequencies

%Cutoff Frequency: 300 to 10kHz
freq_lowerCutOff = 200; % [Hz]
freq_upperCutOff = 800; % [Hz]

%Sampling Rate-double the highest cutoff frequency, generally 40kHz
%Bandpass Butterworth Filter the Signal 
[b,a] = butter(4,[freq_lowerCutOff/(freq_samp/2),freq_upperCutOff/(freq_samp/2)], 'bandpass'); %butter(Order,Wn(lower), wn(upper), 'fn'); Wn=Cutof Frequency

%Obtain the frequencies from the current filter coefficients
freqz(b,a,8192); 

%Filter the Signal
sigfilt=filtfilt(b, a, test_input); 

%Visualize the picture
figure();
plot(time, sigfilt, 'r', time, test_input, 'b');
title('Filtered EMG Signal');
ylabel('Voltage (V)');
xlabel('Time (s)');
xlim(['auto']);
xlim([time(1) time(size(time,1))]); 