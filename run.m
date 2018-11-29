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

csvSize = size(M);
disp('CSV file rows: ');
disp(csvSize(1));
disp('CSV file columns: ');
disp(csvSize(2));

time = M(:,1); % first column is the time series
numTimeSteps = size(time);
disp('number of time steps: ');
disp(numTimeSteps);


freq_samp = (time(2)-time(1))^(-1); % calculate the sample frequecy
disp('sampling frequency: [Hz]');
disp(freq_samp);

numChannels = csvSize(2)-1; % num of channels in the database
channelData = M;
channelData(:,1)= []; % creating matrix with only channel data


% Plotting 
for i=1:numChannels
    str= sprintf('Channel %d',i);
    disp('working on ... ');
    disp(str)
    figure('Name',str,'NumberTitle','off','Color','white');
    p = plot(time, channelData(:,i), 'LineWidth', 1); %plot each channel
    p.Color = [0,0,1,0.125];
    xlabel('seconds');
    title(str);
    xlim([time(1) time(end)]); % label and title each plots
    hold off
end

channel_select = 1; % select channel for testing. channel_select <= channel_number
test_input = channelData(:,channel_select); % test_input will go through all the individual sections

%% Filtering the Signal by Removing Frequencies

y = fft(test_input); % identify the frequency component of the signal
fftPlot = plot(f,power);



%Cutoff Frequency: 300 to 3000 Hz
freq_lowerCutOff = 300; % [Hz]
freq_upperCutOff = 3000; % [Hz]

%Bandpass Butterworth Filter 
%butter(Order,Wn(lower), wn(upper), 'fn'); Wn=Cutof Frequency
[b,a] = butter(4, [freq_lowerCutOff/(freq_samp/2), freq_upperCutOff/(freq_samp/2)], 'bandpass');

%Obtain the frequencies from the current filter coefficients
freqz(b,a,8192); 

%Filter the Signal
filt_sig=filtfilt(b, a, test_input); 

%plotting the filtered signal
figure('Name','Filtered Signal','NumberTitle','off','Color','white');
a = plot(time, filt_sig, time, test_input);
a(1).Color = [0,0,1,0.125];
a(2).Color = [1,0,0,0.125];
title('Filtered EMG Signal');
ylabel('Voltage (V)');
xlabel('Time (s)');
xlim([time(1) time(size(time,1))]); 