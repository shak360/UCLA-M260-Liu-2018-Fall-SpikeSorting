%% Shakthi Visagan 804622954

% Professor Liu, M260 Neuroengineering
% EMG Spike Sorting/Decomposition
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
numTimeSteps_holder = size(time);
numTimeSteps = numTimeSteps_holder(1);
disp('number of time steps: ');
disp(numTimeSteps);


freq_samp = (time(2)-time(1))^(-1); % calculate the sample frequecy
disp('sampling frequency: [Hz]');
disp(freq_samp);

freq_Nyquist = freq_samp/2;
disp('Nyquist frequency: [Hz]');
disp(freq_Nyquist);

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
    xlabel('time [seconds]');
    ylabel('voltage [V]');

    title(str);
    xlim([time(1) time(end)]); % label and title each plots
    hold off;
end

channel_select = 1; % select channel for testing. channel_select <= channel_number
test_input = channelData(:,channel_select); % test_input will go through all the individual sections

%% Filtering the Signal by Removing Frequencies

signal_mean = mean(test_input);
disp('signal mean: ')
disp(signal_mean)

% Fourier analysis
y = fft(test_input-signal_mean); % identify the frequency component of the signal
n = length(test_input);          % number of samples
f = (0:n-1)*(freq_samp/n);       % frequency range
power = abs(y).^2/n;             % power of the DFT

% plotting power spectrum of mean centered signal
figure('Name','Power Spectrum of Mean Centered Signal','NumberTitle','off','Color','white');
fftPlot = plot(f,power);
fftPlot.Color = [0,0,1,1];
xlim([time(1) freq_Nyquist]);
xlabel('Frequency [Hz]');
ylabel('Power');
title('Power Spectrum of Mean Centered Signal');

% Cutoff Frequency: 300 to 3000 Hz
freq_lowerCutOff = 250;             % [Hz]
freq_upperCutOff = 800;  % [Hz]

% Bandpass Butterworth Filter 
[b,a] = butter(4, [freq_lowerCutOff/(freq_Nyquist),freq_upperCutOff/(freq_Nyquist)], 'bandpass');

% Filter the Signal
filt_sig = filtfilt(b, a, test_input); 

% plotting the filtered signal
figure('Name','Filtered Signal','NumberTitle','off','Color','white');
a = plot(time, filt_sig, time, test_input);
a(1).Color = [1,0,0,0.5]; % red
a(2).Color = [0,0,1,0.125]; % blue
legend('Filtered Signal', 'Original Signal')
title('Filtered EMG Signal');
ylabel('voltage [V]')
xlabel('time [seconds]');
xlim([time(1), time(end)]); 

filt_sig_mean = mean(filt_sig)  ;
disp('filtered signal mean: ')
disp(filt_sig_mean)

y_filt = fft(filt_sig-filt_sig_mean);     % identify the frequency component of the signal
n_filt = length(filt_sig);                % number of samples
f_filt = (0:n_filt-1)*(freq_samp/n_filt); % frequency range
power_filt = abs(y_filt).^2/n_filt;       % power of the DFT

% plotting power spectrum of mean centered signal
figure('Name','Power Spectrum of Mean Centered FILTERED Signal','NumberTitle','off','Color','white');
fft_filtPlot = plot(f_filt,power_filt);
fft_filtPlot.Color = [0,0,1,1];
xlim([time(1), freq_Nyquist]);
xlabel('Frequency [Hz]')
ylabel('Power')
title('Power Spectrum of Mean Centered FILTERED Signal') 

%% Detecting Spikes

% Quiroga et. al threshold
std_dev_estimate = median(abs(filt_sig)/0.6745);
MPH_Thr = 5*std_dev_estimate;
MPD = 2.5e-3; % width of a usual pulse
% or the minimum peak spread
 
[peaks,time_locs] = findpeaks(abs(filt_sig),'MinPeakHeight', MPH_Thr, 'MinPeakDistance',MPD);

detected_spike_vis =  zeros(numTimeSteps_holder);
detected_spike_vis(time_locs) = max(abs(filt_sig))/2;
threshold_vis_plus = zeros(numTimeSteps_holder);
threshold_vis_minus = threshold_vis_plus-MPH_Thr;
threshold_vis_plus = threshold_vis_plus+MPH_Thr;


% plotting the spikes that were detected as binary dirac impulses 
% from their time indices
figure('Name','Detecting Spikes','NumberTitle','off','Color','white');
b = plot(time, detected_spike_vis, time, filt_sig, time, threshold_vis_plus, '--', time, threshold_vis_minus, '--');
b(1).Color = [1,0,0,0.5]; % red
b(2).Color = [0,0,1,0.125]; % blue
b(3).Color = [0.25,0.25,0.25,1]; % black
b(4).Color = [0.25,0.25,0.25,1]; % black
legend('Detected Spikes', 'Filtered Original Signal', 'Threshold');
title('Spike Detected EMG Signal');
ylabel('voltage [V]')
xlabel('time [seconds]');
xlim([time(1), time(end)]); 

%% Aligning Spikes

% need a buffer region around the spikes 
point_dist_per_spike = MPD*(freq_samp);

index = 1;
for peak_i=1:length(peaks)
    for k = (loc(peak_i) - (int(point_dist_per_spike/2))):(loc(peak_i) + int(point_dist_per_spike/2))
        % Must sort it in detectedspikes
        if b < length(sigfilt) && b > 1
        %Store data of points collected
        detectedspikes(index, j) = sigfilt(b);
        index = index + 1;
        end 
    end 
    %Reset our pointer value
    index = 1;
end 

% final plotting of spikes
figure; 
plot(detectedspikes);
title('Spike Detected w/ Buffer Region EMG Signal');
ylabel('voltage [V]')
xlabel('time [seconds]');
xlim([time(1), time(end)]); 

