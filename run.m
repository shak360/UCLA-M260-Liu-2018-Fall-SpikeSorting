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

%M = csvread('EMG_example_2_fs_2k.csv'); %read in csv file
%M = csvread('EMG_example_1_90s_fs_2k.csv');
M = csvread('EMG_example_20s_2000Hz-2016.csv');

csvSize = size(M);
disp('CSV file rows: ');
disp(csvSize(1));
disp('CSV file columns: ');
disp(csvSize(2));

% uncomment when time is the first column

%time = M(:,1); % first column is the time series
%numTimeSteps_holder = size(time);
%numTimeSteps = numTimeSteps_holder(1);
%disp('number of time steps: ');
%disp(numTimeSteps);

% uncomment when time is the first column

%freq_samp = (time(2)-time(1))^(-1); % calculate the sample frequecy
freq_samp = 2000;
disp('sampling frequency: [Hz]');
disp(freq_samp);

% comment when time is the first column
numTimeSteps_holder = size(M(:,1));
numTimeSteps = numTimeSteps_holder(1);
time = linspace(0, numTimeSteps/freq_samp, numTimeSteps); % calulate the time manually

freq_Nyquist = freq_samp/2;
disp('Nyquist frequency: [Hz]');
disp(freq_Nyquist);

% uncomment when time is the first column
%numChannels = csvSize(2)-1; % num of channels in the database
% comment when time is the first column
numChannels = 2; % need to rewrite this manually for some csvs
channelData = M;
% uncomment when time is the first column
%channelData(:,1)= []; % creating matrix with only channel data


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

%% Filter Signal: Filtering the Signal by Removing Frequencies

signal_mean = mean(test_input); % calculate the signal mean
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
hold off;

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
hold off;


filt_sig_mean = mean(filt_sig)  ; %calculate the mean of the filtered signal
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
hold off;

%% Detect Spikes: Detecting Spikes

% Quiroga et. al threshold
std_dev_estimate = median(abs(filt_sig)/0.6745);
MPH_Thr = 6*std_dev_estimate;
MPD = 5e-3; % width of a usual pulse
% or the minimum peak spread
 
[peaks,time_locs] = findpeaks(abs(filt_sig),'MinPeakHeight', MPH_Thr, 'MinPeakDistance',MPD);

detected_spike_vis =  zeros(numTimeSteps_holder); %dummy array to visualize binary spikes
detected_spike_vis(time_locs) = max(abs(filt_sig))/2;
threshold_vis_plus = zeros(numTimeSteps_holder);
threshold_vis_minus = threshold_vis_plus-MPH_Thr; %lines representing the threshold we set above
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
hold off;


%% Align Spikes: Aligning Spikes

% need a buffer region around the spikes 
point_dist_per_spike = MPD*(40*freq_samp);
desired_point_dist_per_spike = round(point_dist_per_spike/4);


detected_spikes_holder = zeros(length(peaks), desired_point_dist_per_spike);

for peak_i=1:length(peaks)
    begin = time_locs(peak_i) - point_dist_per_spike/2;
    if begin <1
        begin = 1;
    end
    endin = time_locs(peak_i) + point_dist_per_spike/2;
    if endin > numTimeSteps
        endin = numTimeSteps;
    end
    spike_spline = spline(time(begin:endin), filt_sig(begin:endin), linspace(time(begin),time(endin), desired_point_dist_per_spike));
    detected_spikes_holder(peak_i, :) = spike_spline;
end 

% final plotting of spikes
figure('Name','Detecting and Aligning Spikes','NumberTitle','off','Color','white');
for peak_i=1:length(peaks)
    plot(detected_spikes_holder(peak_i,:))
    xlim([1,desired_point_dist_per_spike]);
    hold on;
end
title('Spike Detected w/ Buffer Region and Cubic Spline Interpolation EMG Signal');
ylabel('voltage [V]')
hold off;


%% Extract Features: t-SNE Latent Space Measurements

rng default
[Y3,loss3] = tsne(detected_spikes_holder,'NumDimensions',3, 'Distance', 'correlation');
disp('t-SNE loss: ')
disp(loss3)

figure('Name','t-SNE','NumberTitle','off','Color','white');
scatter3(Y3(:,1),Y3(:,2),Y3(:,3))
title('3-D Embedding t-SNE for Spikes')
view(-50,8)
axis equal
hold off;


%% Extract Features: PCA Dimensionality Reduction / Feature Extraction

[coeff,score,latent,~,explained] = pca(detected_spikes_holder);
disp('Percent Variance explained: ')
disp(explained)

figure('Name','PCA','NumberTitle','off','Color','white');
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title('PCA Dimensionality Reduction of Spikes')
hold off;

pca_score_holder = score(:,1:3);

%% Cluster Spikes: k-Means Clustering

[k_means_prim_idx] = kmeans(pca_score_holder,2, 'Distance', 'cityblock');
spike_clust1 = detected_spikes_holder(k_means_prim_idx==1,:);
spike_clust2 = detected_spikes_holder(k_means_prim_idx==2,:);

figure('Name','Initial Clustering','NumberTitle','off','Color','white');
for peak_i=1:length(spike_clust1)
    c = plot(spike_clust1(peak_i,:), 'm');
    c.Color= [1,0,1,0.125];
    xlim([1,desired_point_dist_per_spike]);
    hold on;
end
hold on
for peak_i=1:length(spike_clust2)
    d = plot(spike_clust2(peak_i,:), 'g');
    d.Color= [0,1,0,0.125];
    xlim([1,desired_point_dist_per_spike]);
    hold on;
end
title('Initial Cluster Assignments')
legend('Cluster 1', 'Cluster 2');
hold off

figure('Name','Centroids','NumberTitle','off','Color','white');
e = plot(mean(spike_clust1), 'm');
hold on;
f = plot(mean(spike_clust2), 'g');
xlim([1,desired_point_dist_per_spike]);
title('Centroids of initial k-Means')
legend('Cluster 1 centroid', 'Cluster 2 centroid');
hold off;

%% Classify Spikes: k-Means Clustering 2

% elbow method
% evaluate the best number of clusters to take by cutoffing percent
% variance explained
[m,~] = size(spike_clust1);
ToTest=ceil(sqrt(m));
D=zeros(ToTest,1); %initialize the results matrix
for c=1:ToTest %for each sample
    [~,~,dist]=kmeans(spike_clust1,c,'emptyaction','drop','Distance', 'cityblock'); %compute the sum of intra-cluster distances
    tmp=sum(dist); %best so far
    
    for cc=2:3 %repeat the algo
        [~,~,dist]=kmeans(spike_clust1,c,'emptyaction','drop','Distance', 'cityblock');
        tmp=min(sum(dist),tmp);
    end
    D(c,1)=tmp; %collect the best so far in the results vecor
end


Var=D(1:end-1)-D(2:end); %calculate %variance explained
PC=cumsum(Var)/(D(1)-D(end));
[r,~]=find(PC>0.95); %find the best index
K=1+r(1,1); %get the optimal number of clusters
K = round(K/10);
[IDX1,C1,~]=kmeans(spike_clust1,K,'Distance', 'cityblock'); %now rerun one last time with the optimal number of clusters

figure('Name','Centroids of Cluster 1','NumberTitle','off','Color','white');
plotcolormap = jet(K);
for peak_i=1:K
    blah = C1(peak_i,:);
    blah_short = blah(45:55);
    blah_spline = spline(45:55, blah_short, 45:0.1:55);
    g = plot(blah_spline);
    g.Color = plotcolormap(peak_i,:);
    hold on;
end
title('Sorted Spikes in Cluster 1')
hold off;

%% Analysis

figure('Name','Final Clustering','NumberTitle','off','Color','white');
for k=1:K % for each cluster we've identified
    spike_clist_fin = spike_clust1(IDX1==k,:); % find the cluster's spikes
    [sizer,~] = size(spike_clist_fin); % get the length of the cluster
    subplot(3,2,k) % instantiate subplots
    for peak_i=1:sizer
        c = plot(spike_clist_fin(peak_i,:)); %plot the spike according to its color
        c.Color = plotcolormap(k,:);        
        xlim([1,desired_point_dist_per_spike]);
        hold on;
    end
    title('Final Cluster Assignments')
end
hold off;
    

