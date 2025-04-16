%% add paths to field trip and libeep(to read events)
addpath('C:\libeep-3.3.177-matlab');
addpath ('C:\fieldtrip-20240113')
ft_defaults;
%DIR         = 'C:\Users\z5171263\Downloads\plot_ecg';
% DIR= 'C:\Users\sidpa\Downloads\plot_ecg';
DIR = 'D:\plot_ecg';
FILE        = '31_03_25_SIPA_TMS 1.cnt';
fn    = fullfile(DIR,FILE);
info = eepv4_read_info(fn);

%% take events and only keep 10,-10 (start and end of a 5 minute block)
allLabels = {info.triggers.label};
keepIdx = ismember(allLabels, {'10', '-10'});
info.triggers = info.triggers(keepIdx);
%% set up config, only interested in BIP5 (ecg channel)

cfg = [];
cfg.dataset = 'D:\plot_ecg\31_03_25_SIPA_TMS 1.cnt';
cfg.channel = {'BIP5'};        
cfg.bpfilter = 'yes';
cfg.bpfreq   = [0.5 40];
cfg.bsfilter = 'yes';
cfg.bsfreq   = [49 51];
cfg.demean   = 'yes';
data = ft_preprocessing(cfg);
%% plot whole ecg signal vs. time 

ecg_signal = data.trial{1};
t = data.time{1};
plot(t, ecg_signal);
fs= data.fsample;

%% can also just plot a small chunk

%adjust value of n_samples for bigger chunks
n_samples = 1000;
time_vector = (0:n_samples-1)/fs;

ecg_subset = ecg_signal(1:n_samples);
figure;
plot(time_vector, ecg_subset);
xlabel('Time (s)');
ylabel('ECG Amplitude');

%% plot whole ecg signal with start and end blocks marked

figure;
plot(t, ecg_signal, 'b'); 
hold on;
title('Filtered ECG with 10/-10 Trigger Markers');
xlabel('Time (s)');
ylabel('ECG Amplitude');

% Loop through triggers and plot markers
for i = 1:length(info.triggers)
    label = info.triggers(i).label;
    trigger_time = info.triggers(i).seconds_in_file;

    % Only plot for labels '10' and '-10'
    if strcmp(label, '10') || strcmp(label, '-10')
        % Get the sample index corresponding to this time
        trig_sample = round(trigger_time * fs);

        % Choose marker color based on label
        if strcmp(label, '10')
            col = 'g';  % green for '10'
        else
            col = 'm';  % magenta for '-10'
        end

        % Plot the dot at the appropriate sample
        if trig_sample > 0 && trig_sample <= length(ecg_signal)
            plot(t(trig_sample), ecg_signal(trig_sample), 'o', ...
                 'Color', col, 'MarkerSize', 6, 'MarkerFaceColor', col);
        end
    end
end

legend('Filtered ECG', 'Trigger: 10', 'Trigger: -10');

%% segment data based on the blocks

% Create array of the 6 blocks with corresponding time vectors
ecg_blocks = {};
time_blocks = {};

for i = 1:length(info.triggers)
    label = info.triggers(i).label;
    trigger_time = info.triggers(i).seconds_in_file;

    if strcmp(label, '10')  % Start of block
        start_sample = round(trigger_time * fs);

        for k = i+1:length(info.triggers)
            if strcmp(info.triggers(k).label, '-10')  % End of block
                end_sample = round(info.triggers(k).seconds_in_file * fs);

                if start_sample > 0 && end_sample <= length(ecg_signal)
                    ecg_blocks{end+1} = ecg_signal(start_sample:end_sample); % ECG segment
                    time_blocks{end+1} = t(start_sample:end_sample);         % Corresponding time vector
                end

                break; % Stop looking after first '-10'
            end
        end
    end
end

% display duration of each block (should be 300s each= 5minute blocks)
for k = 1:length(time_blocks)
    duration = time_blocks{k}(end) - time_blocks{k}(1);
    fprintf('Block %d duration: %.2f seconds\n', k, duration);
end

%% plot each block just to visualise
%manually adjust which block....need to label blocks later for easier use..
plot(time_blocks{2}, ecg_blocks{2});
title('ECG Block number...');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
%% plot all 6 in a subplot

num_blocks = length(ecg_blocks);
figure;
for i = 1:num_blocks
    subplot(num_blocks, 1, i);  % Rows = total blocks, 1 column
    plot(time_blocks{i}, ecg_blocks{i}, 'b');
    title(sprintf('ECG Block %d', i));
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end
sgtitle('All ECG Blocks');

%% more segmentation (only for the first 100s)

start_time_t = 0;        % Start time (in seconds)
end_time_t = 100;       % End time (in seconds)

% Find the indices corresponding to the start and end times
start_idx_t = find(time_blocks{1} >= start_time_t, 1, 'first');
end_idx_t = find(time_blocks{1} <= end_time_t, 1, 'last');

% Slice the data for the desired time range
ecg_segment_100s = ecg_blocks{1}(start_idx_t:end_idx_t);
time_segment_100s = time_blocks{1}(start_idx_t:end_idx_t);

% Plot the selected segment
figure;
plot(time_segment_100s, ecg_segment_100s);
title('ECG Block Segment (100s)');
xlabel('Time (s)');
ylabel('ECG Amplitude');


%% RR for selected segment (100s)
%amp threshold and time between each peak lines are from stackoverflow, not
%sure how valid need to test more

ecg_segment_100s = detrend(ecg_segment_100s);
 % Auto threshold estimation (e.g., 35% of peak-to-peak amplitude)
amp_thresh = 0.35 * (max(ecg_segment_100s) - min(ecg_segment_100s));
% Detect R-peaks, if two peaks closer than 300 samples then it is not
% considered (again need to check if this is valid)
[~, locs] = findpeaks(ecg_segment_100s, 'MinPeakHeight', amp_thresh,'MinPeakDistance', round(0.6 * fs)); 
% Time of R-peaks
r_times_100s = time_segment_100s(locs);
% RR intervals (in seconds), will be used for the stats later 
rr_intervals_100s = diff(r_times_100s);
% Plot
figure;
plot(time_segment_100s, ecg_segment_100s);
hold on;
plot(time_segment_100s(locs), ecg_segment_100s(locs), 'ro');
title('R-peaks in 100s ECG Segment');
xlabel('Time (s)');
ylabel('Amplitude');

%% R peak detection for all blocks
%amp threshold and time between each peak lines are from stackoverflow, not
%sure how valid need to test more

num_blocks = length(ecg_blocks);
all_rr = cell(1, num_blocks);  % Store RR intervals
for k = 1:num_blocks
    signal = ecg_blocks{k};
    t = time_blocks{k};
    % normalise signal
    signal = detrend(signal);
    amp_thresh = 0.35 * (max(signal) - min(signal));
    [~, locs] = findpeaks(signal, 'MinPeakHeight', amp_thresh,'MinPeakDistance', round(0.6 * fs));  % ~100 bpm max
    % Time of R-peaks
    r_times_all_blocks = t(locs);
    % RR intervals
    rr_intervals_all_blocks = diff(r_times_all_blocks);
    all_rr{k} = rr_intervals_all_blocks;

    figure;
    plot(t, signal);
    hold on;
    plot(t(locs), signal(locs), 'ro');
    title(sprintf('Block %d R-peaks', k));
    xlabel('Time (s)');
    ylabel('Amplitude');
end


%% scatter, doesn't show much but thought maybe we can do sth with this later
for k = 1:num_blocks
    signal = ecg_blocks{k};
    t = time_blocks{k};
    signal = detrend(signal);  % Normalize signal

    % Find R-peaks
    amp_thresh = 0.35 * (max(signal) - min(signal));
    [~, locs] = findpeaks(signal, 'MinPeakHeight', amp_thresh, 'MinPeakDistance', round(0.6 * fs));
    r_times_all_blocks = t(locs);  % Time of R-peaks
    r_amplitudes = signal(locs);  % Amplitude at R-peak locations

    % Plot R-peaks as a scatter plot of R-peak times vs. R-peak amplitudes
    figure;
    scatter(r_times_all_blocks, r_amplitudes, 'r', 'filled');
    title(sprintf('Block %d R-peak Locations', k));
    xlabel('Time (s)');
    ylabel('R-peak Amplitude');
end

%% 
HRV_metrics = struct('meanRR', [], 'SDNN', [], 'RMSSD', []);

for k = 1:num_blocks
    rr = all_rr{k};  
    % 1. Mean RR interval
    meanRR = mean(rr);
    % 2. Standard deviation of RR intervals (SDNN)
    SDRR = std(rr);
    % 3. Root mean square of successive differences (RMSSD)
    % measure of how much the herbeat is changing from one beat to anoter
    diff_rr = diff(rr);
    RMSSD = sqrt(mean(diff_rr.^2));

    % Store metrics
    HRV_metrics(k).meanRR = meanRR;
    HRV_metrics(k).SDNN = SDRR;
    HRV_metrics(k).RMSSD = RMSSD;

    %display metrics...als stored in HRV_metrics structure
    fprintf('Block %d:\n', k);
    fprintf('  Mean RR: %.4f s\n', meanRR);
    fprintf('  SDNN: %.4f s\n', SDRR);
    fprintf('  RMSSD: %.4f s\n\n', RMSSD);
end