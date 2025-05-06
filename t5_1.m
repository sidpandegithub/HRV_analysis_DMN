%% add paths to field trip and libeep(to read events)
addpath ('C:\Users\z5379472\OneDrive - UNSW\Desktop\fieldtrip-20250402')
%addpath('C:\Users\Admin\Documents\TOOLBOX\MarcusVollmer-HRV-58badf9')
addpath('C:\Users\z5379472\OneDrive - UNSW\Desktop\libeep-3.3.177-matlab');
addpath('C:\Users\z5379472\OneDrive - UNSW\Desktop\MarcusVollmer-HRV-58badf9');
ft_defaults;

%% Load data
%addpath('C:\Users\Admin\OneDrive - UNSW\Studies\TMS\TBD - TMS TEP Orientation\src\libeep-3.3.177-matlab');

DIR = 'C:\Users\z5379472\OneDrive - UNSW\Desktop\HRV_WIP\1058_ZHLI';

%load all 4 .cnt data files from 1 participant;
FILES = {'20250407_1058_ZHLI_S1.cnt'}; %; '20250407_1058_ZHLI_S2.cnt'; '20250407_1058_ZHLI_S3.cnt'; '20250407_1058_ZHLI_S4.cnt'};
all_hrv_tables = table();  % Will hold HRV measures from all files

for f = 1:length(FILES)
    FILE = FILES{f};
    fn = fullfile(DIR,FILE);
    info = eepv4_read_info(fn);

    %% take events and only keep 10,-10 (start and end of a 5 minute block)
    allLabels       = {info.triggers.label};
    keepIdx         = ismember(allLabels, {'10', '-10'});
    info.triggers   = info.triggers(keepIdx);

    %% set up config, only interested in BIP5 (ecg channel)
    cfg             = [];
    cfg.dataset     = fn;
    cfg.channel     = {'BIP5'};
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [0.5 40];
    cfg.bsfilter    = 'yes';
    cfg.bsfreq      = [49 51];
    cfg.demean      = 'yes';
    data            = ft_preprocessing(cfg);

    %% plot whole ecg signal vs. time
    ecg_signal  = data.trial{1};
    t           = data.time{1};
    fs          = data.fsample;
    % plot(t, ecg_signal);

    %% plot whole ecg signal with start and end blocks marked
    % figure;
    % plot(t, ecg_signal, 'b');
    % hold on;
    % title('Filtered ECG with 10/-10 Trigger Markers');
    % xlabel('Time (s)');
    % ylabel('ECG Amplitude');

    % Loop through triggers and plot markers
    % for i = 1:length(info.triggers)
    %     label = info.triggers(i).label;
    %     trigger_time = info.triggers(i).seconds_in_file;
    % 
    %     % Only plot for labels '10' and '-10'
    %     if strcmp(label, '10') || strcmp(label, '-10')
    %         % Get the sample index corresponding to this time
    %         trig_sample = round(trigger_time * fs);
    % 
    %         % Choose marker color based on label
    %         if strcmp(label, '10')
    %             col = 'g';  % green for '10'
    %         else
    %             col = 'm';  % magenta/pink for '-10'
    %         end
    % 
    %         % Plot the dot at the appropriate sample
    %         if trig_sample > 0 && trig_sample <= length(ecg_signal)
    %             plot(t(trig_sample), 0, 'o', ...
    %                 'Color', col, 'MarkerSize', 6, 'MarkerFaceColor', col); % ecg_signal(trig_sample)
    %         end
    %     end
    % end
    % 
    % legend('Filtered ECG', 'Trigger: 10', 'Trigger: -10');
    % axis tight

    %% segment data based on the blocks
    ecg_blocks  = {};
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

    % display duration of each block (should be 300s each: 5-minute blocks)
    for k = 1:length(time_blocks)
        duration = time_blocks{k}(end) - time_blocks{k}(1);
        fprintf('Block %d duration: %.2f seconds\n', k, duration);
    end

    %% plot all 6 in a subplot
    % num_blocks = length(ecg_blocks);
    % figure;
    % for i = 1:num_blocks
    %     subplot(num_blocks, 1, i);  % Rows = total blocks, 1 column
    %     plot(time_blocks{i} - time_blocks{i}(1), ecg_blocks{i}, 'b'); % time_blocks{i}
    %     title(sprintf('ECG Block %d', i));
    %     xlabel('Time (s)');
    %     ylabel('Amplitude (mV)');
    %     grid on;
    %     axis tight
    % end
    % sgtitle('All ECG Blocks');

    %% R peak detection for all blocks
    %amp threshold and time between each peak lines are from stackoverflow, not
    %sure how valid need to test more

    num_blocks = length(ecg_blocks);
    all_rr = cell(1, num_blocks);
    for k = 1:num_blocks
        signal  = [];
        t       = [];

        signal  = ecg_blocks{k};
        t       = time_blocks{k};

        % normalise signal
        signal      = detrend(signal);
        amp_thresh  = mean(signal) + 2.5*std(signal); % 0.35 * (max(signal) - min(signal));
        [~, locs]   = findpeaks(signal, 'MinPeakHeight', amp_thresh,'MinPeakProminence', 10);
        % [~, locs]   = findpeaks(signal, 'MinPeakHeight', amp_thresh,'MinPeakDistance', round(0.6 * fs));  % ~100 bpm max

        % Time of R-peaks
        r_times_all_blocks = t(locs);

        % RR intervals
        rr_intervals_all_blocks = diff(r_times_all_blocks);
        all_rr{k} = rr_intervals_all_blocks;

        figure;
        plot(t-t(1), signal);
        hold on;
        plot(t(locs)-t(1), signal(locs), 'ro');
        plot([0, t(end)-t(1)],[amp_thresh, amp_thresh],'linewidth',1,'linestyle',':','color','k')
        title(sprintf('Block %d R-peaks', k));
        xlabel('Time (s)');
        ylabel('Amplitude');
        axis tight
    end

    %% Checks
    num_rr_blocks   = length(all_rr);
    all_rr_outliers = cell(1, num_blocks);

    for r = 1:num_rr_blocks

        all_rr_outliers{r} = find(isoutlier(all_rr{r}));

        if any(all_rr_outliers{r})
            fprintf('Block %d outliers: %d \n', r, length(all_rr_outliers{r}));

            % Generate data for the plot
            signal  = [];
            t       = [];

            signal  = ecg_blocks{r};
            t       = time_blocks{r};

            % normalise signal
            signal      = detrend(signal);
            amp_thresh  = mean(signal) + 3*std(signal); % 0.35 * (max(signal) - min(signal));
            [~, locs]   = findpeaks(signal, 'MinPeakHeight', amp_thresh,'MinPeakDistance', round(0.6 * fs));  % ~100 bpm max

            % timing of outlier signals
            t_out       = t(locs([all_rr_outliers{r},all_rr_outliers{r}+1]));
            signal_out  = signal(locs([all_rr_outliers{r},all_rr_outliers{r}+1]));

            % Plot with outlier gap marked
            figure;
            plot(t-t(1), signal);
            hold on;
            scatter(t_out-t(1), signal_out, 'filled','ro');
            plot([0, t(end)-t(1)],[amp_thresh, amp_thresh],'linewidth',1,'linestyle',':','color','k')
            title(sprintf('Block %d R-peaks', k));
            xlabel('Time (s)');
            ylabel('Amplitude');
            axis tight
        end
    end


    %% In the situation where there are peaks that are not detected in a block
    %can mark peaks after zooming and panning
    %if there is an entire block with no peaks then just change the threshold
    %value in "amp_thresh=..." above

    num_blocks  = length(ecg_blocks);
    all_rr      = cell(1, num_blocks);

    for k = 1:num_blocks
        signal  = [];
        t       = [];

        signal  = ecg_blocks{k};
        t       = time_blocks{k};

        signal      = detrend(signal);
        amp_thresh  = mean(signal) + 3*std(signal); % 0.35 * (max(signal) - min(signal));
        [~, auto_locs] = findpeaks(signal, 'MinPeakHeight', amp_thresh,'MinPeakDistance', round(0.6 * fs));

        % Plot signal and auto peaks
        fig = figure;
        plot(t, signal); hold on;
        plot(t(auto_locs), signal(auto_locs), 'ro');
        title(sprintf('Block %d: Press "m" to mark, "q" to quit marking', k));
        xlabel('Time (s)'); ylabel('Amplitude');
        legend('ECG', 'Auto Peaks');
        axis tight
        zoom on; pan on;

        manual_locs = [];

        % Key press handling
        set(fig, 'KeyPressFcn', @(src, event) setappdata(src, 'key', event.Key));
        disp('Explore freely. Press "m" to mark peaks, "f" to finish/ go to the next figure.');
        while true
            waitfor(fig, 'CurrentCharacter');
            key = get(fig, 'CurrentCharacter');
            set(fig, 'CurrentCharacter', ' ');
            if strcmpi(key, 'm')
                disp('Click on missed peaks, then press Enter...');
                [x, ~] = ginput;
                for i = 1:length(x)
                    [~, idx] = min(abs(t - x(i)));
                    manual_locs = [manual_locs; idx];
                    plot(t(idx), signal(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
                end
            elseif strcmpi(key, 'q')
                break;
            end
        end

        % Combine and sort peaks
        all_locs                = sort(unique([auto_locs; manual_locs]));
        r_times_all_blocks      = t(all_locs);
        rr_intervals_all_blocks = diff(r_times_all_blocks);
        all_rr{k}               = rr_intervals_all_blocks;
    end

    % << STOP THE LOOP HERE AND SAVE THE OUTPUT >> %

    %%

    % << START A NEW LOOP, LOADING CLEANED RR INTERVALS FROM PRECEDING STEP
    % >>
    num_rr_blocks   = length(all_rr);

    hrv_table = table();

    for r = 1:num_rr_blocks
        RR = all_rr{r};

        % Check 1: RR must exist and be long enough
        if isempty(RR) || length(RR) < 10 || any(isnan(RR))
            fprintf('Skipping block %d in file %s: invalid or too short RR series.\n', r, FILE);
            continue;
        end

        % Filter RR series
        RR_filt = HRV.RRfilter(RR, 20);

        % Check 2: Filtered RR still valid
        if isempty(RR_filt) || length(RR_filt) < 10 || any(isnan(RR_filt))
            fprintf('Skipping block %d in file %s: filtered RR invalid.\n', r, FILE);
            continue;
        end

        % Try computation of HRV measures
        try
            RR_filt = HRV.RRfilter(all_rr{r},20);
            RMSSD   = HRV.RMSSD(RR_filt,0)*1000;
            HR      = HRV.HR(RR_filt,0);
            [~,~,LFHFratio,~,LF,HF] = HRV.fft_val_fun(RR_filt,1000); % [pLF,pHF,LFHFratio,VLF,LF,HF]

            % Store values in a table
            hrv_table.File{end+1,1}   = FILE;
            hrv_table.Block(end+1,1)  = r;
            hrv_table.HR(end+1,1)     = HR;
            hrv_table.RMSSD(end+1,1)  = RMSSD;
            hrv_table.LF(end+1,1)     = LF;
            hrv_table.HF(end+1,1)     = HF;
            hrv_table.RATIO(end+1,1)  = LFHFratio;

        catch ME
            fprintf('⚠️ HRV failed for block %d in %s: %s\n', r, FILE, ME.message);
            continue;
        end
    end

    % Add to master table
    all_hrv_tables = [all_hrv_tables; hrv_table];
end

%% Display as Excel sheet

disp(all_hrv_tables);
writetable(all_hrv_tables, fullfile(DIR, 'HRV_Results_All_Files.csv'));
%% EXTRAS
% rrHRV   = HRV.rrHRV(RR_filt,0);
% SDNN    = HRV.SDNN(RR_filt,0)*1000;
% pNN50   = HRV.pNN50(RR_filt,0)*100;
