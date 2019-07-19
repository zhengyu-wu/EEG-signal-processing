% Created by Zhengyu Wu 
% Last modified: 07/15/2019
% Pre-process dataset IIb from http://www.bbci.de/competition/ii/ 
% Dataset description: http://www.bbci.de/competition/ii/albany_desc/albany_desc_ii.pdf

Fs = 240;        % Sample rate (Hz)
dt = round(1000/Fs);    % Delta t (Ms)
nchans = 64;     % number of channels

%%%%%%%%%%%%%%%%%%%%%%% 
%%%%   Filtering    %%%
%%%%%%%%%%%%%%%%%%%%%%%

% % First, let's see the original data from one .mat file
% subplot(4,1,1);
% plot(signal);
% title('Orignal data')
% 
% % Draw a spectrum of original data
% Y = fft(signal(:,1));
% L = size(signal,1);
% P2 = abs(Y/L);
% P1 = P2(1:floor(L/2) + 1);
% P1(2:end-1) = 2*P1(2:end-1);
% fx = Fs*(0:(L/2))/L;
% subplot(4,1,2);
% plot(fx,P1)
% xlim([0 100])
% title('Single-Sided Amplitude Spectrum of original data')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')


% Filter 0.1-30.0
filter_order = 2;
hp = 0.1; % high pass frequency (Hz)
lp = 30;  % low pass frequency (Hz)
[b1,a1] = butter(2,lp/(Fs/2)); % lowpass filter
[b2,a2] = butter(2,hp/(Fs/2),'high'); % highpass filter

% Filt data
% signal_ret = filtfilt(b1, a1, signal);
% signal_ret = filtfilt(b2, a2, signal_ret);

% % Draw filtered data
% subplot(4,1,3);
% plot(signal_ret);
% title('Filtered data')
% 
% % Draw a spectrum of filtered data
% Y = fft(signal_ret(:,1));
% L = size(signal,1);
% P2 = abs(Y/L);
% P1 = P2(1:floor(L/2) + 1);
% P1(2:end-1) = 2*P1(2:end-1);
% fx = Fs*(0:(L/2))/L;
% subplot(4,1,4);
% plot(fx,P1)
% xlim([0 100])
% title('Single-Sided Amplitude Spectrum of filtered data')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%   Epoching     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Variables for epoch and baseline
ep_s = 0;   % epoch start 0ms
ep_e = 500; % epoch end 500ms after stim
bl_s = 0;   % baseline start 0ms
bl_e = 100; % baseline end 100ms after stim
p3_s = 250; % p3 start
p3_e = 450; % p3 end

% Index for times above
i_ep_s = round(ep_s / dt);
i_ep_e = round(ep_e / dt);
i_p3_s = round(p3_s / dt);
i_p3_e = round(p3_e / dt);
i_bl_s = round(bl_s / dt)+ 1;
i_bl_e = round(bl_e / dt)+ 1;

% Find when events happened
% zeroes_idx = find(StimulusCode == 0);
% events_idx = [];
% 
% for i = 1:length(zeroes_idx)
%     if i == length(zeroes_idx)
%         break
%     end
%     if zeroes_idx(i) == (zeroes_idx(i+1)-1)
%         continue;
%     else
%         events_idx = [events_idx, zeroes_idx(i)+1];
%     end 
% end

% Define sub-samling rate for export
k = 10;
dim = round((i_p3_e - i_p3_s) / k); % tot number of samples exported


% Sample_size * nchans * events_num
% epoch_df = zeros(abs(i_ep_s - i_ep_e), nchans, size(events_idx, 2));  % 0- 500ms
% data_df = zeros(dim, nchans, size(events_idx, 2));                    % 250-450ms


% Store epoched data
% for i = 1:size(events_idx, 2)  % 540
%     s_ix = events_idx(i) - i_ep_s;
%     e_ix = events_idx(i) + i_ep_e - 1;
%     epoch_df(:, :, i) = signal_ret(s_ix:e_ix, :);
%     
%     % baseline correct
%     baseline = mean(epoch_df(i_bl_s:i_bl_e, :, i), 1);
%     baseline = repmat(baseline, 125, 1);
%     epoch_df(:, :, i) = epoch_df(:, :, i) - baseline;
% 
%     % Save subsampled data
%     data_df(:, :, i) = epoch_df(i_p3_s:k:(i_p3_e - 1), :, i);
% end


% Divide epoch_df into target (StimulusType == 1) and non-target (StimulusType == 0)
% size(StimulusType(events_idx)) % 540
% target_idx = find(StimulusType(events_idx) == 1); % 90
% target_df = epoch_df(:,:,target_idx);
% nontarget_df = epoch_df;
% nontarget_df(:,:,target_idx) = [];

% Draw all events waves of channel 51 
% subplot(2,1,1);
% plot(squeeze(target_df(:,51,:)));
% subplot(2,1,2);

% Draw the average of all events of channel 51 
% plot(mean(target_df(:,51,:),3));
% hold on
% plot(mean(nontarget_df(:,51,:),3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%   Save processed data in .csv file and .npy file   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Params for saving persistent data
train_data = zeros(dim, nchans, 8000);  % 8000 is big enough to hold the sum of event nums from all files, trim later 
train_label = zeros(8000, 2);
test_data = zeros(dim, nchans, 8000);
test_label = zeros(8000, 1);
train_count = 0;
test_count  = 0;
train_break = [];
test_break = [];


% Iterate through .mat files
pfx = 'AAS0';

for s = 10:12
    for r = 1:8
        % Try loading data, if not possible, continue
        filename = [pfx num2str(s) 'R0' num2str(r)];
        try
            load([filename '.mat'])
        catch e
            continue;
        end
        
        % Filt
        signal_ret = filtfilt(b1, a1, signal);
        signal_ret = filtfilt(b2, a2, signal_ret);
        
        % Epoch
        zeroes_idx = find(StimulusCode == 0);
        events_idx = [];

        for i = 1:length(zeroes_idx)
            if i == length(zeroes_idx)
                break
            end
            if zeroes_idx(i) == (zeroes_idx(i+1)-1)
                continue;
            else
                events_idx = [events_idx, zeroes_idx(i)+1];
            end 
        end
        
        epoch_df = zeros(abs(i_ep_s - i_ep_e), nchans, size(events_idx, 2));  % 0- 500ms
        data_df = zeros(dim, nchans, size(events_idx, 2));                    % 250-450ms
        
        for i = 1:size(events_idx, 2)  % 540
            s_ix = events_idx(i) - i_ep_s;
            e_ix = events_idx(i) + i_ep_e - 1;
            epoch_df(:, :, i) = signal_ret(s_ix:e_ix, :);

            % baseline correct
            baseline = mean(epoch_df(i_bl_s:i_bl_e, :, i), 1);
            baseline = repmat(baseline, 125, 1);
            epoch_df(:, :, i) = epoch_df(:, :, i) - baseline;

            % Save subsampled data
            data_df(:, :, i) = epoch_df(i_p3_s:k:(i_p3_e - 1), :, i);
        end
        
        % If S12 save to test, otherwise train
        if s == 12
            labels = StimulusCode(events_idx);
            for q = 1:size(epoch_df, 3)
                test_count = test_count + 1;
                test_data(:, :, test_count) = data_df(:, :, q);
                test_label(test_count) = labels(q);
            end
            test_break = [test_break test_count-1]; % subtract one for python 0-indexing
        else
            labels = [StimulusType(events_idx) StimulusCode(events_idx)];
            for q = 1:size(epoch_df, 3)
                train_count = train_count + 1;
                train_data(:, :, train_count) = data_df(:, :, q);
                train_label(train_count, :) = labels(q, :);
            end
            train_break = [train_break train_count-1]; % subtract one for python 0-indexing
        end
            
    end
end

% Trim data
train_data = train_data(:, :, 1:train_count);
train_label = train_label(1:train_count, :);
test_data = test_data(:, :, 1:test_count);
test_label = test_label(1:test_count);

% Output path
out_dir = '/Users/wuzhengyu/Desktop/asignment2_data/';

% Export data for numpy (https://github.com/kwikteam/npy-matlab)
writeNPY(train_data, [out_dir 'train_data.npy']);
csvwrite([out_dir 'train_label.csv'], train_label);
writeNPY(test_data, [out_dir 'test_data.npy']);
csvwrite([out_dir 'test_label.csv'], test_label);
