inputdata = OpenBCIRAWdown;

%% Load data and filter BPF and notch at 50Hz 
LFvalue = 30;
HFvalue = 100;
LFBWF = string(LFvalue);
HFBWF = string(HFvalue);

sampleRate = 250;

% Extract data from the table
rawData1 = table2array(inputdata(:,2)); % Extract column 2
rawData2 = table2array(inputdata(:,3)); % Extract column 3
rawData3 = table2array(inputdata(:,4)); % Extract column 4

% Extract the valuable time series 
rawData1 = rawData1(1:end); 
rawData2 = rawData2(1:end);
rawData3 = rawData3(1:end);

% Apply bandpass and notch filter to each channel
filteredData1 = myBandPass(rawData1, LFvalue, HFvalue , 8, sampleRate);
filteredData1 = myNotch(filteredData1, sampleRate);

filteredData2 = myBandPass(rawData2, LFvalue, HFvalue, 8, sampleRate);
filteredData2 = myNotch(filteredData2, sampleRate);

filteredData3 = myBandPass(rawData3, LFvalue, HFvalue, 8, sampleRate);
filteredData3 = myNotch(filteredData3, sampleRate);

% Assign filtered data to WristCh variables
WristCh1 = filteredData1;
WristCh2 = filteredData2;
WristCh3 = filteredData3;

timePoints = (1:length(filteredData1)) / sampleRate;

%% Plot filtered data for each channel
figure
subplot(3,1,1)
plot(timePoints, WristCh1); grid on; title("WristCh1"); axis([0 max(timePoints) -150 150]); xlabel('Time (s)'); ylabel('Voltage (uV)');
subplot(3,1,2)
plot(timePoints, WristCh2); grid on; title("WristCh2"); axis([0 max(timePoints) -150 150]); xlabel('Time (s)'); ylabel('Voltage (uV)');
subplot(3,1,3)
plot(timePoints, WristCh3); grid on; title("WristCh3"); axis([0 max(timePoints) -150 150]); xlabel('Time (s)'); ylabel('Voltage (uV)');
sgtitle('notched only Wrist Channel Data');


PBfreq = 5; % pass band frequency
SBfreq = 10; % stop band frequency

%channel1
LPWristChN1 = abs(filteredData1);
LPWristChN1 = myLowPass(LPWristChN1,PBfreq,SBfreq,sampleRate);
LPWristChN1 = smoothdata(LPWristChN1,"movmedian",200);
x5denChN1 = wdenoise(filteredData1,3,DenoisingMethod="BlockJS"); %wavelet denoised signal
LPx5denChN1 = abs(x5denChN1);
LPx5denChN1 = myLowPass(LPx5denChN1,PBfreq,SBfreq,sampleRate);
LPx5denChN1 = smoothdata(LPx5denChN1,"movmedian",200);
difference1 = LPWristChN1 - LPx5denChN1;
%calculate snr
signal_power1 = mean(LPWristChN1.^2);
noise_power1 = mean(difference1.^2);
snr_dB1 = 10 * log10(signal_power1 / noise_power1);
disp(['SNR: ', num2str(snr_dB1), ' dB']);

%channel2
LPWristChN2 = abs(filteredData2);
LPWristChN2 = myLowPass(LPWristChN2, PBfreq, SBfreq, sampleRate);
LPWristChN2 = smoothdata(LPWristChN2, "movmedian", 200);
x5denChN2 = wdenoise(filteredData2, 3, DenoisingMethod="BlockJS"); % wavelet denoised signal level 5
LPx5denChN2 = abs(x5denChN2);
LPx5denChN2 = myLowPass(LPx5denChN2, PBfreq, SBfreq, sampleRate);
LPx5denChN2 = smoothdata(LPx5denChN2, "movmedian", 200);
difference2 = LPWristChN2 - LPx5denChN2;
signal_power2 = mean(LPWristChN2.^2);
noise_power2 = mean(difference2.^2);
snr_dB2 = 10 * log10(signal_power2 / noise_power2);
disp(['SNR: ', num2str(snr_dB2), ' dB']);

%channel3
LPWristChN3 = abs(filteredData3);
LPWristChN3 = myLowPass(LPWristChN3, PBfreq, SBfreq, sampleRate);
LPWristChN3 = smoothdata(LPWristChN3, "movmedian", 200);
x5denChN3 = wdenoise(filteredData3, 3, DenoisingMethod="BlockJS");
LPx5denChN3 = abs(x5denChN3);
LPx5denChN3 = myLowPass(LPx5denChN3, PBfreq, SBfreq, sampleRate);
LPx5denChN3 = smoothdata(LPx5denChN3, "movmedian", 200);
difference3 = LPWristChN3 - LPx5denChN3;
signal_power3 = mean(LPWristChN3.^2);
noise_power3 = mean(difference3.^2);
snr_dB3 = 10 * log10(signal_power3 / noise_power3);
disp(['SNR: ', num2str(snr_dB3), ' dB']);

%figure channel1
figure
subplot(2,1,1)
plot(timePoints, LPWristChN1,'g'); hold on;
plot(timePoints, LPx5denChN1,'b'); hold on;
plot(timePoints, difference1,'r'); grid on; 
title("WristCh1 enveloped"); 
axis([0 max(timePoints) -10 50]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('enveloped','denoised','difference');
text(0.1 * max(timePoints), 45, sprintf('SNR = %.2f dB', snr_dB1), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
subplot(2,1,2)
plot(timePoints, filteredData1,'b'); hold on;
plot(timePoints, x5denChN1,'r'); grid on; 
title("WristCh1 denoised"); 
axis([0 max(timePoints) -150 150]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('origiinal signal','denoised');

%figure channel2
figure;
subplot(2,1,1);
plot(timePoints, LPWristChN2, 'g'); hold on;
plot(timePoints, LPx5denChN2, 'b'); hold on;
plot(timePoints, difference2, 'r'); grid on; 
title("WristCh2 enveloped"); 
axis([0 max(timePoints) -10 50]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('enveloped', 'denoised', 'difference');
text(0.1 * max(timePoints), 45, sprintf('SNR = %.2f dB', snr_dB2), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
subplot(2,1,2);
plot(timePoints, filteredData2, 'b'); hold on;
plot(timePoints, x5denChN2, 'r'); grid on; 
title("WristCh2 denoised"); 
axis([0 max(timePoints) -150 150]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('origiinal signal','denoised');

%figure channel3 
figure;
subplot(2,1,1);
plot(timePoints, LPWristChN3, 'g'); hold on;
plot(timePoints, LPx5denChN3, 'b'); hold on;
plot(timePoints, difference3, 'r'); grid on; 
title("WristCh3 enveloped"); 
axis([0 max(timePoints) -10 50]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('enveloped', 'denoised', 'difference');
text(0.1 * max(timePoints), 115, sprintf('SNR = %.2f dB', snr_dB3), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
subplot(2,1,2);
plot(timePoints, filteredData3, 'b'); hold on;
plot(timePoints, x5denChN3, 'r'); grid on; 
title("WristCh3 denoised"); 
axis([0 max(timePoints) -150 150]); 
xlabel('Time (s)'); 
ylabel('Voltage (uV)');
legend('origiinal signal','denoised');



function filteredSignal = myLowPass(signal, freq1, freq2, Fs)
    d = designfilt('lowpassiir', ...        % Response type
           'PassbandFrequency', freq1, ...   % Frequency constraints
           'StopbandFrequency', freq2, ...
           'PassbandRipple', 4, ...          % Magnitude constraints
           'StopbandAttenuation', 55, ...
           'DesignMethod', 'butter', ...     % Design method
           'MatchExactly', 'passband', ...   % Design method options
           'SampleRate', Fs);                % Sample rate

    filteredSignal = filtfilt(d, signal);
end

function filteredSignal = myBandPass(signal, freq1, freq2, order, Fs)
    % A function that filters a given signal given the lower and upper
    % frequency limits of interest along with the the order of the filter and
    % the sampling frequency of the data 

    d = designfilt('bandpassiir', 'FilterOrder', order, ...
                   'HalfPowerFrequency1', freq1, 'HalfPowerFrequency2', freq2, ...
                   'DesignMethod', 'butter', 'SampleRate', Fs);

    filteredSignal = filtfilt(d, signal);
end

function filteredSignal = myNotch(signal, fs)
    % 50Hz notch filter 
    Fs = fs;
    Fn = Fs / 2;
    Wp = [48 52] / Fn;     % Stopband Frequency (Normalised)
    Ws = [47 53] / Fn;     % Passband Frequency (Normalised)
    Rp = 1;                % Passband Ripple
    Rs = 50;               % Passband Ripple (Attenuation)
    [n, Wp] = ellipord(Wp, Ws, Rp, Rs);          % Elliptic Order Calculation
    [z, p, k] = ellip(n, Rp, Rs, Wp, 'stop');    % Elliptic Filter Design: Zero-Pole-Gain 
    [sos, g] = zp2sos(z, p, k);                  % Second-Order Section For Stability

    filteredSignal = filtfilt(sos, g, signal);   % Filter Signal
end