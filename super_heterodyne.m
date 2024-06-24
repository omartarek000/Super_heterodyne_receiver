clear
clc
close all
% make the graph of the waveforms of time VERY IMPORTANT
%% reading the sounds 
folder_path = 'sounds/';
% creating an audio data store just like an array
ds = audioDatastore(folder_path);
% reading the files and getting their sample frequncey
% also creating  the arrays to hold the data to mainuplate them later
audiodata = cell(1,numel(ds.Files));
mono_data = cell(1,numel(ds.Files));
paded_audio_data = cell(1,numel(ds.Files));
modulated_signal=cell(1,numel(ds.Files));
sampling_freq = cell(1,numel(ds.Files));
new_samples = cell(1,numel(ds.Files));
% now reading the files 
for i = 1:6
    [audiodata{i},sampling_freq{i}] = read(ds);
end

% turning the wave from stereo aka two channels => mono one channel
for i = 1:numel(ds.Files)
    mono_data{i} = mean(audiodata{i},2);
end


    

% finding the wave with the longest points to make sure the others have the
% same length so we pad them with zeros
max_length = length(mono_data{i});
for i = 1:numel(ds.Files)
   if length(mono_data{i}) > max_length
       max_length = length(mono_data{i});
   end
end

% now we pad the waves 
for i =1:numel(ds.Files)        
    if length(mono_data) < max_length  % padding needed as it is not the longest wave
        paded_audio_data{i} = padarray(mono_data{i}, [max_length - length(mono_data{i}), 0], 'post');
    else
        paded_audio_data{i} = mono_data{i}; % No padding needed as it is the largest wave
    end
end

% plotting the waves in time domain
for i =1:numel(ds.Files)
    subplot(6,1,i)
    plot(paded_audio_data{i});
    title("this is the audio:"+i);
    xlabel("time");
    ylabel("magntuide");
end



%plotting the paded waves:
fs = sampling_freq{i}.SampleRate;
figure
for i= 1:numel(ds.Files)
    N=length(paded_audio_data{i});
    k=-N/2:(N/2)-1;  
    number = k*fs/N;
    signal_fft = fftshift(fft(paded_audio_data{i}));
    signal_AK = (1/N)* abs(signal_fft);
    subplot(3,2,i);
    plot((k*fs/N),signal_AK,'r')
    xlabel("frequncey(HZ))");
    ylabel("magnituide(db)");
    grid on
    title(sprintf("this is the wave number %i",i));
end

% now we want to modulate the signals and make it on a carrier freq
% but first we want to change the freq of the audio to be able to modulate
% it
% we need to make sure that the carrier frequncey whatever its values it
% should be less than the new sampling frequncey so that is why we made the
% sampling freq is 20 times the original sampling freq
%% AM modulator
figure;
for i= 1:6
    N = 740544*20;   % the number of samples increased as the the carrier frequncey increased
    K = -N/2:(N/2)-1;
    carrier_freq = 100000;
    delta_freq = 50000;
    fc = (i-1)* delta_freq + carrier_freq;   % making the formula of the carrier frequncey
    new_sample = interp(paded_audio_data{i},20);
     %disp(length(new_sample));
    % making the sampling frequncey up to 20 times the original 
    t = (0:length(new_sample)-1)./(44100*20);  % adjusting the time to the new sampling frequncey
    carrier = cos(2*pi*t*fc)/1;   % generating the carrier siganl 
    modulated_signal{i} = carrier'.*new_sample;  % assigning the modulated signal  
    fs_sample = 20*44100;    % adjusting the sampling frequncey
    x_axis = (K*fs_sample/N);    % adjusting the parameters to make the plot to make sure that the modulation works
    mod_signal=fftshift(fft(modulated_signal{i}));   % using fourier transform 
    signal_AK = (1/N)* abs(mod_signal);   % plotting all the graphs
    plot(x_axis,signal_AK);
    title("this the is the waveform after modulation");
    xlabel("frequncey");
    ylabel("magntuide");
    legend1 = {'audio1';'audio2';'audio3';'audio4';'audio5';'audio6'};
    legend(legend1);
    grid on;
    hold on ;
    
end

% putting all the waves in one waveform to send them 
sum_waveform = zeros(size(modulated_signal{1}));
for i = 1:6
    sum_waveform = modulated_signal{i} + sum_waveform;
end

%% now we enter the RF stage
% first like a radio we want you to enter your favouritre channel

while true
    fprintf("choose your desired channel\n");
    fprintf("1. Short_BBCArabic2\n");
    fprintf("2. Short_FM9090\n");
    fprintf("3. Short_QuranPalestine\n");
    fprintf("4. Short_RussianVoice\n");
    fprintf("5. also a russainvoice\n");
    fprintf("6. Short_SkyNewsArabia\n");
    fprintf("0.enter zero to exit the radio\n");
    channel = input("Enter the number of the channel from 1 to 6 (or 0 to exit): ");
    switch channel
        case 1
            signal = modulated_signal{1};
            index =1;
            disp("Selected channel 1.");
        case 2
            signal = modulated_signal{2};
            index =2;
            disp("Selected channel 2.");

        case 3
            signal = modulated_signal{3};
            index=3;
            disp("Selected channel 3.");

        case 4
            signal = modulated_signal{4};
            index=4;
            disp("Selected channel 4.");

        case 5
            signal = modulated_signal{5};
            index=5;
            disp("Selected channel 5.");

        case 6
            signal = modulated_signal{6};
            index=6;
            disp("Selected channel 6.");
        case 0
            disp("Exiting program.");
             % Exit the loop
             break;
        otherwise
            disp("Invalid channel number. Please enter a number from 0 to 6.");
    end
    % now as we selected our channel we want to make the bandpass filter

    BW = obw(signal,fs_sample); % getting the bandwidth of the selected wave
    fc_filter = 100000 + (index-1) * 50000;   % picking the carrier frequncey
    
    % the paramters of the filter
    A_stop1 = 60;		% Attenuation in the first stopband = 60 dB
    F_stop1 = fc_filter - 2*BW ;		% Edge of the stopband 
    F_pass1 = fc_filter - BW /2 ;	% Edge of the passband 
    F_pass2 = fc_filter + BW /2 ;	% Closing edge of the passband 
    F_stop2 = fc_filter + 2*BW ;	% Edge of the second stopband 
    A_stop2 = 60;		% Attenuation in the second stopband = 60 dB
    A_pass = 1;         % Amount of ripple allowed in the passband = 1 dB
    
    % making the bandpass filter
    
    band_pass_filter = fdesign.bandpass(F_stop1,F_pass1,F_pass2,F_stop2,A_stop1,A_pass,A_stop2,fs_sample);
    band_pass_filter = design(band_pass_filter,'equiripple');
    signal_RF = filter(band_pass_filter,sum_waveform);
    Signal_RF_fft = fftshift(fft(signal_RF));
    N=length(signal);
    k=-N/2:(N/2)-1;
    number = k*fs_sample/N ;
    figure;
    subplot(3,1,1)
    plot(number,abs(Signal_RF_fft),'g');
    title("this is the RF section of the selected channel:" +index);
    xlabel("frequncey (HZ)");
    ylabel("magntuide (dB)");
    grid on;
    
    
    %% now we enter the IF stage of the superheterdyne receiver 
    % and we want to make an oscillator to demodulate the signal 
    
    fc_demodulate = fc_filter + 25000;
    offset = 25000;   
      % we are adding the offset as the project pdf stated
    t = (0:length(signal)-1)./(44100*20);  % adjusting the time to the new sampling frequncey
    carrier = cos(2*pi*t*fc_demodulate);   % generating the carrier siganl
    signal_IF = carrier' .* signal_RF;
    Signal_IF = fftshift(fft(signal_IF));
    subplot(3,1,2)
    plot(number,abs(Signal_IF),'r');
    title("the signal" +index+ " after converting it to IF frequncey");
    xlabel("frequncey (HZ)");
    ylabel("magntuide (dB)");
    grid on;
    
    % now we complete the IF stage by implmenting the bandpassfilter
    
    A_stop1 = 60;		% Attenuation in the first stopband = 60 dB
    F_stop1 =  offset - BW*1.5 ;		% Edge of the stopband 
    F_pass1 = offset - BW/1 ;	% Edge of the passband 
    F_pass2 = offset + BW /1 ;	% Closing edge of the passband 
    F_stop2 = offset + BW*1.5 ;	% Edge of the second stopband 
    A_stop2 = 60;		% Attenuation in the second stopband = 60 dB
    A_pass = 1;         % Amount of ripple allowed in the passband = 1 dB
    
    band_pass_filter = fdesign.bandpass(F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, A_stop2, 20 * 44100);

    band_pass_filter = design(band_pass_filter,'equiripple');
    
    filter_signal_IF = filter(band_pass_filter,signal_IF);
    Filter_signal_IF = fftshift(fft(filter_signal_IF));
    
    subplot(3,1,3)
    plot(number,abs(Filter_signal_IF),'b');
    title("the IF signal after the bandpass filter");
    xlabel("frequncey(HZ)");
    ylabel("magntuide(dB)");
    grid on ;
    
    %% now we enter the baseband Stage
    % first thing we intiallize the carrier signal to demodulate the
    % remaning signal 
    base_band_freq = 25000;
    carrier = cos(2*pi*t*base_band_freq);
    signal_BB = filter_signal_IF .* carrier';
    figure;
    subplot(2,1,1);
    plot(number,abs(fftshift(fft(signal_BB))));
    title("signal after we return it to the baseband");
    xlabel("frequncey(HZ)");
    ylabel("magntuide(dB)");
    grid on ;

    
    
    %  now we implment the low_pass_filter to make to the baseband
    %  frequncey range
    F_stop = BW / 2 + 30000;                            % Edge of the stopband
    F_pass = BW/2 ;                                        % Edge of the passband
    A_stop = 60;                                        % Attenuation in the second stopband = 60 dB
    A_pass = 1;                                         % Amount of ripple allowed in the passband = 1 dB

    low_pass_filter = fdesign.lowpass(F_pass, F_stop, A_pass, A_stop, 20 * 44100);
    low_pass_filter = design(low_pass_filter, 'butter');
    final_signal = filter(low_pass_filter,signal_BB);
    subplot(2,1,2);
    plot(number,abs(fftshift(fft(final_signal))));
    title("signal after we return it to the baseband");
    xlabel("frequncey(HZ)");
    ylabel("magntuide(dB)");
    grid on ;
    
    %% now the final part we return the signal to its original sampling frequncey and we sound it to see how it performed
    final_signal = 10*resample(final_signal,1,20);
    sound(final_signal,44100);
    
    
end


    



    


