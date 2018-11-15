%{
Created: November 8, 2018
Version: 1.0

Description: 
Extracts features from resp files from pilot data. Extracts the features
for all participants for all sub-movies. The baseline and commentary
periods are considered submovies. Can create visualizations of the
features.

Requirements:
Requires that only the respiration csv files are in the directory

Current Features Extracted:
1. Average Breathing Rate
2. Standard Deviation
3. Low Frequencies from Power Spectrum
4. High Frequencies from Power Spectrum
5. Band Energy Ratio (lf/hf)

Features in Progress of Being Extracted:
1. Breathing Period Estimate
2. Speed of Breathing (Fast/Slow)

Future Features to be Extracted:
1. Median Peak to Peak Time
2. Length of Inhale/Exhale
%}

clear
clc
close all

resp_files = dir;

%calibrations
data_Index = 6; %column containing resp data
time_ms_Index = 5;%column containing ms time data
subMovie_Index = 2;%column containing sub-movie title
totalSubMovies = 25;
windowTime_ms = 30000;
subWindowTime_ms = 15000;
minPeakDistance = 15; %samples; 1 sample is 1/Fs
prominenceScale = 0.05; % 5% of the mean value of the data
debug = true;       %debugging flag to review how well the peak method is doing
plotSubMovieData = true; %plot figures of feature windows

%load the data from the files
idx = 1;
testNum = 1;
while idx <= length(resp_files)
    if contains(resp_files(idx).name,'.csv')
        raw_data(testNum).movie =readtable(resp_files(idx).name, detectImportOptions(resp_files(idx).name));
        testNum = testNum + 1;
    end   
    idx = idx + 1;
end

for participant = 1 : length(raw_data)
    %extract the relevant data
    resp_data = table2array(raw_data(participant).movie(:,data_Index));
    time_data_ms = table2array(raw_data(participant).movie(:,time_ms_Index)) - table2array(raw_data(participant).movie(1,time_ms_Index));
    subMovieTitle = table2cell(raw_data(participant).movie(:,subMovie_Index));
    
    %% window the data for each sub movie & extract features
    curI = 1;
    prevI = curI;
    for subMovieI = 1:totalSubMovies
        %find the sub movie section of the data
        curI = find(strcmp(subMovieTitle, subMovieTitle(curI)),1,'last');
        subData = resp_data(prevI:curI);
        subTime_ms = time_data_ms(prevI:curI); % t1 t2 ... tn 
        Fs = 1000/(subTime_ms(2)-subTime_ms(1));
        curI = curI + 1;
        prevI = curI;
        %find the sub window
        numSampleInWindow = floor(windowTime_ms*Fs/1000);
        numSubSampleInWindow = floor(subWindowTime_ms*Fs/1000);
        %extract features on movie section
        w_start=1;
        windowI=1;
        while true
            %window data
            w_end = w_start+numSampleInWindow;
            if w_end > length(subTime_ms)
                w_end = length(subTime_ms);
            end
            data = subData(w_start:w_end);
            time = subTime_ms(w_start:w_end);
            %% extract features
            %find peaks
            [peakVals,locs,w,prom] = findpeaks(data,'MinPeakProminence',prominenceScale*mean(data), 'MinPeakDistance',minPeakDistance);
            timeBetweenPeaks_ms = time(locs(2:length(locs)))-time(locs(1:length(locs)-1));
            [z,mu,sigma] = zscore(timeBetweenPeaks_ms);
            p = normcdf(z);
            %find band energy ratio
            N = length(data);
            xdft = fft(data);
            xdft = xdft(1:floor(N/2)+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            freq = Fs*(0:(N/2))/N;
            % band energy ratio feature from DEAP
            logLowPsdx = log10(sum(psdx(freq<0.25 & freq > 0.05)));
            logHighPsdx = log10(sum(psdx(freq>=0.25 & freq < 5)));
            f_range = freq(freq>0.05 & freq<5);
            psdx_range = psdx(freq>0.05 & freq<5);
            [~,i] = max(psdx_range);
            %displaying resp data to identify if peaks properly captured
            if debug == true 
                figure(1);
                clf
                hold on
                plot(time/1000, data)
                plot(time(locs)/1000,prom,'rx')
                plot(time(locs)/1000, peakVals, 'o')
                title(strcat('Breathing Data for Window from Participant: ', {' '}, num2str(participant), ' and SubMovie: ', {' '}, (subMovieTitle(curI-1))), 'interpreter', 'none')
                xlabel('Time (s)')
                ylabel('Resp Data')
                legend('Raw Breathing Data', 'Peak Prominence', 'Peaks')
                str = {['Breathing Rate AVG: ',num2str(mu/1000),' (s)'],['Breathing Rate STD: ', num2str(sigma/1000),' (s)']};
                dim = [.6 .1 .1 .1];
                annotation('textbox',dim,'String',str,'FitBoxToText','on');
                %disp(timeBetweenPeaks_ms)
                
                figure(2);
                clf
                hold on
                plot(f_range,10*log10(psdx_range))
                plot(f_range(i), 10*log10(psdx_range(i)),'o')
                legend('Frequency Domain Data', 'Highest Power Frequency')
                str = {['Power Ratio: ',num2str(logLowPsdx - logHighPsdx),' (dB/Hz)'],['Breathing Period EST: ', num2str(1/f_range(i)),' (s)']};
                dim = [.6 .7 .1 .1];
                annotation('textbox',dim,'String',str,'FitBoxToText','on');
                grid on
                title(strcat('FFT of the Window Data for Participant: ',{' '}, num2str(participant), ' and SubMovie: ',{' '}, (subMovieTitle(curI-1))), 'interpreter', 'none')
                xlabel('Frequency (Hz)')
                ylabel('Power/Frequency (dB/Hz)')
                pause
            end
            %store data
            Results(participant).subMovie(subMovieI).window(windowI).p = p;
            Results(participant).subMovie(subMovieI).window(windowI).x_bar = mu;
            Results(participant).subMovie(subMovieI).window(windowI).s = sigma;
            Results(participant).subMovie(subMovieI).window(windowI).logLF = logLowPsdx;
            Results(participant).subMovie(subMovieI).window(windowI).logHF = logHighPsdx;
            Results(participant).subMovie(subMovieI).window(windowI).ratio = logLowPsdx - logHighPsdx;
            Results(participant).subMovie(subMovieI).window(windowI).freqPeak = f_range(i);           
            
            
            %move the index
            w_start = w_start+numSubSampleInWindow;
            windowI = windowI+1;
            %check for finish condition
            if w_start >= length(subTime_ms) || w_end == length(subTime_ms)
                break;
            end
        end
        
        %plot this subMovie data for the current participant
        if plotSubMovieData
            numWindows_movie = 1:length(Results(participant).subMovie(subMovieI).window);
            
            figure(3)
            clf
            hold on
            %plot(numWindows_movie,1000./[Results(participant).subMovie(subMovieI).window(numWindows_movie).freqPeak],'*-')
            plot(numWindows_movie,[Results(participant).subMovie(subMovieI).window(numWindows_movie).x_bar],'*-')
            plot(numWindows_movie,[Results(participant).subMovie(subMovieI).window(numWindows_movie).s],'*-')
            title([strcat('Features for Participant: ', {' '}, num2str(participant), ' and SubMovie: ', {' '},(subMovieTitle(curI-1)))], 'interpreter', 'none')
            xlabel('Windows')
            ylabel('Time (ms)')
            %legend('Breathing Period EST', 'Breathing Rate AVG', 'Breathing Rate STD')
            legend('Breathing Rate AVG', 'Breathing Rate STD')
            
            figure(4)
            clf
            hold on
            plot(numWindows_movie, [Results(participant).subMovie(subMovieI).window(numWindows_movie).ratio],'*-')
            plot(numWindows_movie, [Results(participant).subMovie(subMovieI).window(numWindows_movie).logLF],'*-')
            plot(numWindows_movie, [Results(participant).subMovie(subMovieI).window(numWindows_movie).logHF],'*-')
            title([strcat('Frequency Domain Features for Participant: ', {' '}, num2str(participant), ' and SubMovie: ', {' '}, (subMovieTitle(curI-1)))], 'interpreter', 'none')
            xlabel('Windows')
            ylabel('Power/Frequency (dB/Hz)')
            legend('LF/HF', 'Low Frequency', 'High Frequency')
            pause
        end
    end
end