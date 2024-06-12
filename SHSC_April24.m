%% SHSC-Study
% script to analyse the different conditions of stretch-hold shortening
% study * new April2024*
%--> filters rawdata with Butterworth 4. order
%--> calculates median of curves
%--> checks for different criteria
%--> calculates parameter

%% Part 1- sort conditions in struct with condition name
% Dependencies
% uipickfiles 

%clearvars
clear ; close all;

% load files
if ~exist('files','var')
    clear files
    files(:,1) = uipickfiles('FilterSpec','*.mat','Prompt','Select files to analyse');
elseif files{1} == 0
    clear files
    files(:,1) = uipickfiles('FilterSpec','*.mat','Prompt','Select files to analyse');
end

%define contraction conditions
iso = {'df', 'pf'};
dyn = {'sho', 'shod', 'str'};
ssc = {'ssc', 'sscd'};

conditions = [iso, dyn, ssc];

fsT = 2000; %sampling frequency Isomed
fsE = 5000; %sampling frequency EMG
sstime=11.1*fsT; %same timepoint for history dependent parameters after start of stim (1-s+6.6-s+3.5-s)
%%control on or off
controlofbasline = false;

%load data into different variables depending on contraction condition
for i = 1 : length(conditions) %outer loop through conditions
    x = 1; %set up variable to define trial number
    for j = 1 : length(files) %inner loop through all trials
        if strcmp(conditions{i},files{j}(63:end-5)) %compare condition with the end of the file name
            rawdata.(conditions{i})(x) = load(files{j}); %dynamic struct field by putting round brackets around the reference variable
            x = x +1; %increase trial count for condition by one
        end
    end
end

%% Part2 - Filters
  % Dependencies Wfilt Wcu

%f = f0 * (sqrt(2)-1)^(1/(2N)) see below for correcting f for multiple passes
fcT = 20; %filter cut-off for torque
fcA= 15;  %cut-off for angle

order = 4; %filter order
%Wn = 2*fc/fs/0.8022; %correct for dual-pass of 2nd order Butterworth
WnT = 2*fcT/fsT/0.8957; %correct for dual-pass of 4th order Butterworth--> Bakenecker et al.2019 Butterworth 4. order!!!
[bT,aT] = butter(order,WnT); %determine coefficients for Butterworth filter; if nothing else stated filter =lowpass filter!
WnA = 2*fcA/fsT/0.8957; %correct for dual-pass of 4th order Butterworth--> Bakenecker et al.2019 Butterworth 4. order!!!
[bA,aA] = butter(order,WnA); %determine coefficients for Butterworth filter; if nothing else stated filter =lowpass filter!

for i = 1 : length(conditions)
    for j = 1 : length(rawdata.(conditions{i}))
        tq = Wfilt(rawdata.(conditions{i})(j).Torque.values,20,1/rawdata.(conditions{i})(j).Torque.interval); % 10 Hz low-pass filt.
        tq = tq*-1;
        tqTime = rawdata.(conditions{i})(j).Torque.times;
        ang = Wfilt(rawdata.(conditions{i})(j).Angle.values,6,1/rawdata.(conditions{i})(j).Angle.interval); % 6 Hz low-pass filt.
        ang = interp1(rawdata.(conditions{i})(j).Angle.times,ang,tqTime,'linear','extrap');
        data.(conditions{i})(j).Torquefilt=filtfilt(bT,aT,tq);
        data.(conditions{i})(j).Anglefilt=filtfilt(bA,aA,ang);
        tqF(j,:)= filtfilt(bT,aT,tq);
        angF(j,:)= filtfilt(bA,aA,ang);
        
    end

    mCurves.(conditions{i})(:,1).Torquefiltmean=(mean(tqF));
    mCurves.(conditions{i})(:,1).Anglefiltmean=(mean(angF));

    medCurves.(conditions{i})(:,1).TqfMed=(median(tqF));
    medCurves.(conditions{i})(:,1).AfMed=(median(angF));
   clear var tqF angF
end

%% PART 3 - Correct for baseline torque at final angle
%if controlofbaseline=false without visual inspection otherwise every curve
%is plotted and time-window for baseline could be corrected

for i = 1 : length(conditions)

    %save baseline torque based on 1s-average before and after contraction
    results.(conditions{i}).BlTorqueafter = mean(medCurves.(conditions{i}).TqfMed(14*fsT:14.5*fsT));
    results.(conditions{i}).BlTorquebefore = mean(medCurves.(conditions{i}).TqfMed(0.5*fsT:1*fsT));

    % control if baseline is not accidently influenced by weird data, works only when running the whole script, turn off by setting control to 'false'
    if controlofbasline
        h = figure('units','normalized','position',[0,0,1,1]);
        plot(0.0005:0.0005:length(medCurves.(conditions{i}).TqfMed)/2000,medCurves.(conditions{i}).TqfMed);
        hold on
        plot(12.2:0.0005:13,medCurves.(conditions{i}).TqfMed(12.2*fsT:13*fsT))
        hold on
        plot(0.1:0.0005:1,medCurves.(conditions{i}).TqfMed(0.1*fsT:1*fsT))
        okay = questdlg('Accept baseline before?','Check','Yes','No','Yes');  % input dialog to accept or redo the baseline manually
        while strcmp(okay, 'No')
            [temp, ~] = ginput(1); %define beginning of new 1 second baseline
            results.(conditions{i}).BlTorquebefore = mean(medCurves.(conditions{i}).TqfMed(temp*2000:temp*2000+2000)); %overrite existing baseline
            plot(temp:0.0005:temp+1,medCurves.(conditions{i}).TqfMed(temp*2000:temp*2000+2000));
            okay = questdlg('Accept baseline before now?','Check','Yes','No','Yes');
        end
        okay = questdlg('Accept baseline after?','Check','Yes','No','Yes');  % input dialog to accept or redo the baseline manually
        while strcmp(okay, 'No')
            [temp, ~] = ginput(1); %define beginning of new 1 second baseline
            results.(conditions{i}).BlTorqueafter = mean(medCurves.(conditions{i}).TqfMed(temp*2000:temp*2000+2000)); %overrite existing baseline
            plot(temp:0.0005:temp+1,medCurves.(conditions{i}).TqfMed(temp*2000:temp*2000+2000));
            okay = questdlg('Accept baseline after now?','Check','Yes','No','Yes');
        end
        close(h)
    end
end

clearvars h temp
%% relevant parameters depending on condition
% DF:                   pFE (torque and Angle) (14-14.5s); rFEd (); rFE (sstime:sstime+0.5*fsT)
% STR: endSTR; meanSTR; pFE (torque and Angle) (14-14.5s); rFEd (); rFE (sstime:sstime+0.5*fsT)

% PF:   SS before SHO (50ms);                                                            meanSHO; rFD (sstime:sstime+0.5*fsT)
% SSC:                        endSTR; meanSTR; slope max (lastlocal min:sstime); endSHO; meanSHO; rFD (sstime:sstime+0.5*fsT)
% SSCd: SS before SHO (50ms); endSTR; meanSTR; slope max (lastlocal min:sstime); endSHO; meanSHO; rFD (sstime:sstime+0.5*fsT)
% SHO:  SS before SHO (50ms);                  slope max (lastlocal min:sstime); endSHO; meanSHO; rFD (sstime:sstime+0.5*fsT)
% SHOd: SS before SHO (50ms);                  slope max (lastlocal min:sstime); endSHO; meanSHO; rFD (sstime:sstime+0.5*fsT)


for i = 1 : length(conditions)

    if strcmp (conditions{i},'df')
        results.(conditions{i}).pFEt=(mean(medCurves.(conditions{i}).TqfMed(14*fsT:14.5*fsT)));
        results.(conditions{i}).pFEa=(mean(medCurves.(conditions{i}).AfMed(14*fsT:14.5*fsT)));
        results.(conditions{i}).SS=(mean(medCurves.(conditions{i}).TqfMed(sstime:sstime+0.5*fsT)))-results.(conditions{i}).BlTorqueafter;
        results.(conditions{i}).SSd=(mean(medCurves.(conditions{i}).TqfMed(6.5*fsT:7*fsT)));% steady-state before SHO of SSCd
    
    elseif strcmp (conditions{i},'str')
        results.(conditions{i}).pFEt=(mean(medCurves.(conditions{i}).TqfMed(14*fsT:14.5*fsT)));
        results.(conditions{i}).pFEa=(mean(medCurves.(conditions{i}).AfMed(14*fsT:14.5*fsT)));
        results.(conditions{i}).SS=(mean(medCurves.(conditions{i}).TqfMed(sstime:sstime+0.5*fsT)))-results.(conditions{i}).BlTorqueafter;
        results.(conditions{i}).SSd=(mean(medCurves.(conditions{i}).TqfMed(6.5*fsT:7*fsT)));% steady-state before SHO of SSCd

        %start & end of crank arm movement

        Derivative_a=diff(medCurves.(conditions{i}).AfMed(3*fsT:5*fsT)); %around 3.5 s rotation should appear as it is programmed at 3.5s
        meanbeforeSTR=mean(Derivative_a(1:800));
        stdbeforeSTR=std(Derivative_a(1:800));
        threshold=(meanbeforeSTR+stdbeforeSTR);
        xx=1;

        % search for all values which are higher than threshold during time
        % window
        for k=500:1500 %around 1000 rotation should appear

            if (Derivative_a(k)>threshold)
                idx(xx)=(k);
                xx=xx+1;
                % else
                %     idx(k)=nan;
            end
        end
        % check for 100ms =200 frames if values are increasing
        for kk=1:200 

        if Derivative_a(idx(1)+kk)>Derivative_a(idx(1)+kk-1)
            startSTR=idx(1)+5999; % to use startSTR for whole curves time 0-3s must be added
        else
            
            newStartidx=find (idx>idx(1)+kk);
            newstartSTR=idx(min(newStartidx))+5999;

        end
        end
        if startSTR==newstartSTR
            idxStart=startSTR;
        else
            idxStart=newstartSTR;
        end
      

               meanstartSTR=mean(medCurves.(conditions{i}).AfMed(5800:6800));
               stdstartSTR=std(medCurves.(conditions{i}).AfMed(5800:6800));
               [idx,idy]=find(medCurves.(conditions{i}).AfMed(5800:7400)>(meanstartSTR+stdstartSTR));
               mstartstr=min(idx);
               mstartstr=mstartstr+6800;

               %end
               meanendSTR=mean(data.(conditions{i})(j).Anglefilt(8600:9600));
               stdendSTR=std(data.(conditions{i})(j).Anglefilt(8600:9600));
               [idx,~]=find(data.(conditions{i})(j).Anglefilt(8400:10000)>(meanendSTR-stdendSTR));
               mendstr=min(idx);
               mendstr=mendstr+8400;
               results.(conditions{i})(j).mpeakendstr= data.(conditions{i})(j).Torquefilt(mendstr);
               results.(conditions{i})(j).mAMstr= trapz(rawdata.(conditions{i})(j).Torque.times(mstartstr:mendstr),data.(conditions{i})(j).Torquefilt(mstartstr:mendstr));
               results.(conditions{i})(j).mmeanSTR=mean(data.(conditions{i})(j).Torquefilt(mstartstr:mendstr));
               ind.(conditions{i})(j).start=mstartstr;
               ind.(conditions{i})(j).end=mendstr;

    elseif strcmp (conditions{i},'pf')
        results.(conditions{i})(j).ss500ms=(mean(data.(conditions{i})(j).Torquefilt(11.1*fsT:11.6*fsT)))-data.(conditions{i})(j).BlTorqueafter;



    elseif strcmp (conditions{i},'sho')
        results.(conditions{i})(j).ss500ms=(mean(data.(conditions{i})(j).Torquefilt(11.1*fsT:11.6*fsT)))-data.(conditions{i})(j).BlTorqueafter;

    elseif strcmp (conditions{i},'shod')
        results.(conditions{i})(j).ss500ms=(mean(data.(conditions{i})(j).Torquefilt(11.1*fsT:11.6*fsT)))-data.(conditions{i})(j).BlTorqueafter;

    elseif strcmp (conditions{i},'sscd')
        results.(conditions{i})(j).ss500ms=(mean(data.(conditions{i})(j).Torquefilt(11.1*fsT:11.6*fsT)))-data.(conditions{i})(j).BlTorqueafter;

    else
        results.(conditions{i})(j).ss500ms=(mean(data.(conditions{i})(j).Torquefilt(11.1*fsT:11.6*fsT)))-data.(conditions{i})(j).BlTorqueafter;
    end

end