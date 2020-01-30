function GetTraces_BWPatches_new_ROIs_DL(yearmonthday, mouse, experiment, adata_dir, rdata_dir, bscope2, Matched, Plevels)

% Extracts and saves dF traces and tuning curves of each trial for each 
% patch of the sparse noise stimulation respectively for a single piezo plane.
%
% Input:
% - yearmonthday:           Date of the experiment [cell]
% - mouse:                  Name of the animal [string]
% - experiment:             Experiment #
% - adata_dir:              Directory of the analyzed data [string]
% - rdata_dir:              Directory of the raw data [string]
% - bscope2:                Logical whether data was acquired on BScope2
% - Matched:                Logical whether ROIs are grouped to the 
%                           respective moving grating stimulation experiment
% - Plevels:                Number of Piezo levels
%
% David Laubender


% parameter switchboard
Smoothing = 3;
StimOnsetDelay = 0.1;
StimOvershoot = 0.2;

% define path and load stimulus settings
data_path = [char(rdata_dir) '\' char(mouse) '\ImagingData\' char(yearmonthday) '\'];
cd(data_path);
[aux_data ids ids_new frame_times stimarray] = getauxstim_Ret_DL(experiment,cd, 2, 4, bscope2);
StimSettings = stimarray.StimSettings;

%% Load data
for il = 1:length(experiment)
    rdata_path = [char(rdata_dir) '\' char(mouse) '\ImagingData\' char(yearmonthday) '\'];
    cd(rdata_path)
    
    info_dir = [char(rdata_dir) '\' char(mouse) '\ImagingData\' char(yearmonthday) '\'];
    cd(info_dir)
    info.FileID = experiment;
    files =  dir(['exp' num2str(info.FileID) '*.tif']);
    a= imfinfo(files(1).name);
    info = a(1);
    info.level = Plevels;

    adata_path = [char(adata_dir) '\' char(mouse) '\' char(yearmonthday) '\'];
    cd(adata_path)
    ana_load = dir([char(mouse) '-Adata-' num2str(experiment) '.mat']);
    %ana_load = dir(['*' experiment{1} '*.mat']);
    ana_load = ana_load.name;
    ana_data{il} = load(ana_load, '-mat');
end

ROIs = ana_data{1}.ROIs;

%% Calculate analysis settings
eyes = size(ids,2);

% Define output variables
Black_Traces_contra=[]; Black_Avr_tst_contra=[]; TC_1D_Black_contra=[]; B_test_contra=[];  B_test_trace_contra=[];
White_Traces_contra=[]; White_Avr_tst_contra=[]; TC_1D_White_contra=[];W_test_contra=[]; W_test_trace_contra=[];
Black_Traces_ipsi=[]; Black_Avr_tst_ipsi=[]; TC_1D_Black_ipsi=[]; B_test_ipsi=[];  B_test_trace_ipsi=[];
White_Traces_ipsi=[]; White_Avr_tst_ipsi=[]; TC_1D_White_ipsi=[];W_test_ipsi=[]; W_test_trace_ipsi=[];
        
for eyerun = 1:eyes
    StimOnsetFrames = ids(eyerun).StimOnsetFrames;
    StimOffsetFrames = ids(eyerun).StimOffsetFrames;
    
    SamplingFreq = regexp(info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
    SamplingFreq = str2num(SamplingFreq{1}) / info.level;
    
    NumTrials = length(StimOnsetFrames{1});
    NumPatches = length(StimSettings.Stimulus);
    NumROIs = size(ROIs,2);
    
    BaselineWindow = round((-0.7*StimSettings.ITIlength)*SamplingFreq):round(0.5*StimOnsetDelay*SamplingFreq);
    TestingWindow = round(StimOnsetDelay*SamplingFreq):round((StimSettings.StimDur+StimOvershoot)*SamplingFreq);
    DisplayWindow = round((-0.7*StimSettings.ITIlength)*SamplingFreq):round((StimSettings.StimDur+(15*StimSettings.ITIlength))*SamplingFreq);
    
    disp(' ');
    disp(['Baseline period  = ' sprintf('%5.2f',min(BaselineWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(BaselineWindow)./SamplingFreq) ...
        's (frames: ' num2str(min(BaselineWindow)) ' .. ' num2str(max(BaselineWindow)) ')']);
    disp(' ');
    disp(['Display period   = ' sprintf('%5.2f',min(DisplayWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(DisplayWindow)./SamplingFreq) ...
        's (frames: ' num2str(min(DisplayWindow)) ' .. ' num2str(max(DisplayWindow)) ')']);
    disp(['Testing period   = ' sprintf('%5.2f',min(TestingWindow)./SamplingFreq) ' to ' sprintf('%5.2f',max(TestingWindow)./SamplingFreq) ...
        's (frames: ' num2str(min(TestingWindow)) ' .. ' num2str(max(TestingWindow)) ')']);
    disp(' ');
    
    
    %% Loop ROI's ang get PSTH's and tuning curves  
    for nr = 1:NumROIs
        
        % define F
        raw_data = ROIs(nr).activity(logical(ids(eyerun).eye_open));
        F = raw_data;
        F = smooth(F,Smoothing);
        F(isnan(raw_data))= 1;
        B_trace_ipsi = [];
        W_trace_ipsi = [];
        B_trace_contra = [];
        W_trace_contra = [];

        % Get tuning curves
        for P = 1:NumPatches
            for t = 1:NumTrials
                Y = StimSettings.Stimulus(P).Y;
                X = StimSettings.Stimulus(P).X;    
                C = StimSettings.Stimulus(P).Color;
                
                if length(StimOnsetFrames{P}) < t || length(StimOffsetFrames{P}) < t
                    disp(['Warning: Patch (x=' num2str(X) ', y=' num2str(Y) ', trial ' num2str(t) ') was not in recording..']);
                    if C == 0
                        if eyerun == 1
                            Black_Traces_ipsi{nr,Y,X}(t,:) = NaN;
                            Black_Avr_tst_ipsi(nr,Y,X,t) = NaN;
                            TC_1D_Black_ipsi(nr, ((Y-1)*NumTrials)+t, X) =NaN;
                        elseif eyerun == 2
                            Black_Traces_contra{nr,Y,X}(t,:) = NaN;
                            Black_Avr_tst_contra(nr,Y,X,t) = NaN;
                            TC_1D_Black_contra(nr, ((Y-1)*NumTrials)+t, X) =NaN;
                        end
                    else
                        if eyerun == 1
                            White_Traces_ipsi{nr,Y,X}(t,:) = NaN;
                            White_Avr_tst_ipsi(nr,Y,X,t) = NaN;
                            TC_1D_White_ipsi(nr, ((Y-1)*NumTrials)+t, X) =NaN;
                        elseif eyerun == 2
                            White_Traces_contra{nr,Y,X}(t,:) = NaN;
                            White_Avr_tst_contra(nr,Y,X,t) = NaN;
                            TC_1D_White_contra(nr, ((Y-1)*NumTrials)+t, X) =NaN;
                        end
                    end
                else
                    DISPix = StimOnsetFrames{P}(t)+DisplayWindow;
%                     if max(DISPix) > NumFrames
%                         if nr == 1
%                             disp(['Warning: Restricted test-index for (patch x=' num2str(X) ', y=' num2str(Y) ', trial ' num2str(t) ')']);
%                         end
%                         DISPix(DISPix>NumFrames) = NumFrames;
%                     end
                     TESTix = StimOnsetFrames{P}(t)+TestingWindow;
%                     if max(TESTix) > NumFrames
%                         if nr == 1
%                             disp(['Warning: Restricted test-index for (patch x=' num2str(X) ', y=' num2str(Y) ', trial ' num2str(t) ')']);
%                         end
%                         TESTix(TESTix>NumFrames) = NumFrames;
%                     end
                    if C == 0
                        if eyerun == 1
                            Black_Traces_ipsi{nr,Y,X}(t,:) =F(DISPix );
                            temp=F(TESTix);
                            B_test_ipsi(nr,Y,X,t,:)=temp;
                            Black_Avr_tst_ipsi(nr,Y,X,t)=nanmean(temp);
                            temp2=F(DISPix);
                            TC_1D_Black_ipsi(nr, ((Y-1)*NumTrials)+t, X) =nanmean(temp2); %calculated over the whole display window to calculate significance of response to patches over time
                        elseif eyerun == 2
                            Black_Traces_contra{nr,Y,X}(t,:) =F(DISPix );
                            temp=F(TESTix);
                            B_test_contra(nr,Y,X,t,:)=temp;
                            Black_Avr_tst_contra(nr,Y,X,t)=nanmean(temp);
                            temp2=F(DISPix);
                            TC_1D_Black_contra(nr, ((Y-1)*NumTrials)+t, X) =nanmean(temp2);
                        end
                    else
                        if eyerun == 1
                            White_Traces_ipsi{nr,Y,X}(t,:) =F(DISPix );
                            temp=F(TESTix);
                            W_test_ipsi(nr,Y,X,t,:)=temp;
                            White_Avr_tst_ipsi(nr,Y,X,t)=nanmean(temp);
                            temp2=F(DISPix);
                            TC_1D_White_ipsi(nr, ((Y-1)*NumTrials)+t, X) = nanmean(temp2);
                        elseif eyerun == 2
                            White_Traces_contra{nr,Y,X}(t,:) =F(DISPix );
                            temp=F(TESTix);
                            W_test_contra(nr,Y,X,t,:)=temp;
                            White_Avr_tst_contra(nr,Y,X,t)=nanmean(temp);
                            temp2=F(DISPix);
                            TC_1D_White_contra(nr, ((Y-1)*NumTrials)+t, X) = nanmean(temp2);
                        end
                    end  
                end
            end
            
            if C==0
                if eyerun == 1
                    temp3=(squeeze(nanmean(B_test_ipsi(nr,Y,X,:,:),4)))';
                    B_trace_ipsi=[B_trace_ipsi temp3];
                elseif eyerun == 2
                    temp3=(squeeze(nanmean(B_test_contra(nr,Y,X,:,:),4)))';
                    B_trace_contra=[B_trace_contra temp3];
                end
            else
                if eyerun == 1
                    temp4=(squeeze(nanmean(W_test_ipsi(nr,Y,X,:,:),4)))';
                    W_trace_ipsi=[W_trace_ipsi temp4];
                elseif eyerun == 2
                    temp4=(squeeze(nanmean(W_test_contra(nr,Y,X,:,:),4)))';
                    W_trace_contra=[W_trace_contra temp4];
                end
            end
        end
        
        if eyerun == 1
            B_test_trace_ipsi(nr,:)=B_trace_ipsi;
            W_test_trace_ipsi(nr,:)=W_trace_ipsi;
        elseif eyerun == 2
            B_test_trace_contra(nr,:)=B_trace_contra;
            W_test_trace_contra(nr,:)=W_trace_contra;
        end
    end
end
Average_B_ipsi= median(B_test_trace_ipsi,1);
Average_W_ipsi = median(W_test_trace_ipsi,1);
Average_B_contra= median(B_test_trace_contra,1);
Average_W_contra =  median(W_test_trace_contra,1);
close all

save(['exp' num2str(experiment) 'TestTraces_newROIs_matched.mat'],'B_test_trace_ipsi','W_test_trace_ipsi','Average_B_ipsi','Average_W_ipsi','B_test_trace_contra','W_test_trace_contra','Average_B_contra','Average_W_contra');

save(['exp' num2str(experiment) 'TuningCurves_new_ROIs_matched.mat'], 'White_Traces_ipsi','Black_Traces_ipsi','TC_1D_Black_ipsi','TC_1D_White_ipsi','White_Avr_tst_ipsi','Black_Avr_tst_ipsi',...
    'White_Traces_contra','Black_Traces_contra','TC_1D_Black_contra','TC_1D_White_contra','White_Avr_tst_contra','Black_Avr_tst_contra',...
    'SamplingFreq', 'DisplayWindow', 'TestingWindow');


    
    
    
    
    
    
    