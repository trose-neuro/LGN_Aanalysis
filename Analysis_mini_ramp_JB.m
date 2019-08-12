%Analysis of minis/failures and ramps for dLGN experiments (created by SW181211)



%CHECK FLAG SECTION (line 14-21) AS WELL AS DIRECTORY SECTION (line 35-39)
%NESTED FUNCTIONS:
%parseExperimentsXls_dLGN
%rampanalysis
%minianalyis,
%dLGN_plot_analysis

%%INITIATION

% clear all;%delete all current variables in workspace
% close all;%close all open windows/figures

%%%%%%IMPORTANT FLAGS, PLEASE CHANGE HERE%%%%%%%%%%%%%%%%
filterephys=1;%TR2019 filter ephys traces (filter specs see nested functions)
analyze_mini=1;%flag if either mini only and/or ramp should be analyzed (1 or 0)
analyze_ramp=1;
analyze_spon=0;
morpho=0;%run morphology script
scracm=0;%run 2scracm script/analysis; REMEMBER TO SET RIGHT FLAGS
fanalysis=0;
factor=4;%std threshold factor
display=0;%flag to display plot (1 or 0)
ramp_rtrace=1;%save raw ephystraces or not (1 or 0)
savefile=1;%save file at the end or not
scracm_ana=0;%run 2scracm further analysis; ODI etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if analyze_mini==1 || analyze_ramp==1 || morpho==1 || scracm==1
    disp('dLGN Analysis');
    
    sent = 'Which user data will be analyzed? type in 0 for SW or 1 for MF\n';%text appears in command window
    user = input(sent);%waiting for input which is either 0 or 1
    
    experimentator = 'SW';%default SW data
    if user==1
        experimentator = 'MF';%MF data
    end
    
    dLGN_ephys={};%empty structure for saving variables
    
    if ispc
        %     %%%%%%DIRECTORIES%%%%%%%
        rdata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\RawData';%data directory of raw data;change accordingly
        %     adata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';%data directory of extracted date;change accordingly
        adata_dir         = 'D:\LGN project';%data directory of extracted date;change accordingly
        ExpXls            = 'I:\Martin Fernholz\LGN _Project_Common\dLGN_ephys_analysis_excel_spread_sheet\Experiments_dLGN_2019_testing.xlsx';%directory where excel batch file is located;change accordingly
        %     %%%%%%%%%%%%%%%%%%%%%%%%
    else
        
        %%%%DIRECTORIES - TR2019 %%%%%%%
        %     % TR2019: MAC: mount smb shares to /Volumes/first (easiest: do command+K in
        %     % finder)
        %     % 'smb://10G.ISI01.neuro.mpg.de/archive_bonhoeffer_group$/Martin Fernholz/'
        %     % 'smb://S15.neuro.mpg.de/R-bonhoe/Share/Simon/LGN_2019_SW_MF_JB_TR/dLGN_ephys_analysis_excel spread sheet/'
        rdata_dir         = '/Volumes/Martin Fernholz/LGN _Project_Common/RawData/';
        adata_dir         = '~/Analysis/dLGN_ephys_Analysis/';
        ExpXls            = '/Volumes/Martin Fernholz/LGN _Project_Common/dLGN_ephys_analysis_excel_spread_sheet/Experiments_dLGN_2019_testing.xlsx';
        %%%%%%%%%%%%%%%%%%%%%%
    end
    %% parse Experiments XLS database
    batchopt          = parseExperimentsXls_dLGN(ExpXls,user);%calls the nested function parseExperimentsXls_dLGN and considers the user flag (1 or 0)
    nummice           = length(batchopt.mouse);%length of experiments to be analyzed
    %%
    
    
    adder=1;%counting variable
    for i=1:nummice%for loop over experiments across days
        datapath=fullfile(rdata_dir, batchopt.mouse{i}, filesep);%directory and name of experiments (from excel sheet)
        cd(char(datapath));%go to directory
        
        for k=1:length(batchopt.exp_ids{i})%loop in bigger loop for each cell per experimental day
            if batchopt.exp_ids{i}(k)<10%for cells with id less then XX0010, e.g., XX0001-XX0009
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));
            else
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));%for cells with id mor then XX0010, e.g., XX0011-XX0099
            end
            fold_name=[experimentator n_str];%complete cell folder name such as SW0001 or MF0001
            exp_folder=fullfile(datapath,fold_name);%complete folder and directory
            list=dir([char(exp_folder) filesep '*.xsg']);%xsg files per cell
            len=length(list);%number of xsg files per cell
            for j=1:len
                load([char(exp_folder) filesep list(j).name],'-mat');%load each xsg file
                iterations(:,j)=header.loopGui.loopGui.iterations;%find out whether mini or ramp recording
                clamp(:,j)=mean(data.ephys.trace_1(1:100,:))<100;
                
            end
            
            ramp=find(iterations==11);%ramp recordings
            failure1=find(iterations==50);%mini recordings
            failure2=find(iterations==100);%mini recordings
            spon=find(iterations==5);
            red_rep=find(iterations==3);
            mapping=find(iterations==64);
            
            
            
            disp(['CURRRENT EXPERIMENT is ', char(batchopt.mouse{i}), fold_name]);
            if user==0
                disp([num2str(length(ramp)/11),' ramp recordings']);
                disp([num2str(length(red_rep)/3),' red rep recordings']);
            else
                disp([num2str(length(ramp)),' ramp recordings']);
                disp([num2str(length(red_rep)),' red rep recordings']);
            end
            disp([num2str(length(failure1)),' failure recordings with 50 reps']);
            disp([num2str(length(failure2)),' failure recordings with 100 reps']);
            disp([num2str(length(spon)),' spontaneous recordings']);
            iterations=[];
            
            %% RAMP ANALYSIS
            if analyze_ramp==1
                [blue_ramp, red_ramp ODI Ratio_AMPA_NMDA]=rampanalysis_JB(list, ramp, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir,batchopt.injection_order{i}(k));%use nested function rampanalysis
                ramp=[];%clear variables for next iteration
                %list=[];%clear variables for next iteration
            end
            
            %% MINI ANALYSIS
            %             if analyze_mini==1
            %                 if length(failure1)>=1 & length(failure2)>=1% Why do you run it twice!?!
            %                     [Ipsi, Contra]=minianalysis_JB(list, failure1, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir, batchopt.injection_order{i}(k));%call minianalysis
            %                     [Ipsi, Contra]=minianalysis_JB(list, failure2, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir, batchopt.injection_order{i}(k));%call minianalysis
            %                     failure1=[];
            %                     failure2=[];
            %                     %list=[];
            %                 elseif length(failure1)>=1 & length(failure2)==0%
            %                     [Ipsi, Contra]=minianalysis_JB(list, failure1, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir, batchopt.injection_order{i}(k));
            %                     failure1=[];
            %                     %list=[];
            %                 elseif length(failure2)>=1 & length(failure1)==0%
            %                     [Ipsi, Contra]=minianalysis_JB(list, failure2, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir, batchopt.injection_order{i}(k));
            %                     failure2=[];
            %                 else
            %                     disp('No failure recording');
            %                     Ipsi=[];
            %                     Contra=[];
            %                 end
            %             end
            
            % extraction version by JB
            if analyze_mini
                if length(failure1)
                    [red_failure_AMPA, blue_failure_AMPA, red_failure_NMDA, blue_failure_NMDA] = ...
                        minianalysis_JB_new(list, failure1, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys);
                    failure1=[];
                    failure2=[];
                elseif length(failure2)
                    %                     [red_failure_AMPA, blue_failure_AMPA, red_failure_NMDA, blue_failure_NMDA] = ...
                    %                         minianalysis_JB_new(list, failure2, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys);
                    failure1=[];
                    failure2=[];
                    red_failure_AMPA = [];
                    blue_failure_AMPA = [];
                    red_failure_NMDA = [];
                    blue_failure_NMDA = [];
                    
                else
                    disp('No failure recording');
                    %                     Ipsi=[];
                    %                     Contra=[];
                    red_failure_AMPA = [];
                    blue_failure_AMPA = [];
                    red_failure_NMDA = [];
                    blue_failure_NMDA = [];
                    
                end
            end
            
            
            
            %% %% Morphology extraction and analysis
            if morpho==1
                morpho_folder=[char(exp_folder) filesep 'SingleTraces' char(batchopt.mouse{i}), fold_name];
                morpho_data=dir([char(exp_folder) filesep 'SingleTraces' char(batchopt.mouse{i}), fold_name filesep '*.swc']);
                if isempty(morpho_data)==0;
                    morphology=morpho_LGN(morpho_folder,morpho_data,batchopt.morphox{i}(k),batchopt.morphoy{i}(k),batchopt.morphoz{i}(k));
                else
                    disp('No morphology');
                    morphology=[];
                end
            end
            
            %% 2SCRACM
            if scracm==1
                [blue_scracm, red_scracm]=LGN_scracm(list, mapping, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir,0);%use nested function rampanalysis
            end
            
            
            clamp=[];
            %%  prepare structure for all cells
            
            %             if batchopt.slice{i}(k)==1 %%% WHAT THE FUCK! what if slice=3
            %
            %                 slice_nr = 1 ;
            %             else
            %                 slice_nr = 2 ;
            %             end
            %
            %                 dLGN_ephys.data{1,1}='Animal ID';
            %                 dLGN_ephys.data{1,2}='Experimental ID';
            %                 dLGN_ephys.data{1,3}='Slice';
            %                 dLGN_ephys.data{1,4}='Injection order';
            %                 dLGN_ephys.data{1,5}='Category';
            %                 dLGN_ephys.data{1,6}='MDside';
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 dLGN_ephys.data{1,7}='Ramp blue';
            %                 dLGN_ephys.data{1,8}='Ramp red';
            %                 dLGN_ephys.data{1,9}='Failure Ipsi';
            %                 dLGN_ephys.data{1,10}='Failure Contra';
            %                 dLGN_ephys.data{1,11}='Morphology';
            %                 dLGN_ephys.data{1,12}='Scracm Blue';
            %                 dLGN_ephys.data{1,13}='Scracm Red';
            %                 dLGN_ephys.data{1,14}='ODI_AMPA';
            %                 dLGN_ephys.data{1,15}='ODI_NMDA';
            %                 dLGN_ephys.data{1,16}='Ratio_A_N_r';
            %                 dLGN_ephys.data{1,17}='Ratio_A_N_b';
            %                 dLGN_ephys.data{1,18}='Morph DOi';
            %%%%%     DATA          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %dLGN_ephys.data{adder+1,1}=[char(batchopt.mouseID{i})];
            LGN(adder).animal_name=[char(batchopt.mouseID{i})];
            LGN(adder).patching_date=[char(batchopt.mouse{i})];
            LGN(adder).experimentator=experimentator;
            LGN(adder).cellname=n_str;
            LGN(adder).eye_inj_ord=batchopt.eye_inj_order{i};
            LGN(adder).slice_nr=batchopt.slice{i}(k);
            LGN(adder).brain_contra_ipsi=batchopt.injection_order{i}(k);
            LGN(adder).ocular_category=batchopt.category{i}(k);
            LGN(adder).MD=batchopt.MD{i}(k);
            LGN(adder).clear_sl_nr=batchopt.clear_sl_nr{i}(k);
            hem=[2 5 8 11 14 17 20 23 26 29];
            LGN(adder).hemisphere=batchopt.hemisphere{i}(hem(k):hem(k)+1);
            LGN(adder).photodiode_flag=batchopt.photodiode_flag{i}(k);
            
            %LGN(adder).imaged_slice_nr=
            %LGN(adder).hemisphere=
            %LGN(adder).eye_inj_ord=batchopt.eye_inj_ord{i}(k);
            %LGN(adder).eye_lid=
            
            %  dLGN_ephys.data{adder+1,2}=[char(batchopt.mouse{i}), fold_name];
            %                 dLGN_ephys.data{adder+1,3}=slice_nr;
            %                 dLGN_ephys.data{adder+1,4}=batchopt.injection_order{i}(k);
            %                 dLGN_ephys.data{adder+1,5}=batchopt.category{i}(k);
            %                 dLGN_ephys.data{adder+1,6}=batchopt.MD{i}(k);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ramp
            if analyze_ramp==1
                LGN(adder).step_red=red_ramp;
                LGN(adder).step_blue=blue_ramp;
                
                
                %                 % change these to using a stats test and recalc ODI and
                %                 % ratio (16 Juli 2019)!!!
                %                 LGN(adder).ODI_AMPA_step=ODI.ODI_AMPA;
                %                 LGN(adder).ODI_NMDA_step=ODI.ODI_NMDA;
                %                 LGN(adder).AMPA_NMDA_r_red=Ratio_AMPA_NMDA.R_r;
                %                 LGN(adder).AMPA_NMDA_r_blue=Ratio_AMPA_NMDA.R_b;
                
                % JB: recalculating ODI
                % perform nonparametric tests
                testAMPAblue = signrank(...
                    LGN(adder).step_blue.neg_mean2(2,end-5:end),...
                    LGN(adder).step_blue.neg_base_mean2(2,end-5:end))<0.05;
                
                testAMPAred = signrank(...
                    LGN(adder).step_red.neg_mean1(2,end-5:end),...
                    LGN(adder).step_red.neg_base_mean1  (2,end-5:end))<0.05;
                
                if ~isempty(LGN(adder).step_blue.pos_peak2) % if this is empty it is assumed +40 exp were not performed
                    testNMDAblue = signrank(...
                        LGN(adder).step_blue.pos_mean2(2,end-5:end),...
                        LGN(adder).step_blue.pos_base_mean2(2,end-5:end))<0.05;
                    
                    testNMDAred = signrank(...
                        LGN(adder).step_red.pos_mean1(2,end-5:end),...
                        LGN(adder).step_red.pos_base_mean1(2,end-5:end))<0.05;
                else
                    testNMDAblue = nan;
                    testNMDAred = nan;
                end
                
                LGN(adder).step_red.resp_tests.testAMPAblue = testAMPAblue;
                LGN(adder).step_red.resp_tests.testAMPAred = testAMPAred;
                LGN(adder).step_red.resp_tests.testNMDAblue = testNMDAblue;
                LGN(adder).step_red.resp_tests.testNMDAred = testNMDAred;
                
                %                 try;clf(20);end
                
                % AMPA ODI
                if ~testAMPAred && ~testAMPAblue % no red or blue
                    LGN(adder).ODI_AMPA_step=nan;
                    RespBA =nan;
                    RespRA = nan;
                    RespRA_peak = nan;
                    RespBA_peak = nan;
                else
                    RespB = abs(LGN(adder).step_blue.neg_mean2(2,:));
                    RespR = abs(LGN(adder).step_red.neg_mean1(2,:));
                    RespB_peak = abs(LGN(adder).step_blue.neg_peak2(2,:));
                    RespR_peak = abs(LGN(adder).step_red.neg_peak1(2,:));
                    BaseB = abs(LGN(adder).step_blue.neg_base_mean2(2,:));
                    BaseR = abs(LGN(adder).step_red.neg_base_mean1(2,:));
                    
                    
                    steps_used = find([RespR-BaseR] == max(RespR-BaseR));
                    RespRA = RespR(steps_used)-BaseR(steps_used);
                    RespBA = RespB(steps_used)-BaseB(steps_used); % baseline subtract (could also use the average across a few steps)
                    RespRA_peak = RespR_peak(steps_used)-BaseR(steps_used); % baseline subtracted peak of response (nessesary for some calculations)
                    RespBA_peak = RespB_peak(steps_used)-BaseB(steps_used);
                    
                    % if ether of the responses are negative
                    if RespRA < 0; RespRA = 0; end
                    if RespBA < 0; RespBA = 0; end
                    if RespRA_peak < 0; RespRA_peak = 0; end
                    if RespBA_peak < 0; RespBA_peak = 0; end
                    
                    %                     figure(20);
                    %                     subplot(2,1,1); plot(RespR-BaseR,'r'); hold on; plot(RespB-BaseB,'b');
                    %                     plot([steps_used steps_used], [RespB(steps_used)-BaseB(steps_used) RespR(steps_used)-BaseR(steps_used)],'k--')
                    
                    % ODI is based on response mean as this should be more stable
                    if testAMPAred && ~testAMPAblue % Contra
                        LGN(adder).ODI_AMPA_step=1;
                        RespBA = 0; % RespRA (stays the same)
                        RespBA_peak = 0; % RespRA_peak (stays the same)
                        
                    elseif  ~testAMPAred && testAMPAblue % no red
                        LGN(adder).ODI_AMPA_step=-1;
                        RespRA = 0;
                        RespRA_peak = 0;
                        steps_used = 6:11; % Taking the average here because there is no rational way of chooseing a single step. instead we take the average from the trials used for the responseiveness test testNMDAred
                        RespBA = nanmean(RespB(steps_used)-BaseB(steps_used));
                        RespBA_peak = nanmean(RespB_peak(steps_used)-BaseB(steps_used));
                        
                    elseif testAMPAred && testAMPAblue % Binocular
                        LGN(adder).ODI_AMPA_step = (nanmean(RespRA)-nanmean(RespBA))/(nanmean(RespRA)+nanmean(RespBA));
                        LGN(adder).ODI_AMPA_step_peak = (nanmean(RespRA_peak)-nanmean(RespBA_peak))/(nanmean(RespRA_peak)+nanmean(RespBA_peak));
                        
                        if LGN(adder).ODI_AMPA_step>1
                            LGN(adder).ODI_AMPA_step = 1;
                            LGN(adder).ODI_AMPA_step_peak = 1;
                            
                        elseif LGN(adder).ODI_AMPA_step<-1
                            LGN(adder).ODI_AMPA_step = -1;
                            LGN(adder).ODI_AMPA_step_peak = 1;
                            
                        end
                    end
                    
                    
                end
                LGN(adder).step_red.red_resp_AMPA = RespRA;
                LGN(adder).step_red.blue_resp_AMPA = RespBA;
                LGN(adder).step_red.red_resp_AMPA_peak = RespRA_peak;
                LGN(adder).step_red.blue_resp_AMPA_peak = RespBA_peak;
                LGN(adder).step_red.steps_use_AMPA = steps_used;
                
                % NMDA ODI
                if any(isnan([testNMDAblue, testNMDAred])) || all([~testNMDAred, ~testNMDAblue]) % no red or blue
                    LGN(adder).ODI_NMDA_step=nan;
                    RespBN = nan;
                    RespRN = nan;
                    RespBN_peak = nan;
                    RespRN_peak = nan;
                    steps_used = nan;
                    
                else
                    RespB = abs(LGN(adder).step_blue.pos_mean2(2,:));
                    RespR = abs(LGN(adder).step_red.pos_mean1(2,:));
                    RespB_peak = abs(LGN(adder).step_blue.pos_peak2(2,:));
                    RespR_peak = abs(LGN(adder).step_red.pos_peak1(2,:));
                    BaseB = abs(LGN(adder).step_blue.pos_base_mean2(2,:));
                    BaseR = abs(LGN(adder).step_red.pos_base_mean1(2,:));
                    
                    steps_used = find([RespR-BaseR] == max(RespR-BaseR));
                    RespRN = RespR(steps_used)-BaseR(steps_used);
                    RespBN = RespB(steps_used)-BaseB(steps_used); % baseline subtract
                    RespRN_peak = RespR_peak(steps_used)-BaseR(steps_used); % baseline subtracted peak of response (nessesary for some calculations)
                    RespBN_peak = RespB_peak(steps_used)-BaseB(steps_used);
                    
                    % if ehter of the responses a
                    if RespRN < 0; RespRN = 0; end
                    if RespBN < 0; RespBN = 0; end
                    if RespRN_peak < 0; RespRN_peak = 0; end
                    if RespBN_peak < 0; RespBN_peak = 0; end
                    
                    %                     subplot(2,1,2); plot(RespR-BaseR,'r'); hold on; plot(RespB-BaseB,'b');
                    %                     plot([steps_used steps_used], [RespB(steps_used)-BaseB(steps_used) RespR(steps_used)-BaseR(steps_used)],'k--')
                    
                    if testNMDAred && ~testNMDAblue % Contra
                        LGN(adder).ODI_NMDA_step=1;
                        RespBN = 0;% RespRN (stays the same)
                        RespBN_peak = 0;% RespRN_peak (stays the same)
                        
                    elseif ~testNMDAred && testNMDAblue % no red
                        LGN(adder).ODI_NMDA_step=-1;
                        RespRN = 0;
                        RespRN_peak = 0;
                        steps_used = 6:11; % Taking the average here because there is no rational way of chooseing a single step. instead we take the average from the trials used for the responseiveness test testNMDAred
                        RespBN = nanmean(RespB(steps_used)-BaseB(steps_used));
                        RespBN_peak = nanmean(RespB_peak(steps_used)-BaseB(steps_used));
                        
                    elseif testNMDAred && testNMDAblue % Binocular
                        LGN(adder).ODI_NMDA_step = (nanmean(RespRN)-nanmean(RespBN))/(nanmean(RespRN)+nanmean(RespBN));
                        LGN(adder).ODI_NMDA_step_peak = (nanmean(RespRN_peak)-nanmean(RespBN_peak))/(nanmean(RespRN_peak)+nanmean(RespBN_peak));
                        if LGN(adder).ODI_NMDA_step>1
                            LGN(adder).ODI_NMDA_step = 1;
                            LGN(adder).ODI_NMDA_step_peak  = 1;
                        elseif LGN(adder).ODI_NMDA_step<-1
                            LGN(adder).ODI_NMDA_step = -1;
                            LGN(adder).ODI_NMDA_step_peak  = -1;
                        end
                    end
                end
                LGN(adder).step_red.red_resp_NMDA = RespRN;
                LGN(adder).step_red.blue_resp_NMDA = RespBN;
                LGN(adder).step_red.red_resp_NMDA_peak = RespRN_peak;
                LGN(adder).step_red.blue_resp_NMDA_peak = RespBN_peak;
                LGN(adder).step_red.steps_use_NMDA = steps_used;
                
                % switch ODIs if contra~=red
                if ~LGN(adder).brain_contra_ipsi
                    LGN(adder).ODI_AMPA_step = LGN(adder).ODI_AMPA_step*-1;
                    LGN(adder).ODI_NMDA_step = LGN(adder).ODI_NMDA_step*-1;
                    LGN(adder).ODI_AMPA_step_peak = LGN(adder).ODI_AMPA_step_peak*-1;
                    LGN(adder).ODI_NMDA_step_peak = LGN(adder).ODI_NMDA_step_peak*-1;
                end
                
                % AMPA NMDA ratio
                if any(isnan([testAMPAblue, testAMPAred, testNMDAblue, testNMDAred])) || all([~testNMDAred, ~testNMDAblue])
                    LGN(adder).AMPA_NMDA_r_red=nan;
                    LGN(adder).AMPA_NMDA_r_blue=nan;
                else
                    LGN(adder).AMPA_NMDA_r_red=nanmean(RespRA_peak)/nanmean(RespRN_peak);
                    LGN(adder).AMPA_NMDA_r_blue=nanmean(RespBA_peak)/nanmean(RespBN_peak);
                end
                
            else
                LGN(adder).step_blue=[];
                LGN(adder).step_red=[];
            end
            
            
            
            
            
            
            %failures
            %             if analyze_mini==1
            %                 LGN(adder).ipsi_fail=Ipsi;
            %                 LGN(adder).contra_fail=Contra;
            %             else
            %                 LGN(adder).ipsi_fail=[];
            %                 LGN(adder).contra_fail=[];
            %             end
            if analyze_mini==1
                LGN(adder).red_failure_AMPA = red_failure_AMPA;
                LGN(adder).blue_failure_AMPA = blue_failure_AMPA;
                LGN(adder).red_failure_NMDA = red_failure_NMDA;
                LGN(adder).blue_failure_NMDA = blue_failure_NMDA;
            else
                LGN(adder).red_failure_AMPA = [];
                LGN(adder).blue_failure_AMPA = [];
                LGN(adder).red_failure_NMDA = [];
                LGN(adder).blue_failure_NMDA = [];
            end
            
            
            %Morphology
            if morpho==1
                LGN(adder).morphology=morphology;
                if ~isempty(morphology)==1
                    LGN(adder).DOi=morphology.DOi;
                    LGN(adder).DOi_min_pl=morphology.min_planes_nu;
                    LGN(adder).DOi_max_pl=morphology.max_planes_denom;
                else
                    LGN(adder).morphology=[];
                end
            end
            
            %scracm
            if scracm==1
                LGN(adder).scracm_blue=blue_scracm;
                LGN(adder).scracm_red=red_scracm;
            else
                LGN(adder).scracm_blue=[];
                LGN(adder).scracm_red=[];
            end
            %                 dLGN_ephys.data{adder+1,7}=blue_ramp;
            %                 dLGN_ephys.data{adder+1,8}=red_ramp;
            %                 dLGN_ephys.data{adder+1,14}=ODI.ODI_AMPA;
            %                 dLGN_ephys.data{adder+1,15}=ODI.ODI_NMDA;
            %                 dLGN_ephys.data{adder+1,16}=Ratio_AMPA_NMDA.R_r;
            %                 dLGN_ephys.data{adder+1,17}=Ratio_AMPA_NMDA.R_b;
            %                 else
            %                 dLGN_ephys.data{adder+1,7}=[];
            %                 dLGN_ephys.data{adder+1,8}=[];
            %                 end
            %                 %failures
            %                 if analyze_mini==1
            %                 dLGN_ephys.data{adder+1,9}=Ipsi;
            %                 dLGN_ephys.data{adder+1,10}=Contra;
            %                 else
            %                 dLGN_ephys.data{adder+1,9}=[];
            %                 dLGN_ephys.data{adder+1,10}=[];
            %                 end
            %                 %Morphology
            %                 if morpho==1
            %                 dLGN_ephys.data{adder+1,11}=morphology;
            %                 if ~isempty(morphology)==1;
            %                 dLGN_ephys.data{adder+1,18}=morphology.DOi;
            %                 else
            %                 dLGN_ephys.data{adder+1,18}=[];
            %                 end
            %                 end
            
            %                 %scracm
            %                 if scracm==1
            %                 dLGN_ephys.data{adder+1,12}=blue_scracm;
            %                 dLGN_ephys.data{adder+1,13}=red_scracm;
            %                 else
            %                 dLGN_ephys.data{adder+1,12}=[];
            %                 dLGN_ephys.data{adder+1,13}=[];
            %                 end
            
            
            adder=adder+1;
            
        end
        list=[];
        
    end
    LGN=LGN';
    % SAVE in analyzed directory
    if savefile==1
        cd(adata_dir);
        FileName=['Data_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName,'LGN','-v7.3');
        disp('FILE SAVED');
    else
        disp('FILE NOT SAVED');
    end
    
    %     if savefile==1
    %         cd(adata_dir);
    %         FileName=['Data_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
    %         save(FileName,'-struct','dLGN_ephys');
    %         disp('FILE SAVED');
    %     else
    %         disp('FILE NOT SAVED');
    %     end
    %% ALternative structure
    
    
    
    
    
end
%% FURTHER ANALYSIS/PLOTTING, Relation Morpho to other parameters

% if fanalysis==1
%     disp('dLGN data plotting of extracted parameters and calculation of ODI and AMPA/NMDA RATIOS');
%     adata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData\';%data directory of saved data
%
%     [ramps_peak ODI com]=dLGN_plot_analysis(adata_dir,0,1,0,0,0)
%
% end
%%

if fanalysis==1
    disp('dLGN data plotting of extracted Parameters and Morphology');
    adata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData\';%data directory of saved data
    
    LGN_fanalyis(adata_dir);
    
end

%%

if scracm_ana==1;
    disp('dLGN data plotting of extracted SCRACM parameters and calculation of ODI and AMPA/NMDA RATIOS');
    adata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData\';%data directory of saved data
    
    [ODI Morpho]=scracm_analysis(adata_dir,0);
    
end


%-------------------------------------------
