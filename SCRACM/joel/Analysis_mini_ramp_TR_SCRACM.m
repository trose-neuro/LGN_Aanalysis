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
analyze_mini=0;%flag if either mini only and/or ramp should be analyzed (1 or 0)
analyze_ramp=0;
analyze_spon=0;
morpho=0;%run morphology script
scracm=1;%run 2scracm script/analysis; REMEMBER TO SET RIGHT FLAGS
fanalysis=0;
factor=7;%std threshold factor
display=0;%flag to display plot (1 or 0)
ramp_rtrace=1;%save raw ephystraces or not (1 or 0)
savefile=1;%save file at the end or not
scracm_ana=0;%run 2scracm further analysis; ODI etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear LGN
if exist('user')
    user = 0;
end

if analyze_mini==1 || analyze_ramp==1 || morpho==1 || scracm==1
    disp('dLGN Analysis');
    
    sent = 'Which user data will be analyzed? type in 0 for SW or 1 for MF\n';%text appears in command window
    %      user = input(sent);%waiting for input which is either 0 or 1
    
    experimentator = 'SW';%default SW data
    if user==1
        experimentator = 'MF';%MF data
    end
    
    dLGN_ephys={};%empty structure for saving variables
    
    %     %%%%%%DIRECTORIES%%%%%%%
    rdata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\RawData';%data directory of raw data;change accordingly
    %     adata_dir         = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';%data directory of extracted date;change accordingly
    adata_dir         = 'C:\temp\LGN project';%data directory of extracted date;change accordingly
    %     ExpXls            = 'I:\Martin Fernholz\LGN _Project_Common\dLGN_ephys_analysis_excel_spread_sheet\Experiments_dLGN_2019_testingSW.xlsx';%directory where excel batch file is located;change accordingly
    ExpXls            = 'I:\Martin Fernholz\LGN _Project_Common\dLGN_ephys_analysis_excel_spread_sheet\Experiments_dLGN_2019_tr_scracm_testing.xlsx';%directory where excel batch file is located;change accordingly
    
    %     %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%DIRECTORIES - TR2019 %%%%%%%
    %     % TR2019: MAC: mount smb shares to /Volumes/first (easiest: do command+K in
    %     % finder)
    %     % 'smb://10G.ISI01.neuro.mpg.de/archive_bonhoeffer_group$/Simon Weiler/EXPLORER ONE'
    %     % 'smb://S15.neuro.mpg.de/R-bonhoe/Share/Simon/LGN_2019_SW_MF_JB_TR/dLGN_ephys_analysis_excel spread sheet/'
    %     rdata_dir         = '/Volumes/EXPLORER ONE/dLGN_rawDATA/';
    %     adata_dir         = '~/Analysis/dLGN_ephys_Analysis/';
    %     ExpXls            = '/Volumes/dLGN_ephys_analysis_excel spread sheet/Experiments_dLGN_SW.xlsx';
    %%%%%%%%%%%%%%%%%%%%%%
    
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
            %             if isempty(list)
            %                 try
            %                     dirinfo = dir(exp_folder{1});
            %                     dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
            %                     dirinfo(cellfun(@(x,c) str(),{dirinfo.name})
            %                     subdirinfo = cell(length(dirinfo));
            %                     for K = 1 : length(dirinfo)
            %                         thisdir = dirinfo(K).name;
            %                         subdirinfo{K} = dir(fullfile(thisdir, '*.mat'));
            %                     end
            %                 end
            %             end
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
            mapping=find(iterations==64 | iterations==49);
            
            
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
                morpho_folder=[char(exp_folder) '\' 'SingleTraces' char(batchopt.mouse{i}), fold_name];
                morpho_data=dir([char(exp_folder) '\' 'SingleTraces' char(batchopt.mouse{i}), fold_name '\*.swc']);
                if isempty(morpho_data)==0;
                    morphology=morpho_LGN(morpho_folder,morpho_data,batchopt.morphox{i}(k),batchopt.morphoy{i}(k),batchopt.morphoz{i}(k));
                else
                    disp('No morphology');
                    morphology=[];
                end
            end
            
            %% 2SCRACM
            if scracm==1
                %                 [blue_scracm, red_scracm]=LGN_scracm_new(list, mapping, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir,1);%use nested function rampanalysis
                [resp_maps, traces, illumination, resp_windows] ...
                    = LGN_scracm_new(list, mapping, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir,1); %use nested function rampanalysis
                %                 [blue_scracm, red_scracm]=LGN_scracm(list, mapping, clamp, exp_folder, factor,display,ramp_rtrace,user, filterephys,adata_dir,1);%use nested function rampanalysis
            end
            
            clamp=[];
            %%  prepare structure for all cells
            LGN(adder).animal_name=[char(batchopt.mouseID{i})];
            LGN(adder).patching_date=[char(batchopt.mouse{i})];
            LGN(adder).experimentator=experimentator;
            LGN(adder).cellname=n_str;
            LGN(adder).eye_inj_ord=batchopt.eye_inj_order{i};
            LGN(adder).slice_nr=batchopt.slice{i}(k);
            LGN(adder).brain_contra_ipsi=batchopt.injection_order{i}(k);
            try;LGN(adder).ocular_category=batchopt.category{i}(k);catch;LGN(adder).ocular_category=[];end
            try;LGN(adder).MD=batchopt.MD{i}(k);catch;LGN(adder).MD=[];end
            try;LGN(adder).clear_sl_nr=batchopt.clear_sl_nr{i}(k);catch;LGN(adder).clear_sl_nr=[];end
            hem=[2 5 8 11 14 17 20 23 26 29];
            try;LGN(adder).hemisphere=batchopt.hemisphere{i}(hem(k):hem(k)+1);catch;LGN(adder).hemisphere=[];end
            try;LGN(adder).photodiode_flag=batchopt.photodiode_flag{i}(k);catch;LGN(adder).photodiode_flag=[];end
            
            % plot scracm
            if scracm==1
                Map_trace_m = traces.Map_trace_m;
                figure('Name','30','units','normalized','outerposition',[0 0 1 1]); clf
                Rect    = [0.1, 0.1, 0.8, 0.8];
                AxisPos = myPlotPos(size(Map_trace_m,1), size(Map_trace_m,2), Rect);
                for j = 1:size(Map_trace_m,1)*size(Map_trace_m,1)
                    axes('Position', AxisPos(j, :))
                    row = floor((j-1)/size(Map_trace_m,1))+1;
                    col = rem((j-1),size(Map_trace_m,1))+1;
                    ylim([min(Map_trace_m(:)), max(Map_trace_m(:))]); hold on
                    xlim([0, size(Map_trace_m,3)])
                    xlimit = xlim; ylimit = ylim;
                    h = patch([0,xlimit(2)/2,xlimit(2)/2,0],[ylimit(1),ylimit(1),ylimit(2),ylimit(2)],[1,0,0],...
                        'FaceAlpha',0.2,'EdgeColor','none'); hold on
                    h = patch([xlimit(2)/2,xlimit(2),xlimit(2),xlimit(2)/2],[ylimit(1),ylimit(1),ylimit(2),ylimit(2)],[0,0,1],...
                        'FaceAlpha',0.2,'EdgeColor','none'); hold on
                    %         plot([1000 1000],ylim,'--k')
                    plot(squeeze(Map_trace_m(row,col,:)),'k'); axis off
                end
                suptitle([LGN(adder).patching_date(1:6) LGN(adder).experimentator num2str(n_str)])
                set(gcf,'color','w');
                
%                 figure;
%                 subplot(1,2,1); h = imagesc(resp_maps.red_scracm_AMPA_mean); title('red'); set(h,'AlphaData',~isnan(resp_maps.red_scracm_AMPA_mean));
%                 subplot(1,2,2); h = imagesc(resp_maps.blue_scracm_AMPA_mean); title('blue'); set(h,'AlphaData',~isnan(resp_maps.blue_scracm_AMPA_mean));
                fileloc = 'C:\temp\LGN project\scracm cell plots';
                saveas(gcf,[ fileloc '\' LGN(adder).patching_date(1:6) LGN(adder).experimentator num2str(n_str) '.png'])
                
            end
            
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
                % perform nonparametric tests % check if signrank is
                % correct or if ranksum should be used.
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
                    LGN(adder).ODI_AMPA_step_peak=nan;
                    RespBA =nan;
                    RespRA = nan;
                    RespRA_peak = nan;
                    RespBA_peak = nan;
                else
                    %                     RespB = abs(LGN(adder).step_blue.neg_mean2(2,:));
                    %                     RespR = abs(LGN(adder).step_red.neg_mean1(2,:));
                    %                     RespB_peak = abs(LGN(adder).step_blue.neg_peak2(2,:));
                    %                     RespR_peak = abs(LGN(adder).step_red.neg_peak1(2,:));
                    %                     BaseB = abs(LGN(adder).step_blue.neg_base_mean2(2,:));
                    %                     BaseR = abs(LGN(adder).step_red.neg_base_mean1(2,:));
                    
                    RespB = (LGN(adder).step_blue.neg_mean2(2,:))*-1;
                    RespR = (LGN(adder).step_red.neg_mean1(2,:))*-1;
                    RespB_peak = (LGN(adder).step_blue.neg_peak2(2,:))*-1;
                    RespR_peak = (LGN(adder).step_red.neg_peak1(2,:))*-1;
                    BaseB = (LGN(adder).step_blue.neg_base_mean2(2,:))*-1;
                    BaseR = (LGN(adder).step_red.neg_base_mean1(2,:))*-1;
                    
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
                        LGN(adder).ODI_AMPA_step_peak=1;
                        RespBA = 0; % RespRA (stays the same)
                        RespBA_peak = 0; % RespRA_peak (stays the same)
                        
                    elseif  ~testAMPAred && testAMPAblue % no red
                        LGN(adder).ODI_AMPA_step=-1;
                        LGN(adder).ODI_AMPA_step_peak=-1;
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
                    LGN(adder).ODI_NMDA_step_peak=nan;
                    RespBN = nan;
                    RespRN = nan;
                    RespBN_peak = nan;
                    RespRN_peak = nan;
                    steps_used = nan;
                    
                else
                    %                     RespB = abs(LGN(adder).step_blue.pos_mean2(2,:));
                    %                     RespR = abs(LGN(adder).step_red.pos_mean1(2,:));
                    %                     RespB_peak = abs(LGN(adder).step_blue.pos_peak2(2,:));
                    %                     RespR_peak = abs(LGN(adder).step_red.pos_peak1(2,:));
                    %                     BaseB = abs(LGN(adder).step_blue.pos_base_mean2(2,:));
                    %                     BaseR = abs(LGN(adder).step_red.pos_base_mean1(2,:));
                    
                    RespB = (LGN(adder).step_blue.pos_mean2(2,:));
                    RespR = (LGN(adder).step_red.pos_mean1(2,:));
                    RespB_peak = (LGN(adder).step_blue.pos_peak2(2,:));
                    RespR_peak = (LGN(adder).step_red.pos_peak1(2,:));
                    BaseB = (LGN(adder).step_blue.pos_base_mean2(2,:));
                    BaseR = (LGN(adder).step_red.pos_base_mean1(2,:));
                    
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
                        LGN(adder).ODI_NMDA_step_peak=1;
                        RespBN = 0;% RespRN (stays the same)
                        RespBN_peak = 0;% RespRN_peak (stays the same)
                        
                    elseif ~testNMDAred && testNMDAblue % no red
                        LGN(adder).ODI_NMDA_step=-1;
                        LGN(adder).ODI_NMDA_step_peak=-1;
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
                %                 LGN(adder).scracm_blue=blue_scracm;
                %                 LGN(adder).scracm_red=red_scracm;
                
                LGN(adder).scracm_red.resp_maps.ampa=resp_maps.red_scracm_AMPA_mean;
                LGN(adder).scracm_red.resp_maps.nmda=resp_maps.red_scracm_NMDA_mean;
                LGN(adder).scracm_blue.resp_maps.ampa=resp_maps.blue_scracm_AMPA_mean;
                LGN(adder).scracm_blue.resp_maps.nmda=resp_maps.blue_scracm_NMDA_mean;
                LGN(adder).scracm_red.traces = traces;
                LGN(adder).scracm_red.illumination = illumination;
                LGN(adder).scracm_red.resp_windows =  resp_windows;
                
                figure;
                subplot(2,2,1); imagesc(resp_maps.red_scracm_AMPA_mean);
                title('red'); colorbar; ylabel('AMPA'); pbaspect([1,1,1]);
                subplot(2,2,3); imagesc(resp_maps.red_scracm_NMDA_mean);
                colorbar; ylabel('NMDA'); pbaspect([1,1,1]);
                subplot(2,2,2); imagesc(resp_maps.blue_scracm_AMPA_mean);
                title('blue'); colorbar; pbaspect([1,1,1]);
                subplot(2,2,4); imagesc(resp_maps.blue_scracm_NMDA_mean);
                colorbar; pbaspect([1,1,1]);
                suptitle('Average response maps, across reps')
                
                figure;
                subplot(2,2,1); imagesc(mean(illumination.PD1_r_AMPA,3)./max(max(mean(illumination.PD1_r_AMPA,3)))*100);
                title('red'); colorbar; ylabel('AMPA'); pbaspect([1,1,1]);
                subplot(2,2,3); imagesc(mean(illumination.PD1_r_NMDA,3)./max(max(mean(illumination.PD1_r_NMDA,3)))*100);
                colorbar; ylabel('NMDA'); pbaspect([1,1,1]);
                subplot(2,2,2); imagesc(mean(illumination.PD2_b_AMPA,3)./max(max(mean(illumination.PD2_b_AMPA,3)))*100);
                title('blue'); colorbar; pbaspect([1,1,1]);
                subplot(2,2,4); imagesc(mean(illumination.PD2_b_NMDA,3)./max(max(mean(illumination.PD2_b_NMDA,3)))*100);
                colorbar; pbaspect([1,1,1]);
                suptitle('Average photodiode signal, across reps')
                
%                 autoArrangeFigures; pause(5)
                close all
            else
                %                 LGN(adder).scracm_blue=[];
                %                 LGN(adder).scracm_red=[];
                
                LGN(adder).scracm_red.resp_maps.ampa=[];
                LGN(adder).scracm_red.resp_maps.nmda=[];
                LGN(adder).scracm_blue.resp_maps.ampa=[];
                LGN(adder).scracm_blue.resp_maps.nmda=[];
                LGN(adder).scracm_red.traces=[];
                LGN(adder).scracm_red.illumination=[];
                LGN(adder).scracm_red.map_stim_orders=[];
                
            end

            
            adder=adder+1;
            
        end
        list=[];
        
    end
    LGN=LGN';
    % SAVE in analyzed directory
    if savefile==1
        cd(adata_dir);
        if scracm & ~analyze_mini & ~analyze_ramp
            FileName=['Data_S2CRACM_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
        else
            FileName=['Data_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
        end
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
