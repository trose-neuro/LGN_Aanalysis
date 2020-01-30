function [ana_data roi_data] = load_tseries_batch_ax(yearmonthday, mouse,  exp, reanalyze, datapath, varargin)

% AAHHH DAMN. I dont have time to improve the code. So I'm just doing it again!

% control + '=' to fold all
% control + shift + '.' to unfold level
%tesst
%% switchboard



close all;

loadflag = '_TR_';
adder ='_TR_';


delmishaps = 1;
delcrazies = 1;
cra = 3000;
darkframe_subtract_2nd = 1; %also includes bleedthrough correction;
green2red_bleed = 0.05;
neuropilfct = 0.7;

savedata = 0;
if ~savedata
    disp('DATA WILL NOT BE SAVED!!!!')
end

refit = 0;
fixfit = 0;
oopsi = 0 ;
ratio = 0;

% [xls_num,xls_txt]=xlsread('R:\Share\Juliane\Experiments_TR.xlsx');
%         [xls_num,xls_txt]=xlsread('R:\Share\trose\PROJECTS\Experiments_awake.xlsx');
% ExpXls            = 'R:\Share\Juliane\Experiments_BaselineMuscimol.xlsx';
ExpXls            = 'R:\Share\Juliane\Experiments_MDMuscimol.xlsx';

% reanalyze = 1;
if ratio ==  0
    disp('RATIO SWITCHED OFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('RATIO SWITCHED OFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('RATIO SWITCHED OFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end


overwrite =1;

if isempty(varargin)
    firstana = 0; %!!!ONLY RUN FOR VERY FIRST ANALYSIS OF ROIS!
    batchcall = 0;
    fitcaran = 0;
    PSTHs = 0;
    runexclude = 0;
else
    firstana = varargin{1};
    batchcall = varargin{2};
    fitcaran = varargin{3};
    PSTHs =  varargin{4};
    runexclude = varargin{5};
    %     plotdata = 1; % plot single timepoint PSTH during reanalysis (in make_ori_fct_batch call)
end

if exp(1) > 60000
    bscope2 = 1;
    runextract =1;
else
    bscope2 = 0;
    runextract = 0;
    runexclude = 0;
end


plotdata = 1; % plot single timepoint PSTH during reanalysis (in make_ori_fct_batch call)

loadaux = 1;

if firstana
    loadaux = 1;
end
if PSTHs
    loadaux = 1;
end
if refit
    loadaux = 1;
end
% reanalyze =1;
explorez =0;
% loadaux = 0;
% firstana =1;

%  fitcaran = 0;
pixelmaps = 1; %make cell-sorted maps
sortmap = 1;
responly = 1; % if 1, include only cells responsive throughout for pixelmap!
baselinerois = 0;% if 1, include all cells responsive during baseline for pixelmap!
sigrespmap = 0;

roimaps = 0;
plotarray_OD_map = 1;
baselinerois_map = 0;
baselinerois_all_map = 0;
baselinerois_md_all_map = 0;
collapse = 1; %collapse volume stack to one level for ROImaps

% PSTHs = 0; %show and save longitudinal PSTHs

recalc_circ  = 1; %Just % in hack to recalc the circul'ar variance
bootstrap_var = 0;

fitamp  = 0;

noplot = 0;
justload = 0;

try
    closefigs;
end

%% prepare batch-load variables
if isempty(yearmonthday)
    %date filter not fully implemented yet
    yearmonthday = {''};
end

if isempty(exp)
    startexp =1;
    endexp = 99999;
    selected = 0;
elseif length(exp) > 1;
    selected = 1;
else
    startexp = exp-1;
    endexp = exp+1;
    selected = 0;
end

PC = getenv('computername');

adata_dir       = 'I:\Juliane Jaepel-Schael\AnalyzedData\';

if isempty(mouse)
    %multiple mice filter not fully implemented yet
    mouse = 'JJ_'
end


try
    cd(adata_dir)
catch
    pause(1)
    cd(adata_dir)
end

str = mouse;

%% load excel sheet with experiment parameters  (baselinse sessions,
%number of MDs etc.)


[xls_num,xls_txt]=xlsread(ExpXls);

mouse_row = find(strcmp(xls_txt(:,1), mouse)==1);
exp_row = find(xls_num(:,3)==exp(1)); %assuming that the first experiment in the list is the samesite_id
baseline = xls_num(exp_row,1);
recovery1 = xls_num(exp_row,2);
baseline2 = 0;
if isnan(baseline2)
    baseline2 = 0;
end

if ~isempty(exp_row)
    
    excludeids= xls_txt(exp_row+1,6);
    
    if ~isempty(excludeids{1})
        exclude_sessions = eval((excludeids{1}));
    else
        exclude_sessions =[];
    end
    
    baseids= xls_txt(exp_row+1,12);
    baseline_pair = eval((baseids{1}));
    
    baseids14 = xls_txt(exp_row+1,13);
    baseline_pair14 = eval((baseids14{1}));
    
    loaddrive = xls_txt(exp_row+1,15);
    
    recovery_incl= xls_num(exp_row,4);
    md_incl = xls_num(exp_row,7);
    % baseids= xls_txt(exp_row+1,12);
    
    
    l = xls_num(exp_row,15);
    
    if ~isempty(exclude_sessions)
        disp(['EXCLUDING SESSION ' ])
        disp( exp(exclude_sessions))
        orig_exp = exp;
        exp(exclude_sessions) =[];
    else orig_exp = exp;
    end
    
else
    
    baseids= 1;
    baseline_pair = 1;
    
    baseids14 = [];
    baseline_pair14 = [];
    
    loaddrive = 'I';
    
    recovery_incl= 0;
end
% baseids= xls_txt(exp_row+1,12);


%% batch load data
h = waitbar(0,'Loading timeseries data, auxdata and analysis...');
for reprun = 1:length(str)
    
    try
        a =dir (['*' mouse{1} ]);
    catch
        a =dir (['*' mouse ]);
    end
    s =[];
    
    for il = 1:1:length(a);
        s{il} =rdir([cd '\' a(il).name '\' '*\*.mat']);
    end
    
    for kl = 1:1:length(s);% length(s):-1:1
        
        if ~isempty(s{kl})
            expnumbs = regexp({s{kl}.name},'(?<=Adata-)[0-9]*', 'match', 'once')'; %note the 'once option to get rid of cell arrays of cells and make cellfun work
            [idxs b] = unique(expnumbs);
            
            b = b(ismember(idxs,num2str(exp))); %limit exp indices to selected experiments...
            idxs = idxs(ismember(idxs,num2str(exp))); %limit exp indices to selected experiments...
            
            if isempty(idxs{1});
                idxs(1) = [];
            end
            
            ct = 1;
            for il =1:length(idxs)
                waitbar(il / length(idxs));
                
                if ~selected
                    if str2num(idxs{il}) <= startexp || str2num(idxs{il}) >= endexp;
                        continue
                    end
                elseif sum(str2num(idxs{il}) == exp) == 0;
                    continue
                    %                 elseif
                end
                
                try
                    cd(fileparts(s{kl}(b(il)).name));
                catch
                    pause(1)
                    cd(fileparts(s{kl}(b(il)).name));
                end
                roi_load = dir(['*' loadflag '*' idxs{il} '.mat']);
                roi_load = roi_load.name;
                disp(' - - - loading roi_data')
                roi_data{il} = load(roi_load, '-mat');
                
                %check for duplicate droups (non zero)
                check_duplicate_group(roi_data{il}.ROIs);
                
                % Spine data?
                spinect = sum(strcmp({roi_data{1}.ROIs(:).typename}', 'Spine'));
                if spinect > .25 * length(roi_data{il}.ROIs);
                    spine_ana = 1;
                    disp('More than 1/4th of the ROIs are Spine ROIS. Switching to spine extraction!')
                else
                    %                     disp('Less than 1/4th of the ROIs are Spine ROIS. Switching to cell extraction!')
                    spine_ana = 0;
                end
                disp(' ')
                disp(' - - - - - - - - - - - ');
                disp(['EXP ' idxs{il}]);
                disp(' - - - - - - - - - - - ');
                
                
                % load further ROI info if not in structure
                dirparts{il} = regexp(s{kl}(b(il)).name,['[^\\]*'], 'match');
                dirparts2{il} = regexp(s{kl}(b(il)).name,['[^\\]*'], 'match');
                
%                 if ~strcmp(PC, 'P1-364')
%                     dirparts{il}(5) = dirparts{il}(6) ;
%                 end
                
                if loadaux && firstana
                    if ~overwrite
                        
                        
                        %cheap trick to start off with the  primary
                        %ana_file generation where we stopped the last time
                        %ORIanalyis_z8_oopsi_full
                        %ORIanalyis_z8_meds_trace
                        %ORIanalyis_z8_meds
                        %ORIanalyis_z8_meds_trace_fit_mean
                        %                         ORIanalyis_z8_meds_trace_meanrun_fitfix
                        %                         ORIanalyis_z8_meds_trace_fit_mean.mat
                        ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit.mat']);
                        
                        %                         if oopsi
                        %                             ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_oopsi_full']);
                        %                         end
                        %                         ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z12_full.mat']);
                        ana_load = [cd '\Analysis\exp' idxs{il} '\' ana_load.name];
                        ana_data{il} = load(ana_load, '-mat');
                        disp(['Loading: ' ana_load '.mat']);
                        
                    else
                        
                        try
                            auxstartdir = [datapath  a(1).name '\ImagingData\' dirparts{il}{5}];
                            cd(auxstartdir);
                        catch
                            auxstartdir = [datapath  a(1).name '\ImagingData\' dirparts{il}{6}];
                            cd(auxstartdir);
                        end
                        %get eye number
                        try
                            diri = cd;
                            dat = diri(end-9:end);
                            stimdir=['..\..\_stim\' dat '\'];
                            stimfile = dir([stimdir '\*' idxs{il} '*.mat*']);
                            stimarrayT = load([stimdir stimfile(1).name], '-mat');
                            eyenum = sum([stimarrayT.calc_eyes stimarrayT.binocular])+1;
                            
                            if stimarrayT.awake_eyes
                                eyenum = 2;
                            end
                        catch
                            eyenum = 2;
                        end
                        [auxdata{il}, ids{il}, ids_new{il}, frame_times{il}, stimarray{il}] = getauxstim(str2num(idxs{il}),cd, eyenum, roi_data{1}.info.level, bscope2); %disp('assuming 4 levels');
                        
                        if ~isfield(roi_data{il}, 'info')
                            files_t = dir([diri '\*' idxs{il} '*.tif*']);
                            tio = Tiff([diri '\' files_t(1).name], 'r');
                            roi_data{il}.info.ImageDescription = tio.getTag('ImageDescription');
                            clear tio
                        end
                        SamplingFreq2 = regexp(roi_data{il}.info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
                        roi_data{il}.info.SamplingFreq1 = str2num(SamplingFreq2{1}) / roi_data{1}.info.level;
                        
                        
                        if runextract
                            [Position, Speed, running_dist, running_maxspeed, runbinary_dist, runbinary_maxspeed ] = getrunningtrials(auxdata{il}, ids_new{il}, frame_times{il}, roi_data{1}.info.level,  roi_data{il}.info.SamplingFreq1);
                            ids_new{il}(1).running_dist = running_dist(:,:,1);
                            ids_new{il}(1).running_maxspeed = running_maxspeed(:,:,1);
                            ids_new{il}(2).running_dist = running_dist(:,:,2);
                            ids_new{il}(2).running_maxspeed = running_maxspeed(:,:,2);
                            ids_new{il}(1).Position = Position;
                            ids_new{il}(1).Speed = Speed;
                            ids_new{il}(2).Position = Position;
                            ids_new{il}(2).Speed = Speed;
                            
                            
                        end
                        
                        % run initial analysis get auxdata;
                        disp(' - - running initial analysis');
                        
                        if delmishaps
                            roi_data{il}.nan_idxs = {roi_data{il}.nan_idxs}
                            tp = correctroidata(roi_data(il), roi_data{il}.nan_idxs,roi_data{il}.info.level);
                            roi_data{il} = tp{1};
                        end
                        
                        make_ori_figures_fct(adata_dir, mouse, dat, exp(il), roi_data{il}.ROIs, roi_data{il}.np, neuropilfct, ids_new{il} , 'none', roi_data{il}.info, plotdata,0, ratio, runexclude);
                                                 disp(' - - returning to master batch');
                        continue
                    end
                end
                
                cd(fileparts(s{kl}(b(il)).name));
                
                %CHANGE FOR DIFFERENT PRESELECTION CRITERION
                
                if ~runexclude
                    ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit.mat']);
                    %                                         ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg.mat']);
                else
                    ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_runexclude.mat']);
                end
                if spine_ana
                    ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_runexclude_SPINE.mat']);
                end
                if oopsi
                    ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_oopsi_full.mat']);
                end
                %                 ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z12_full.mat']);
                ana_load = [cd '\Analysis\exp' idxs{il} '\' ana_load.name];
                ana_data{il} = load(ana_load, '-mat');
                disp(['Loading: ' ana_load '.mat']);
                
                
                if fixfit
                    
                    ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit.mat']);
                    
                    %                 ana_load = dir([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z12_full.mat']);
                    ana_load = [cd '\Analysis\exp' idxs{il} '\' ana_load.name];
                    ana_data2{il} = load(ana_load, '-mat');
                    ana_data2{il}.Fit = ana_data{il}.Fit
                    peaks = ana_data2{il}.peaks;
                    Fit = ana_data2{il}.Fit;
                    info = ana_data2{il}.info;
                    disp(['FIXING FIT OF: ' ana_load '.mat']);
                    save([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_runexclude_fitfix.mat'], 'peaks',  'Fit', 'info')
                    %                     return
                end
                %                 Fit= ana_data{il}.Fit;
                %                                         save([cd '\Analysis\exp' idxs{il} '\ORIanalyis_z8_meds_trace.mat' ], 'Fit', '-append')
                %                         disp(['Appended Fit: ' ana_load '.mat']);
                % %                         return
                %
                
                if ~isfield(roi_data{il}.ROIs, 'indices')
                    %extract date from adata path
                    
                    ROIname = [datapath  a(1).name '\ImagingData\' dirparts{il}{5} '\Pixelmaps\exp' idxs{il} '\templates\ROIs\ROI.mat'];
                    load(ROIname, '-mat')
                    
                    disp(['Loading :' ROIname ]);
                    
                    for pp = 1:length(ROI)
                        roi_data{il}.ROIs(pp).indices = sub2ind([1024 1024],ROI(pp).body(:,2),ROI(pp).body(:,1));
                        roi_data{il}.ROIs(pp).perimeter = ROI(pp).perimeter;
                    end
                end
                
                % get imageinfo if necessary
                if ~isfield(ana_data{il}, 'info')
                    disp(['had to grab image info for recording ' num2str(il)])
                    ainfof = dir([datapath  a(1).name '\ImagingData\' dirparts{il}{5} '\exp' num2str(idxs{il}) '*001.tif']);
                    ainfo = imfinfo([datapath  a(1).name '\ImagingData\' dirparts{il}{5} '\' ainfof.name]);
                    ana_data{il}.info = ainfo(1);
                end
                
                %extract triggertime and level info
                trigtime = regexp(ana_data{il}(1).info.ImageDescription,'(?<=triggerClockTimeFirst = )\S+[\w]+[ ]+[\w]+\S+[\n]','match');
                level = str2num(cell2mat(regexp(ana_data{il}.info.ImageDescription,'(?<=stackNumSlices = )[0-9]*', 'match')));
                ana_data{il}.triggertime = datenum(trigtime{1}(2:end), 'dd-mm-yyyy HH:MM:SS.FFF');
                ana_data{il}.level = level;
                
                % get auxdata;
                if loadaux
                    
                    auxstartdir = [datapath  a(1).name '\ImagingData\' dirparts{il}{5}];
                    try
                        pause(1);
                        cd(auxstartdir);
                    catch
                        pause(1)
                        auxstartdir = [datapath  a(1).name '\ImagingData\' dirparts{il}{6}];
                        cd(auxstartdir)
                    end
                    try
                        diri = cd;
                        dat = diri(end-9:end);
                        stimdir=['..\..\_stim\' dat '\'];
                        try
                            stimfile = dir([stimdir '\*' idxs{il} '*.mat*']);
                            stimarrayT = load([stimdir stimfile(1).name], '-mat');
                        catch
                            stimfile = dir([stimdir '\*' idxs{il} '*.']);
                            stimarrayT = load([stimdir stimfile(1).name], '-mat');
                        end
                        eyenum = sum([stimarrayT.calc_eyes stimarrayT.binocular])+1;
                        try
                            if stimarrayT.awake_eyes
                                eyenum = 2;
                            end
                        end
                    catch
                        eyenum = 2;
                    end
                    disp(' - - Loading AUX Data')
                    try
                        [auxdata{il}, ids{il}, ids_new{il},  frame_times{il}, stimarray{il}] = getauxstim(str2num(idxs{il}),cd, eyenum, level, bscope2);
                    catch
                        [auxdata{il}, ids{il}, ids_new{il},  frame_times{il}, stimarray{il}] = getauxstim(str2num(idxs{il}),cd, eyenum, level, bscope2);
                    end
                    SamplingFreq2 = regexp(roi_data{il}.info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
                    try
                        roi_data{il}.info.SamplingFreq1 = str2num(SamplingFreq2{1}) / roi_data{1}.info.level;
                    catch
                        roi_data{1}.info.level = level
                        roi_data{il}.info.SamplingFreq1 = str2num(SamplingFreq2{1}) / roi_data{1}.info.level;
                    end
                    
                    if runextract
                        disp(' - - Extracting Running info')
                        [Position, Speed, running_dist, running_maxspeed, runbinary_dist, runbinary_maxspeed ] = getrunningtrials(auxdata{il}, ids_new{il}, frame_times{il}, roi_data{1}.info.level,  roi_data{il}.info.SamplingFreq1);
                        ids_new{il}(1).running_dist = running_dist(:,:,1);
                        ids_new{il}(1).running_maxspeed = running_maxspeed(:,:,1);
                        ids_new{il}(2).running_dist = running_dist(:,:,2);
                        ids_new{il}(2).running_maxspeed = running_maxspeed(:,:,2);
                        ids_new{il}(1).Position = Position;
                        ids_new{il}(1).Speed = Speed;
                        ids_new{il}(2).Position = Position;
                        ids_new{il}(2).Speed = Speed;
                        
                        try
                            close(854853)
                        catch
                            figure(854853);
                        end
                        %                             plot(auxdata{il}(15,frame_times{il}(1:4:end)))
                        hold all
                        plot(auxdata{il}(8,frame_times{il}(1:4:end)))
                        plot(Position/30)
                        plot(Speed*2)
                        plot(runbinary_maxspeed,'-k')
                        plot(runbinary_dist/2,':k')
                        xlabel('sample #');ylabel('arbitrary units');
                        dist_crit = 1;  %distance covered (a.u.) to count as running trial
                        speed_crit = 0.4; %corresponds to 3.5 cm/s
                        mintrials = 4;
                        rejct = sum(running_maxspeed>speed_crit,2);
                        overmin = sum(rejct(:)>mintrials);
                        rejct = sum(rejct(:));
                        prct = rejct * 100 / length(running_maxspeed(:));
                        
                        title([mouse{1} ' - ' num2str(prct) '% rejected - ' num2str(overmin) ' oris > 4 trials rejected'], 'Interpreter', 'none');
                        disp(['Excluded ' num2str(prct) '% (' num2str(rejct) ' trials out of ' num2str(length(running_maxspeed(:))) ') because of animal running' ]);
                        disp(['Of these ' num2str(overmin) ' orientations had <= 4 trials'])
                        saveas(gcf, [ fileparts(s{kl}(b(il)).name) '\exp' num2str(idxs{il}) 'run_extraction.fig']);
                        
                        ids_new{il}(1).exclprc = prct;
                        ids_new{il}(2).exclprc = prct;
                    end
                end
                
                roi_data{il}.mishap_idxs = find_mishap(roi_data{il}.nan_idxs, 150);
                
                if ~isempty(roi_data{il}.mishap_idxs)
                    set(gcf, 'Name', num2str(idxs{il}), 'NumberTitle', 'off');
                    try
                        disp(['Shutter mishap at exp ' num2str(idxs{il}) ' idx: ' num2str(roi_data{il}.mishap_idxs(:,1)) ' to ' num2str(roi_data{il}.mishap_idxs(:,2))])
                    catch
                        disp(['MULTIPLE MISHAPS! CHECK EXP: ' num2str(idxs{il})])
                    end
                end
                
                if ~isempty(roi_data{il}.mishap_idxs) && delmishaps
                    kill =  round(roi_data{il}.mishap_idxs./level);
                    for loop = 1:length(roi_data{il}.ROIs)
                        roi_data{il}.ROIs(loop).activity(kill(1):kill(2)) = NaN;
                        roi_data{il}.np(loop).activity(kill(1):kill(2)) = NaN;
                        roi_data{il}.ROIs(loop).activity_r(kill(1):kill(2)) = NaN;
                        roi_data{il}.np(loop).activity_r(kill(1):kill(2)) = NaN;
                    end
                end
                
                if darkframe_subtract_2nd
                    disp(' - - Secondary Darkframe extraction and bleedthrough correction')
                    medgreen = median(reshape([roi_data{il}.ROIs(:).activity],size(roi_data{il}.ROIs(1).activity,2),size(roi_data{il}.ROIs(:),1)),2);
                    if ratio
                        medred = median(reshape([roi_data{il}.ROIs(:).activity_r],size(roi_data{il}.ROIs(1).activity_r,2),size(roi_data{il}.ROIs(:),1)),2);
                    end
                    
                    for loop = 1:length(roi_data{il}.ROIs)
                        if isfield(roi_data{il}.info, 'darkframes')
                            Dframes = [ceil(roi_data{il}.info.darkframes(1)/4)+3: floor(roi_data{il}.info.darkframes(end)/4)-3];
                            
                            %take just the last 5 frames
                            Dframes = Dframes(end-4:end);
                            
                            DCg = mean(medgreen(Dframes));
                            %                             NPg = mean(roi_data{il}.np(loop).activity(Dframes));
                            if ratio
                                DCr = mean(medred(Dframes));
                                %                                 NPr = mean(roi_data{il}.np(loop).activity_r(Dframes));
                            end
                            roi_data{il}.ROIs(loop).activity = roi_data{il}.ROIs(loop).activity - DCg;
                            roi_data{il}.np(loop).activity = roi_data{il}.np(loop).activity - DCg;
                            
                            if ratio
                                roi_data{il}.ROIs(loop).activity_r = roi_data{il}.ROIs(loop).activity_r - DCr;
                                roi_data{il}.np(loop).activity_r = roi_data{il}.np(loop).activity_r - DCr;
                            end
                            
                        else %disp('Detecting darkframe in mean fluorescence trace');
                            try
                                dfenabled = str2num(cell2mat(regexp(roi_data{il}.info.ImageDescription,'(?<=scanimage.SI4.userFunctionsCfg__3.Enable = )[a-z]*', 'match')));
                                if dfenabled
                                    dframes = str2num(cell2mat(regexp(roi_data{il}.info.ImageDescription,'(?<=scanimage.SI4.userFunctionsCfg__3.Arguments = {)[0-9]*', 'match')));
                                    dframes = round(dframes/level)+1;
                                    
                                    dfstart = find(abs(diff(medgreen)) > 4* nanstd(abs(diff(medgreen))));
                                    dfstart = dfstart(dfstart>1);
                                    Dframes = dfstart(1)+1:dfstart(1)+1+dframes*.95;
                                    %take just the last 5 frames
                                    try
                                        Dframes = Dframes(end-4:end);
                                    end
                                    
                                    DCg = mean(medgreen(Dframes));
                                    %                             NPg = mean(roi_data{il}.np(loop).activity(Dframes));
                                    if ratio
                                        DCr = mean(medred(Dframes));
                                        %                                 NPr = mean(roi_data{il}.np(loop).activity_r(Dframes));
                                    end
                                    roi_data{il}.ROIs(loop).activity = roi_data{il}.ROIs(loop).activity - DCg;
                                    roi_data{il}.np(loop).activity = roi_data{il}.np(loop).activity - DCg;
                                    if ratio
                                        roi_data{il}.ROIs(loop).activity_r = roi_data{il}.ROIs(loop).activity_r - DCr;
                                        roi_data{il}.np(loop).activity_r = roi_data{il}.np(loop).activity_r - DCr;
                                    end
                                else disp('NO DARKFRAME CORRECTION!');
                                end
                            catch
                                disp('NO DARKFRAME CORRECTION!');
                            end
                            
                            
                        end
                        if loop == 1;
                            llhd = figure(246246264);
                            plot(roi_data{il}.ROIs(loop).activity, 'g'); hold all;
                            plot(roi_data{il}.ROIs(loop).activity_r, 'r'); hold all;
                        end
                        %                         hline(0)
                        %green to red bleedthough correction
                        try
                            roi_data{il}.ROIs(loop).activity_r = roi_data{il}.ROIs(loop).activity_r - green2red_bleed * roi_data{il}.ROIs(loop).activity;
                            roi_data{il}.np(loop).activity_r = roi_data{il}.np(loop).activity_r - green2red_bleed * roi_data{il}.np(loop).activity;
                        end
                        
                    end
                end
                ct = ct+1;
            end
        end
    end
end
close(h);


if firstana
    return
end


try
    try
        
        saveas(llhd, [adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) '_df_traces.png'], 'png');
    catch
        try
            saveas(llhd, [adata_dir '\Summaries\' mouse{1} '_site' num2str(exp(1)) '_df_traces.png'], 'png');
        catch
            saveas(llhd, [adata_dir '\Summaries\' mouse '_site' num2str(exp(1)) '_df_traces.png'], 'png');
            %         mouse{1} = mouse;
        end
    end
    % hl
    
    try
        close(llhd);
    end
end

% remove empty cells
ana_data = ana_data(~cellfun('isempty',ana_data));
roi_data = roi_data(~cellfun('isempty',roi_data));
dirparts = dirparts(~cellfun('isempty',dirparts));

if loadaux
    ids = ids(~cellfun('isempty',ids));
    ids_new = ids_new(~cellfun('isempty',ids_new));
    
    stimarray = stimarray(~cellfun('isempty',stimarray));
    assignin('base', 'ids',ids);
    assignin('base', 'ids_new',ids_new);
    assignin('base', 'stimarray',stimarray);
    
end

% assign to workspace
assignin('base', 'roi_data',roi_data);
assignin('base', 'ana_data',ana_data);
assignin('base', 'baseline',baseline);
assignin('base', 'recovery1',recovery1);


if justload || firstana
    return
end
if fixfit
    return
end

if reanalyze
    loadaux = 1;
    %     h = waitbar(0,'reanalyzing data...');
end

%% [+] - - - - ANALYSIS MASTERLOOP through sessions - - - - -
cl  = 1;
colors =([0 0 1;1 0 0; 0 1 0]);

if ~noplot
    figure(352253)
end

% THE MASTERLOOP STARTS HERE!
for i = 1:length(roi_data) %!!!!!!!
    try
        yearmonthday = [];
        for pop = 1:length(dirparts);
            yearmonthday{pop} = strcat( [dirparts{pop}{5}]);
        end
        
    end
    SamplingFreqt = regexp(ana_data{i}.info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
    SamplingFreq1(i) = str2num(SamplingFreqt{1}) / level;
    [~,~,~,~,px] = getzoom([],16,ana_data{i}.info);
    
    %% reanalysis?
    if reanalyze && fitcaran
        disp(' - - Fitting Orientation Data')
        if ~refit
            try
                FitIpsi = [ana_data{i}.Fit.ipsi];
                ana_data{i}.peaks(1).Tune_Anova_maxAmpDelAve_contra
                disp(['skipping reanalysis of exp' num2str(exp(i)) ' -> fit and tune already done'])
            end
        else
            %             waitbar(i / length(roi_data));
            disp(['performing reanalysis of exp' num2str(exp(i))])
            try
                make_ori_figures_fct(adata_dir, mouse, yearmonthday{i}, exp(i), roi_data{i}.ROIs, roi_data{i}.np, neuropilfct, ids_new{i} , 'none', ana_data{i}.info, plotdata, fitcaran, ratio, runexclude);
            catch
                make_ori_figures_fct(adata_dir, mouse, yearmonthday(i,:), exp(i), roi_data{i}.ROIs, roi_data{i}.np, neuropilfct, ids_new{i} , 'none', ana_data{i}.info, plotdata, fitcaran, ratio, runexclude);
            end
        end
    elseif reanalyze
        %         try
        make_ori_figures_fct(adata_dir, mouse, yearmonthday{i}, exp(i), roi_data{i}.ROIs, roi_data{i}.np, neuropilfct, ids_new{i} , 'none', ana_data{i}.info, plotdata, fitcaran, ratio, runexclude);
        
        %         catchr
        %             make_ori_figures_fct(adata_dir, mouse, yearmonthday(i,:), exp(i), roi_data{i}.ROIs, roi_data{i}.np, neuropilfct, ids_new{i} , 'none', ana_data{i}.info, plotdata, fitcaran, ratio);
        %         end
    end
    
    
    %% Figure 18: explore SNR dependencies
    if explorez
        explore_z(ana_data, roi_data, baseline,  recovery1)
    end
    
    if batchcall
        continue
    end
    
    %% reacalculate circular variance?
    if recalc_circ
        disp(' - - Recalculating circular variance')
        for y = 1:size(ana_data{i}.peaks,2)
            
            ana_data{i}.Fit(y).circvar_dir_ipsi = TT_CircularVariance(ana_data{i}.peaks(y).deltapeaks_averagetrace_ipsi');
            ana_data{i}.Fit(y).circvar_dir_contra =  TT_CircularVariance(ana_data{i}.peaks(y).deltapeaks_averagetrace_contra');
            [ana_data{i}.Fit(y).circvar_ipsi, ~, ana_data{i}.Fit(y).resultant_ori_angle_ipsi] = TT_CircularVariance_ORI(ana_data{i}.peaks(y).deltapeaks_averagetrace_ipsi');
            [ana_data{i}.Fit(y).circvar_contra, ~, ana_data{i}.Fit(y).resultant_ori_angle_contra] =  TT_CircularVariance_ORI(ana_data{i}.peaks(y).deltapeaks_averagetrace_contra');
            %                 if bino; ana_data{i}.Fit(y).circvar_bino = TT_CircularVariance(peaks(y).deltapeaks_averagetrace_bino');end
        end
    end
    
    
    
    %% SELECTION CTITERIA - CROSSECTIONAL  overlap (full, base + MD, n, n+1 etc.) - Full Session to session overlap with responsiveness criteria comes later!
    %
    %  SEE:
    %  http://www.evernote.com/shard/s4/sh/c5f66af5-c6ce-410e-9be3-1597c7bb21fa/586379acbd03d989985942d9b8257830
    %
    %  GROUPS! (repsonsiveness over timepoints is scored later)
    %  negative traces are always excluded.
    
    
    % RECALCULATE RESPONDER CRITERION
    z_thresh = 8; %8
    z_thresh2 = 3; %8
    z_thresh_fraction = 0.5; %0.5
    ipsi =1;
    contra = 1;
    and_or = 0;
    %     usemean = 1;
    
    zeropeak = 0;
    hillcorrect = 0;
    
    disp(' - - Extracting response selectors')
    [responder_new contra_crossers ipsi_crossers contra_resp ipsi_resp  contra_resp_low ipsi_resp_low]= response_selector(ana_data(i), roi_data(i), z_thresh, z_thresh_fraction, ipsi, contra, and_or);
    responder_new = [responder_new{:}]';
    contra_resp= [contra_resp{:}]';
    ipsi_resp= [ipsi_resp{:}]';
    contra_resp_low= [contra_resp_low{:}]';
    ipsi_resp_low= [ipsi_resp_low{:}]';
    
    roict = size([ana_data{i}.peaks(:).responder],2);
    %replace the first run data
    for rs = 1: roict
        %         zaverage_c = ana_data{i}.peaks(rs).deltapeaks_averagetrace_contra./ ana_data{i}.peaks(rs).sigma_averagetrace_contra;
        %         zaverage_i = ana_data{i}.peaks(rs).deltapeaks_averagetrace_ipsi./ ana_data{i}.peaks(rs).sigma_averagetrace_ipsi;
        %         ana_data{i}.peaks(rs).zscore_peaks_contra
        zover_c =  ana_data{i}.peaks(rs).zscore_peaks_contra > z_thresh2;
        zover_i =  ana_data{i}.peaks(rs).zscore_peaks_ipsi > z_thresh2;
        if zeropeak
            ana_data{i}.peaks(rs).deltapeaks_contra(~zover_c) = 0;
            ana_data{i}.peaks(rs).deltapeaks_ipsi(~zover_i) = 0;
            ana_data{i}.peaks(rs).maxAmpDelAve_contra = max(mean( ana_data{i}.peaks(rs).deltapeaks_contra,2));
            ana_data{i}.peaks(rs).maxAmpDelAve_ipsi =  max(mean( ana_data{i}.peaks(rs).deltapeaks_ipsi,2));
            %             ana_data{i}.peaks(rs).maxAmpDelAve_contra = max(ana_data{i}.peaks(rs).deltapeaks_averagetrace_contra);
            %             ana_data{i}.peaks(rs).maxAmpDelAve_ipsi = max(ana_data{i}.peaks(rs).deltapeaks_averagetrace_ipsi);
        end
        if hillcorrect
            Fsat = 800 * 100/85;
            P.n     = 2.9 %Hill coefficient
            P.k_d   = 1 %normalized .144e-6  %dissociation constant
            
            ana_data{i}.peaks(rs).deltapeaks_contra =  invHill_v2(P, ana_data{i}.peaks(rs).deltapeaks_contra ./ Fsat) .* Fsat;
            ana_data{i}.peaks(rs).deltapeaks_ipsi =  invHill_v2(P, ana_data{i}.peaks(rs).deltapeaks_ipsi ./ Fsat) .* Fsat;
            %             ana_data{i}.peaks(rs).deltapeaks_averagetrace_contra = invHill_v2(P, ana_data{i}.peaks(rs).deltapeaks_contra ./ Fsat) .* Fsat;
            %             ana_data{i}.peaks(rs).deltapeaks_averagetrace_ipsi = invHill_v2(P, ana_data{i}.peaks(rs).deltapeaks_ipsi./ Fsat) .* Fsat;
            ana_data{i}.peaks(rs).maxAmpDelAve_contra == max(mean( ana_data{i}.peaks(rs).deltapeaks_contra,2));
            ana_data{i}.peaks(rs).maxAmpDelAve_ipsi = max(mean( ana_data{i}.peaks(rs).deltapeaks_ipsi,2));
            
            
        end
        
        ana_data{i}.peaks(rs).responder = responder_new(rs);
        ana_data{i}.peaks(rs).contra_responder = contra_resp(rs);
        ana_data{i}.peaks(rs).ipsi_responder = ipsi_resp(rs);
        ana_data{i}.peaks(rs).contra_low_responder = contra_resp_low(rs);
        ana_data{i}.peaks(rs).ipsi_low_responder = ipsi_resp_low(rs);
    end
    %set empty groupcounts to 0
    eidx = cellfun('isempty',{roi_data{i}.ROIs.groupsize});
    [roi_data{i}.ROIs(eidx).groupsize] = deal(0);
    
    %ROI: Anova responsiveness score - don't use!
    %         fullidx_anov{i} = find([roi_data{i}.ROIs(:).groupsize] >= size(orig_exp,1) &  [ana_data{i}.peaks(:).responder_anov] ==1 & [ana_data{i}.peaks(:).excluded] == 0 );
    fullidx_anov{i} = find([roi_data{i}.ROIs(:).groupsize] >= size(orig_exp,1) &  [ana_data{i}.peaks(:).responder_anov] ==1 & [ana_data{i}.peaks(:).excluded] == 0 );
    
    %ROI: z-score responders
    fullidx_z_rois{i} = find([roi_data{i}.ROIs(:).groupsize] >= size(orig_exp,1) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
    
    %ROI: all responders, regardless of groups refound
    respidx_z_rois{i} = find([ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0 );
    
    %ROI: full groups - no responsiveness criterion, just rejection of
    %negative-going crazies
    fullidx_rois{i} = find([roi_data{i}.ROIs(:).groupsize] >= size(orig_exp,1) & [ana_data{i}.peaks(:).excluded] == 0);
    
    %  ROIS: tuned cells (ipsi or contra ): ALL RECS -> to get the
    %  population Delta ORI!
    p = 0.01;
    try
        tuned_both_rois{i} = find([roi_data{i}.ROIs(:).group]   &  [ana_data{i}.peaks(:).responder] ==1 &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
        tuned_both_groups{i} = [roi_data{i}.ROIs(tuned_both_rois{i}).group];
        FitIpsi = [ana_data{i}.Fit(tuned_both_rois{i}).ipsi];
        FitContra =[ana_data{i}.Fit(tuned_both_rois{i}).contra];
        tuned_both_PrefOriContra{:,i} = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        tuned_both_PrefOriIpsi{:,i} = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
        
        [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(tuned_both_PrefOriContra{:,i}*2, tuned_both_PrefOriIpsi{:,i}*2) ;
        tuned_both_PrefDeltaOri{:,i} = MinAngDiff /2;
        tuned_both_PrefDeltaOri_bidi{:,i} = AngDiff /2;
        tuned_both_PrefDeltaOri_bidi_offset{:,i} = median(tuned_both_PrefDeltaOri_bidi{:,i});
        
        tuned_both_bidi_offset_corr{:,i} = mod([tuned_both_PrefDeltaOri_bidi{:,i} + 90] -  tuned_both_PrefDeltaOri_bidi_offset{:,i},180)-90;
    catch
        tuned_both_PrefDeltaOri{:,i} = [];
        tuned_both_PrefDeltaOri_bidi{:,i} = [];
        tuned_both_PrefDeltaOri_bidi_offset{:,i} = [];
    end
    %GROUP: refound (responsive) groups of the chosen 4d baseline pair
    base_n_groups = intersect([roi_data{baseline_pair(1)}.ROIs(:).group], [roi_data{baseline_pair(2)}.ROIs(:).group]);
    base_n_groups(base_n_groups==0) = [];
    base_n_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
    base_n_morph_rois{i} =  find(ismember([roi_data{i}.ROIs(:).group], base_n_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
    
    %  ROIS: tuned cells (ipsi or contra ): baseline n n+1
    p = 0.01;
    
    base_n_tuned_contra_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n_groups) &  [ana_data{i}.peaks(:).responder] ==1&   [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p) ;
    base_n_tuned_contra_groups{i} = [roi_data{i}.ROIs(base_n_tuned_contra_rois{i}).group];
    
    base_n_tuned_ipsi_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
    base_n_tuned_ipsi_groups{i} = [roi_data{i}.ROIs(base_n_tuned_ipsi_rois{i}).group];
    
    base_n_tuned_both_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n_groups)  &  [ana_data{i}.peaks(:).responder] ==1 &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
    base_n_tuned_both_groups{i} = [roi_data{i}.ROIs(base_n_tuned_both_rois{i}).group];
    
    %GROUP: refound (responsive) groups of the chosen 14d baseline pair
    if ~isempty(baseline_pair14)
        base_n14_groups = intersect([roi_data{baseline_pair14(1)}.ROIs(:).group], [roi_data{baseline_pair14(2)}.ROIs(:).group]);
        base_n14_groups(base_n14_groups==0) = [];
        base_n14_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
        base_n14_morph_rois{i} =  find(ismember([roi_data{i}.ROIs(:).group], base_n14_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
        %  ROIS: tuned cells (ipsi or contra ): baseline n n+14
        p = 0.01;
        base_n14_tuned_contra_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p) ;
        base_n14_tuned_contra_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_contra_rois{i}).group];
        
        base_n14_tuned_ipsi_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
        base_n14_tuned_ipsi_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_ipsi_rois{i}).group];
        
        base_n14_tuned_both_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
        base_n14_tuned_both_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_ipsi_rois{i}).group];
        
        
        
        base_n14_AND_n4_groups =  intersect(base_n_groups, base_n14_groups);
        base_n14_AND_n4_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_AND_n4_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
        base_n14_AND_n4_morph_rois{i} =  find(ismember([roi_data{i}.ROIs(:).group], base_n14_AND_n4_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
        p = 0.01;
        base_n14_AND_n4_tuned_contra_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_AND_n4_groups)&  [ana_data{i}.peaks(:).responder] ==1 &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p) ;
        base_n14_AND_n4_tuned_contra_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_contra_rois{i}).group];
        
        base_n14_AND_n4_tuned_ipsi_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_AND_n4_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
        base_n14_AND_n4_tuned_ipsi_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_ipsi_rois{i}).group];
        
        base_n14_AND_n4_tuned_both_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n14_AND_n4_groups)&  [ana_data{i}.peaks(:).responder] ==1 &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <=p) ;
        base_n14_AND_n4_tuned_both_groups{i} = [roi_data{i}.ROIs(base_n14_tuned_ipsi_rois{i}).group];
    end
    
    %GROUP: refound (responsive) groups of 3 baseline sessions
    try
        try
            base_n3_groups = intersect([roi_data{baseline-2}.ROIs(:).group], [roi_data{baseline-1}.ROIs(:).group]);
            
        catch
            base_n3_groups = intersect([roi_data{baseline-1}.ROIs(:).group], [roi_data{baseline-1}.ROIs(:).group]);
        end
        base_n3_groups = intersect(base_n3_groups, [roi_data{baseline}.ROIs(:).group]);
        base_n3_groups(base_n3_groups==0) = [];
        base_n3_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n3_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
        base_n3_morph_rois{i} =  find(ismember([roi_data{i}.ROIs(:).group], base_n3_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
    end
    
    %% electors only necessary for MD experiments!
    try
        if l % selectors only necessary for MD experiments!
            %GROUP: refound (responsive) groups of 3 baseline sessions plus MD
            try
                base_n3_md_groups = intersect([roi_data{baseline-2}.ROIs(:).group], [roi_data{baseline-1}.ROIs(:).group]);
            catch
                base_md_n3_groups = intersect([roi_data{baseline-1}.ROIs(:).group], [roi_data{baseline-1}.ROIs(:).group]);
            end
            base_n3_md_groups = intersect(base_n3_md_groups, [roi_data{baseline}.ROIs(:).group]);
            base_n3_md_groups = intersect(base_n3_md_groups, [roi_data{baseline+1}.ROIs(:).group]);
            base_n3_md_groups(base_n3_md_groups==0) = [];
            base_n3_md_resp_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n3_md_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
            base_n3_md_morph_rois{i} =  find(ismember([roi_data{i}.ROIs(:).group], base_n3_md_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
            
            p = 0.01;
            try
                base_n3_md_tuned_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_n3_md_groups) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p & [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <= p & [ana_data{i}.peaks(:).excluded] == 0) ;
            catch
                base_n3_md_tuned_rois{i} = [];
            end
            %     base_n3_md_tuned_groups{i} = [roi_data{i}.ROIs(base_md_tuned_rois{i}).group];
            
            
            %GROUP: refound cells (NO RESP_CRIT) at last baseline and first MD.
            base_md_morph_groups = intersect([roi_data{baseline}.ROIs(:).group], [roi_data{baseline+1}.ROIs(:).group] );
            base_md_morph_groups(base_md_morph_groups==0) = [];
            %ROIs
            base_md_morph_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
            
            %GROUP: refound cells (NO RESP_CRIT) at last baseline and first MD.
            pre_base_md_morph_groups = intersect([roi_data{baseline}.ROIs(:).group], [roi_data{baseline+1}.ROIs(:).group] );
            pre_base_md_morph_groups = intersect(pre_base_md_morph_groups, [roi_data{baseline-1}.ROIs(:).group] );
            
            pre_base_md_morph_groups(pre_base_md_morph_groups==0) = [];
            %ROIs
            pre_base_md_morph_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], pre_base_md_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
            
            %GROUP refound cells (RESPONSIVE) at last baseline, pre-baseline and first MD.
            pre_base_md_resp_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], pre_base_md_morph_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
            pre_base_md_resp_groups{i} = [roi_data{i}.ROIs(pre_base_md_resp_rois{i}).group];
            pre_base_md_resp_groups(pre_base_md_resp_groups{i}==0) = [];
            
            %GROUP:  refound cells (NO RESP_CRIT)  baseline, pre-baseline and first MD.
            base_md_rec_morph_groups = intersect(base_md_morph_groups, [roi_data{baseline+recovery1}.ROIs(:).group] );
            base_md_rec_morph_groups(base_md_rec_morph_groups==0) = [];
            %ROIs
            base_md_rec_morph_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_rec_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
            base_md_rec_resp_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_rec_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0 &  [ana_data{i}.peaks(:).responder] ==1) ;
            
            try
                %GROUP:  refound cells (NO RESP_CRIT)  baseline first MD recovery MD2
                base_md_rec_md2_morph_groups = intersect(base_md_rec_morph_groups, [roi_data{baseline+baseline2+1}.ROIs(:).group] );
                base_md_rec_md2_morph_groups(base_md_rec_md2_morph_groups==0) = [];
                %ROIs
                base_md_rec_md2_morph_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_rec_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0) ;
                base_md_rec_md2_resp_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_rec_morph_groups) & [ana_data{i}.peaks(:).excluded] == 0 &  [ana_data{i}.peaks(:).responder] ==1 );
            end
            
            
            
            %GROUP refound cells (RESPONSIVE) at last baseline and first MD.
            base_md_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_morph_groups) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0) ;
            %ROIs
            base_md_z_groups{i} = [roi_data{i}.ROIs(base_md_z_rois{i}).group];
            base_md_z_groups(base_md_z_groups{i}==0) = [];
            
            %   ROIS: tuned cells (both eyes): baseline + MD
            p = 0.01;
            base_md_tuned_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_morph_groups)&  [ana_data{i}.peaks(:).responder] ==1 &  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p & [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <= p ) ;
            base_md_tuned_groups{i} = [roi_data{i}.ROIs(base_md_tuned_rois{i}).group];
            
            base_md_tuned_rois_v2{i} = find(ismember([roi_data{i}.ROIs(:).group], base_md_morph_groups) &  [ana_data{i}.peaks(:).Tune_Anova_orth_maxAmpDelAve_contra] <=p & [ana_data{i}.peaks(:).Tune_Anova_orth_maxAmpDelAve_ipsi] <= p ) ;
            base_md_tuned_groups_v2{i} = [roi_data{i}.ROIs(base_md_tuned_rois_v2{i}).group];
            
            % refound (responsive) group at last baseline and last recovery1
            base_rec_groups_t = intersect([roi_data{baseline}.ROIs(:).group], [roi_data{baseline+recovery1}.ROIs(:).group]);
            base_rec_groups_t(base_rec_groups_t==0) = [];
            base_rec_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_groups_t) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
            %           just groups
            base_rec_groups{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_groups_t) & [ana_data{i}.peaks(:).excluded] == 0);
            
            %   ROIS: tuned cells (both eyes): baseline + recovery
            p = 0.01;
            base_rec_tuned_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_groups{i}) &  [ana_data{i}.peaks(:).responder] ==1&  [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_contra] <=p & [ana_data{i}.peaks(:).Tune_Anova_maxAmpDelAve_ipsi] <= p ) ;
            base_rec_tuned_groups{i} = [roi_data{i}.ROIs(base_rec_tuned_rois{i}).group];
            
            base_rec_tuned_rois_v2{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_groups{i}) &  [ana_data{i}.peaks(:).Tune_Anova_orth_maxAmpDelAve_contra] <=p & [ana_data{i}.peaks(:).Tune_Anova_orth_maxAmpDelAve_ipsi] <= p ) ;
            base_rec_tuned_groups_v2{i} = [roi_data{i}.ROIs(base_rec_tuned_rois_v2{i}).group];
            
            % refound (responsive) group at last baseline and MDand last recovery1 and MD2
            base_rec_md_groups_a = intersect(base_rec_groups{i}, [roi_data{baseline+1}.ROIs(:).group]);
            base_rec_md_groups_a(base_rec_md_groups_a==0) = [];
            try
                base_rec_md_groups_t= intersect(base_rec_md_groups_a, [roi_data{baseline+baseline2+1}.ROIs(:).group]);
                base_rec_md_z_rois{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_md_groups_t) &  [ana_data{i}.peaks(:).responder] ==1 & [ana_data{i}.peaks(:).excluded] == 0);
                %           just groups
                base_rec_md_groups{i} = find(ismember([roi_data{i}.ROIs(:).group], base_rec_md_groups_t) &  [ana_data{i}.peaks(:).excluded] == 0);
            end
            %% Figure 1: CDFplot of all session-wise (crosssectional) responders
            
            % Figure 1: CDFplot of all session-wise (crosssectional) responders
            % (respidx_z). Regardless of chronic group. NO crazies excluded!
            if ~fitamp
                if ~noplot
                    if i == baseline || i == baseline + 1  || i == baseline + recovery1
                        h(cl) = cdfplot([ana_data{i}.peaks(respidx_z_rois{i}).ODscore_average_max]);hold all
                        k(cl) = vline(median([ana_data{i}.peaks(respidx_z_rois{i}).ODscore_average_max]));
                        set(h(cl), 'Color', colors(cl,:), 'LineWidth', 2);
                        set(k(cl), 'Color', colors(cl,:), 'LineWidth', 2);
                        if cl <3
                            cl = cl+1;
                        end
                    end
                end
            else
                if i == baseline || i == baseline + 1  || i == baseline + recovery1
                    FitIpsi = [ana_data{i}.Fit(respidx_z_rois{i}).ipsi];
                    FitContra = [ana_data{i}.Fit(respidx_z_rois{i}).contra];
                    
                    disp('using summed fit-amplitudes throughout as basis for OD calculations')
                    amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
                    amp_ipsi = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]);
                    ODplot= (amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
                    if ~noplot
                        h(cl) = cdfplot([ODplot]);hold all
                        k(cl) = vline(median([ODplot]))
                        set(h(cl), 'Color', colors(cl,:), 'LineWidth', 2);
                        set(k(cl), 'Color', colors(cl,:), 'LineWidth', 2);
                        if cl <3
                            cl = cl+1;
                        end
                    end
                end
            end
        end
    end
    if ~noplot
        title('population ODI (FIT!) - all ungrouped z-responders [respidx_z] (>8 Z-SCORE IN 50% OF THE TRIALS FOR AT LEAST ONE ORI)')
        legend('baseline', 'after MD', 'recovery');
    end
    %     try
    %% Data extraction: obtain POPULATION OD values (amplitudes averaged _before_ OD calculation)
    if md_incl %selectors only necessary for MD experiments!
        
        if ~fitamp
            %z-responders all timepoints refound`
            pop_amp_ODI(i) = ( mean([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_contra]) -  mean([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_ipsi]) ) /( mean([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_contra]) +  mean([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_ipsi]) );
            pop_amp_ODI_med(i) = ( median([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_contra]) -  median([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_ipsi]) ) /( median([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_contra]) +  median([ana_data{i}.peaks(fullidx_z_rois{i}).maxAmpDelAve_ipsi]) );
            
            %z-responders single TP (crosssecitonal responders - regardless of group)
            pop_amp_ODI_allresp(i) = ( mean([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_contra]) -  mean([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_ipsi]) ) /( mean([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_contra]) +  mean([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_ipsi]) );
            pop_amp_ODI_med_allresp(i) = ( median([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_contra]) -  median([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_ipsi]) ) /( median([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_contra]) +  median([ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_ipsi]) );
            
            %z-responders single TP (cross-secitonal responders)
            %groups bsaeline and MD
            OD_base_md_z_rois =  ( ([ana_data{i}.peaks(base_md_z_rois{i}).maxAmpDelAve_contra]) -  ([ana_data{i}.peaks(base_md_z_rois{i}).maxAmpDelAve_ipsi]) ) ./( ([ana_data{i}.peaks(base_md_z_rois{i}).maxAmpDelAve_contra]) +  ([ana_data{i}.peaks(base_md_z_rois{i}).maxAmpDelAve_ipsi]) );
            
            %all full groups, regardless of responsiveness (and same number of
            %cells per crtossection [-> negative rejects!)
            pop_amp_ODI_all_full(i) =   ( mean([ana_data{i}.peaks(fullidx_rois{i}).maxAmpDelAve_contra]) -  mean([ana_data{i}.peaks(fullidx_rois{i}).maxAmpDelAve_ipsi]) ) /( mean([ana_data{i}.peaks(fullidx_rois{i}).maxAmpDelAve_contra]) +  mean([ana_data{i}.peaks(fullidx_rois{i}).maxAmpDelAve_ipsi]) );;
        else
            %z-responders all timepoints refound (full morpho groups cross-sectional responders - not necessarily responsive throughout)
            FitIpsi = [ana_data{i}.Fit(fullidx_z_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(fullidx_z_rois{i}).contra];
            
            amp_contra = mean(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi =   mean(([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]));
            amp_contra_md = median(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi_md = median(([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]));
            
            pop_amp_ODI(i) =(amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            pop_amp_ODI_med(i) =  (amp_contra_md - amp_ipsi_md) ./ (amp_contra_md + amp_ipsi_md);
            
            %z-responders single TP (cross-secitonal responders - regardless of group)
            FitIpsi = [ana_data{i}.Fit(respidx_z_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(respidx_z_rois{i}).contra];
            
            amp_contra = mean(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi =   mean(([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]));
            amp_contra_md = median(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi_md = median(([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]));
            
            pop_amp_ODI_allresp(i) =(amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            pop_amp_ODI_med_allresp(i) = (amp_contra_md - amp_ipsi_md) ./ (amp_contra_md + amp_ipsi_md);
            
            %z-responders single TP (cross-secitonal responders)
            %groups bsaeline and MD
            
            FitIpsi = [ana_data{i}.Fit(base_md_z_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(base_md_z_rois{i}).contra];
            
            amp_contra = mean(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi =   mean(([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]));
            amp_contra_md = median(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi_md = median(([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]));
            
            amp_all_contra =  ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_all_ipsi = ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            OD_base_md_z_rois =  (amp_all_contra - amp_all_ipsi) ./ (amp_all_contra + amp_all_ipsi);
            
            
            %all full groups, regardless of responsiveness (and same number of
            %cells per crtossection [-> negative rejects!)
            FitIpsi = [ana_data{i}.Fit(fullidx_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(fullidx_rois{i}).contra];
            
            amp_contra = mean(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi =   mean(([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]));
            amp_contra_md = median(([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
            amp_ipsi_md = median(([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]));
            
            pop_amp_ODI_all_full(i) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            %         pop_amp_ODI_med_all_full(i) =  (amp_contra_md - amp_ipsi_md) ./ (amp_contra_md + amp_ipsi_md);
        end
        
        %% Data extraction: binary classification of OD types over 1st MD I
        
        
        ODclassifier = 0.25
        
        if i == baseline;
            Contra_resp_base_md_rois = find(OD_base_md_z_rois>ODclassifier);
            Contra_resp_base_md_goups= base_md_z_groups{i}(Contra_resp_base_md_rois);
            Contra_resp_base_md_goups(Contra_resp_base_md_goups == 0) = [];
            
            Bino_resp_base_md_rois = find(OD_base_md_z_rois<ODclassifier & OD_base_md_z_rois>-ODclassifier);
            Bino_resp_base_md_goups= base_md_z_groups{i}(Bino_resp_base_md_rois);
            Bino_resp_base_md_goups(Bino_resp_base_md_goups == 0) = [];
            
            Ipsi_BASE_resp_base_md_rois = find(OD_base_md_z_rois<-ODclassifier);
            Ipsi_BASE_resp_base_md_goups= base_md_z_groups{i}(Ipsi_BASE_resp_base_md_rois);
            Ipsi_BASE_resp_base_md_goups(Ipsi_BASE_resp_base_md_goups == 0)=[];
            
            all_BASE_resp_base_md_rois = find(OD_base_md_z_rois);
            all_BASE_resp_base_md_goups= base_md_z_groups{i}(all_BASE_resp_base_md_rois);
            all_BASE_resp_base_md_goups(all_BASE_resp_base_md_goups == 0)=[];
            
        elseif i == baseline +1;
            
            Contra_MD_resp_base_md_rois = find(OD_base_md_z_rois>ODclassifier);
            Contra_MD_resp_base_md_goups= base_md_z_groups{i}(Contra_MD_resp_base_md_rois);
            Contra_MD_resp_base_md_goups(Contra_MD_resp_base_md_goups == 0) = [];
            
            Ipsi_resp_base_md_rois = find(OD_base_md_z_rois<-ODclassifier);
            Ipsi_resp_base_md_goups= base_md_z_groups{i}(Ipsi_resp_base_md_rois); % IPSI POST HIOC (AFTE MD!_)
            Ipsi_resp_base_md_goups(Ipsi_resp_base_md_goups == 0)=[];
            
            Bino_MD_resp_base_md_rois = find(OD_base_md_z_rois<ODclassifier & OD_base_md_z_rois>-ODclassifier);
            Bino_MD_resp_base_md_goups= base_md_z_groups{i}(Bino_MD_resp_base_md_rois);
            Bino_MD_resp_base_md_goups(Bino_MD_resp_base_md_goups == 0) = [];
            
        end
    end
    %     end
    %% Data extraction: imaging timepoints (precise to the second)
    TB(i) = (ana_data{i}.triggertime - ana_data{1}.triggertime);
    
end
% THE MASTERLOOP ENDS HERE!
%
% if batchcall ||reanalyze
% %     return
% end

%% bootstrap single session response distribution based on sinlge trila stapping?

if  bootstrap_var
    runs = 500;
    tic
    disp(['bootsrapping tuning parameter variance from single session - ' num2str(runs) 'draws with replacement'])
    bs_ana_data = bootstrap_corrs(ana_data{baseline_pair(1)}, runs,1);
    toc
end

% try
if md_incl %selectors only necessary for MD experiments!
    
    %% Data extraction: binary classification of OD types over 1st MD II - ipsi contra MD survivors
    
    
    %bidirectional survivors
    total_cells = length(base_md_morph_groups);
    
    %pre-MD - > post-MD
    contra_base = length(find(ismember([roi_data{baseline}.ROIs(:).group], Contra_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1 ));
    contra_lost = length(find(ismember([roi_data{baseline+1}.ROIs(:).group], Contra_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==0) );
    contra_md = length(find(ismember([roi_data{baseline+1}.ROIs(:).group], Contra_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1 ));
    
    bino_base = length(find(ismember([roi_data{baseline}.ROIs(:).group], Bino_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1 ));
    bino_lost = length(find(ismember([roi_data{baseline+1}.ROIs(:).group], Bino_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==0) );
    bino_md = length(find(ismember([roi_data{baseline+1}.ROIs(:).group], Bino_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1 ));
    
    ipsi_base = length( find(ismember([roi_data{baseline}.ROIs(:).group], Ipsi_BASE_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1 ));
    ipsi_lost =length( find(ismember([roi_data{baseline+1}.ROIs(:).group], Ipsi_BASE_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==0));
    ipsi_md = length( find(ismember([roi_data{baseline+1}.ROIs(:).group], Ipsi_BASE_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1));
    
    %post-MD - > pre-MD
    ipsi_MD_md = length( find(ismember([roi_data{baseline+1}.ROIs(:).group], Ipsi_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1 ));
    ipsi_MD_gained =length( find(ismember([roi_data{baseline}.ROIs(:).group], Ipsi_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==0));
    ipsi_MD_base = length( find(ismember([roi_data{baseline}.ROIs(:).group], Ipsi_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1));
    
    contra_MD_md = length( find(ismember([roi_data{baseline+1}.ROIs(:).group], Contra_MD_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1 ));
    contra_MD_gained =length( find(ismember([roi_data{baseline}.ROIs(:).group], Contra_MD_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==0));
    contra_MD_base = length( find(ismember([roi_data{baseline}.ROIs(:).group], Contra_MD_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1));
    
    bino_MD_md = length( find(ismember([roi_data{baseline+1}.ROIs(:).group], Bino_MD_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1 ));
    bino_MD_gained =length( find(ismember([roi_data{baseline}.ROIs(:).group], Bino_MD_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==0));
    bino_MD_base = length( find(ismember([roi_data{baseline}.ROIs(:).group], Bino_MD_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] ==1));
    
    %EXPLANATION
    %read:
    %1. base: ODI group present & responsive pre-MD
    %2. md: how many of these are responsive post-md (REGARDLESS of ODI)
    %3. MD_md: ODI group present & responsive post-MD
    %4. MD_base: how many of these are responsive pre-md (REGARDLESS of ODI)
    
    contra_base_frac = contra_base/total_cells;
    contra_md_frac =  contra_md/total_cells;
    contra_MD_md_frac =  contra_MD_md/total_cells;
    contra_MD_base_frac =  contra_MD_base/total_cells;
    
    
    bino_base_frac = bino_base/total_cells;
    bino_md_frac =  bino_md/total_cells;
    bino_MD_md_frac =  bino_MD_md/total_cells;
    bino_MD_base_frac =  bino_MD_base/total_cells;
    
    ipsi_base_frac = ipsi_base/total_cells
    ipsi_md_frac = ipsi_md/total_cells
    ipsi_MD_md_frac =  ipsi_MD_md/total_cells;
    ipsi_MD_base_frac =  ipsi_MD_base/total_cells;
    
    %amplitudes of pre-MD responders pre and post md
    %bidirectional DeltaODI
    
    %contra (find last baseline responder in contra group and look at the same
    %cells after MD)
    
    % contra_base_md_rois{1} = find(ismember([roi_data{baseline}.ROIs(:).group], Contra_resp_base_md_goups));
    % contra_base_md_rois{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], Contra_resp_base_md_goups));
    %
    % [~, sorter{1}] = sort([roi_data{baseline}.ROIs(contra_base_md_rois{1}).group]);
    % [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(contra_base_md_rois{2}).group]);
    % disp_groups([contra_base_md_rois{1}(sorter{1}); contra_base_md_rois{2}(sorter{2})]', [3 4], [], mouse, exp ,roi_data, 25);
    
    contra_base_md_rois_resp{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], Contra_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1)
    [~,~, contra_base_md_rois_resp{1}] = intersect([roi_data{baseline+1}.ROIs(contra_base_md_rois_resp{2}).group], [roi_data{baseline}.ROIs(:).group]);
    contra_base_md_rois_resp{1} = contra_base_md_rois_resp{1}'
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(contra_base_md_rois_resp{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(contra_base_md_rois_resp{2}).group]);
    contra_base_md_rois_resp{1} = contra_base_md_rois_resp{1}(sorter{1});
    contra_base_md_rois_resp{2} = contra_base_md_rois_resp{2}(sorter{2});
    % disp_groups(contra_base_md_rois_resp, [baseline baseline+1], [], mouse, exp ,roi_data, 25);
    
    %ipsi (find last baseline responder in bino group and look at the same
    %cells after MD)
    
    % ipsi_base_md_rois{1} = find(ismember([roi_data{baseline}.ROIs(:).group], Ipsi_resp_base_md_goups));
    % ipsi_base_md_rois{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], Ipsi_resp_base_md_goups));
    
    ipsi_base_md_rois_resp{1} = find(ismember([roi_data{baseline}.ROIs(:).group], Ipsi_resp_base_md_goups) &  [ana_data{baseline}.peaks(:).responder] == 1)
    [~,~, ipsi_base_md_rois_resp{2}] = intersect([roi_data{baseline}.ROIs(ipsi_base_md_rois_resp{1}).group], [roi_data{baseline+1}.ROIs(:).group]);
    ipsi_base_md_rois_resp{2} = ipsi_base_md_rois_resp{2}'
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(ipsi_base_md_rois_resp{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(ipsi_base_md_rois_resp{2}).group]);
    ipsi_base_md_rois_resp{1} = ipsi_base_md_rois_resp{1}(sorter{1});
    ipsi_base_md_rois_resp{2} = ipsi_base_md_rois_resp{2}(sorter{2});
    
    % disp_groups(ipsi_base_md_rois_resp, [baseline baseline+1], [], mouse, exp ,roi_data, 25);
    
    %bino (find last baseline responder in bino group and look at the same
    %cells after MD)
    %
    % Bino_base_md_rois{1} = find(ismember([roi_data{baseline}.ROIs(:).group], Bino_resp_base_md_goups));
    % Bino_base_md_rois{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], Bino_resp_base_md_goups));
    
    bino_base_md_rois_resp{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], Bino_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1)
    [~,~, bino_base_md_rois_resp{1}] = intersect([roi_data{baseline+1}.ROIs(bino_base_md_rois_resp{2}).group], [roi_data{baseline}.ROIs(:).group]);
    bino_base_md_rois_resp{1} = bino_base_md_rois_resp{1}'
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(bino_base_md_rois_resp{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(bino_base_md_rois_resp{2}).group]);
    bino_base_md_rois_resp{1} = bino_base_md_rois_resp{1}(sorter{1});
    bino_base_md_rois_resp{2} = bino_base_md_rois_resp{2}(sorter{2});
    
    % disp_groups(bino_base_md_rois_resp, [baseline baseline+1], [], mouse, exp ,roi_data, 25);
    
    % all
    % all_base_md_rois_resp{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_z_groups{baseline}));
    % all_base_md_rois_resp{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_z_groups{baseline}));
    
    all_base_md_rois_resp{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], all_BASE_resp_base_md_goups) &  [ana_data{baseline+1}.peaks(:).responder] ==1)
    [~,~, all_base_md_rois_resp{1}] = intersect([roi_data{baseline+1}.ROIs(all_base_md_rois_resp{2}).group], [roi_data{baseline}.ROIs(:).group]);
    all_base_md_rois_resp{1} = all_base_md_rois_resp{1}'
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(all_base_md_rois_resp{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(all_base_md_rois_resp{2}).group]);
    all_base_md_rois_resp{1} = all_base_md_rois_resp{1}(sorter{1});
    all_base_md_rois_resp{2} = all_base_md_rois_resp{2}(sorter{2});
    % disp_groups(all_base_md_rois_resp, [baseline baseline+1], [], mouse, exp ,roi_data, 25,datapath);
    
    base_md_morph_groups_matched = intersect([roi_data{baseline}.ROIs(base_md_morph_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_morph_rois{baseline+1}).group]);
    base_md_morph_rois_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_morph_groups_matched)) ;
    base_md_morph_rois_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_morph_groups_matched)) ;
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_morph_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_morph_rois_matched{2}).group]);
    base_md_morph_rois_matched{1} = base_md_morph_rois_matched{1}(sorter{1});
    base_md_morph_rois_matched{2} = base_md_morph_rois_matched{2}(sorter{2});
    
    %pre baseline basliene md morph
    pre_base_md_morph_groups_matched = intersect([roi_data{baseline}.ROIs(pre_base_md_morph_rois{baseline}).group], [roi_data{baseline+1}.ROIs(pre_base_md_morph_rois{baseline+1}).group]);
    pre_base_md_morph_groups_matched = intersect(pre_base_md_morph_groups_matched, [roi_data{baseline-1}.ROIs(pre_base_md_morph_rois{baseline-1}).group]);
    pre_base_md_morph_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], pre_base_md_morph_groups_matched)) ;
    pre_base_md_morph_rois_matched{2} = find(ismember([roi_data{baseline}.ROIs(:).group], pre_base_md_morph_groups_matched)) ;
    pre_base_md_morph_rois_matched{3} = find(ismember([roi_data{baseline+1}.ROIs(:).group], pre_base_md_morph_groups_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(pre_base_md_morph_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline}.ROIs(pre_base_md_morph_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline+1}.ROIs(pre_base_md_morph_rois_matched{3}).group]);
    
    pre_base_md_morph_rois_matched{1} = pre_base_md_morph_rois_matched{1}(sorter{1});
    pre_base_md_morph_rois_matched{2} = pre_base_md_morph_rois_matched{2}(sorter{2});
    pre_base_md_morph_rois_matched{3} = pre_base_md_morph_rois_matched{3}(sorter{3});
    % disp_groups(pre_base_md_morph_rois_matched, [baseline-1 baseline baseline+1], [], mouse, exp ,roi_data, 25,datapath);
    
    %pre baseline basliene md responders
    pre_base_md_resp_groups_matched = intersect([roi_data{baseline}.ROIs(pre_base_md_resp_rois{baseline}).group], [roi_data{baseline+1}.ROIs(pre_base_md_resp_rois{baseline+1}).group]);
    pre_base_md_resp_groups_matched = intersect(pre_base_md_resp_groups_matched, [roi_data{baseline-1}.ROIs(pre_base_md_resp_rois{baseline-1}).group]);
    pre_base_md_resp_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], pre_base_md_resp_groups_matched)) ;
    pre_base_md_resp_rois_matched{2} = find(ismember([roi_data{baseline}.ROIs(:).group], pre_base_md_resp_groups_matched)) ;
    pre_base_md_resp_rois_matched{3} = find(ismember([roi_data{baseline+1}.ROIs(:).group], pre_base_md_resp_groups_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(pre_base_md_resp_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline}.ROIs(pre_base_md_resp_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline+1}.ROIs(pre_base_md_resp_rois_matched{3}).group]);
    
    pre_base_md_resp_rois_matched{1} = pre_base_md_resp_rois_matched{1}(sorter{1});
    pre_base_md_resp_rois_matched{2} = pre_base_md_resp_rois_matched{2}(sorter{2});
    pre_base_md_resp_rois_matched{3} = pre_base_md_resp_rois_matched{3}(sorter{3});
    % disp_groups(pre_base_md_resp_rois_matched, [baseline-1 baseline baseline+1], [], mouse, exp ,roi_data, 25,datapath);
    
    %tuned base+md base only tuned
    base_md_morph_groups_base_tuned_matched = [roi_data{baseline}.ROIs(base_md_tuned_rois{baseline}).group];
    base_md_morph_rois_base_tuned_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_morph_groups_base_tuned_matched)) ;
    base_md_morph_rois_base_tuned_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_morph_groups_base_tuned_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_morph_rois_base_tuned_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_morph_rois_base_tuned_matched{2}).group]);
    
    base_md_morph_rois_base_tuned_matched{1} = base_md_morph_rois_base_tuned_matched{1}(sorter{1});
    base_md_morph_rois_base_tuned_matched{2} = base_md_morph_rois_base_tuned_matched{2}(sorter{2});
    % disp_groups(base_md_morph_rois_base_tuned_matched,[baseline baseline+1], [], mouse, exp ,roi_data, 25, datapath);
    %  make_sortmap_from_PSTH(ana_data(baseline:baseline+1),base_md_morph_rois_base_tuned_matched, 1);
    %tuned base+md both
    base_md_morph_groups_tuned_matched = intersect([roi_data{baseline}.ROIs(base_md_tuned_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_tuned_rois{baseline+1}).group]);
    base_md_morph_rois_tuned_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_morph_groups_tuned_matched)) ;
    base_md_morph_rois_tuned_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_morph_groups_tuned_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_morph_rois_tuned_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_morph_rois_tuned_matched{2}).group]);
    
    base_md_morph_rois_tuned_matched{1} = base_md_morph_rois_tuned_matched{1}(sorter{1});
    base_md_morph_rois_tuned_matched{2} = base_md_morph_rois_tuned_matched{2}(sorter{2});
    
    % make_sortmap_from_PSTH(ana_data(baseline:baseline+1),base_md_morph_rois_matched, 1);
    % disp_groups(base_md_morph_rois_tuned_matched, [baseline baseline+1], [], mouse, exp ,roi_data, 25, datapath);
    
    %tuned base+rec both
    base_rec_morph_groups_tuned_matched = intersect([roi_data{baseline}.ROIs(base_rec_tuned_rois{baseline}).group], [roi_data{baseline+recovery1}.ROIs(base_rec_tuned_rois{baseline+recovery1}).group]);
    base_rec_morph_rois_tuned_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_rec_morph_groups_tuned_matched)) ;
    base_rec_morph_rois_tuned_matched{2} = find(ismember([roi_data{baseline+recovery1}.ROIs(:).group], base_rec_morph_groups_tuned_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_rec_morph_rois_tuned_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+recovery1}.ROIs(base_rec_morph_rois_tuned_matched{2}).group]);
    
    base_rec_morph_rois_tuned_matched{1} = base_rec_morph_rois_tuned_matched{1}(sorter{1});
    base_rec_morph_rois_tuned_matched{2} = base_rec_morph_rois_tuned_matched{2}(sorter{2});
    
    
    
    
    % disp_groups(all_base_md_rois_resp, [baseline baseline+1], [], mouse, exp ,roi_data, 25,datapath);
    
    % baseline md and recovery
    base_md_rec_morph_groups_matched_t = intersect([roi_data{baseline}.ROIs(base_md_rec_morph_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_rec_morph_rois{baseline+1}).group]);
    base_md_rec_morph_groups_matched = intersect(base_md_rec_morph_groups_matched_t, [roi_data{baseline+recovery1}.ROIs(base_md_rec_morph_rois{baseline+recovery1}).group]);
    
    
    base_md_rec_morph_rois_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_rec_morph_groups_matched)) ;
    base_md_rec_morph_rois_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_rec_morph_groups_matched)) ;
    base_md_rec_morph_rois_matched{3} = find(ismember([roi_data{baseline+recovery1}.ROIs(:).group], base_md_rec_morph_groups_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_rec_morph_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_rec_morph_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline+recovery1}.ROIs(base_md_rec_morph_rois_matched{3}).group]);
    
    
    base_md_rec_morph_rois_matched{1} = base_md_rec_morph_rois_matched{1}(sorter{1});
    base_md_rec_morph_rois_matched{2} = base_md_rec_morph_rois_matched{2}(sorter{2});
    base_md_rec_morph_rois_matched{3} = base_md_rec_morph_rois_matched{3}(sorter{3});
    
    % disp_groups(base_md_morph_rois_matched, [baseline baseline+1 baseline+recovery1], [], mouse, exp ,roi_data, 25, datapath);
    
    % baseline md and recovery responder
    base_md_rec_resp_groups_matched_t = intersect([roi_data{baseline}.ROIs(base_md_rec_resp_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_rec_resp_rois{baseline+1}).group]);
    base_md_rec_resp_groups_matched = intersect(base_md_rec_resp_groups_matched_t, [roi_data{baseline+recovery1}.ROIs(base_md_rec_resp_rois{baseline+recovery1}).group]);
    
    
    base_md_rec_resp_rois_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_rec_resp_groups_matched)) ;
    base_md_rec_resp_rois_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_rec_resp_groups_matched)) ;
    base_md_rec_resp_rois_matched{3} = find(ismember([roi_data{baseline+recovery1}.ROIs(:).group], base_md_rec_resp_groups_matched)) ;
    
    [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_rec_resp_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_rec_resp_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline+recovery1}.ROIs(base_md_rec_resp_rois_matched{3}).group]);
    
    
    base_md_rec_resp_rois_matched{1} = base_md_rec_resp_rois_matched{1}(sorter{1});
    base_md_rec_resp_rois_matched{2} = base_md_rec_resp_rois_matched{2}(sorter{2});
    base_md_rec_resp_rois_matched{3} = base_md_rec_resp_rois_matched{3}(sorter{3});
    
    %  disp_groups(base_md_rec_resp_rois_matched, [baseline baseline+1 baseline+recovery1], [], mouse, exp ,roi_data, 25, datapath);
    
    
    
    
    
    % baseline md and recovery and md2
    % try
    if baseline2 %BORKEN
        base_md_rec_md2_morph_groups_matched_t = intersect([roi_data{baseline}.ROIs(base_md_rec_md2_morph_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_rec_md2_morph_rois{baseline+1}).group]);
        base_md_rec_md2_morph_groups_matched = intersect(base_md_rec_md2_morph_groups_matched_t, [roi_data{baseline+recovery1}.ROIs(base_md_rec_md2_morph_rois{baseline+recovery1}).group]);
        base_md_rec_md2_morph_groups_matched = intersect(base_md_rec_md2_morph_groups_matched, [roi_data{baseline+baseline2}.ROIs(base_md_rec_md2_morph_rois{baseline+baseline2}).group]);
        base_md_rec_md2_morph_groups_matched = intersect(base_md_rec_md2_morph_groups_matched, [roi_data{baseline+baseline2+1}.ROIs(base_md_rec_md2_morph_rois{baseline+baseline2+1}).group]);
        
        base_md_rec_md2_morph_rois_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_rec_md2_morph_groups_matched)) ;
        base_md_rec_md2_morph_rois_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_rec_md2_morph_groups_matched)) ;
        base_md_rec_md2_morph_rois_matched{3} = find(ismember([roi_data{baseline+recovery1}.ROIs(:).group], base_md_rec_md2_morph_groups_matched)) ;
        base_md_rec_md2_morph_rois_matched{4} = find(ismember([roi_data{baseline+baseline2}.ROIs(:).group], base_md_rec_md2_morph_groups_matched)) ;
        base_md_rec_md2_morph_rois_matched{5} = find(ismember([roi_data{baseline+baseline2+1}.ROIs(:).group], base_md_rec_md2_morph_groups_matched)) ;
        
        [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_rec_md2_morph_rois_matched{1}).group]);
        [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_rec_md2_morph_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline+recovery1}.ROIs(base_md_rec_md2_morph_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+baseline2}.ROIs(base_md_rec_md2_morph_rois_matched{4}).group]);
        [~, sorter{5}] = sort([roi_data{baseline+baseline2+1}.ROIs(base_md_rec_md2_morph_rois_matched{5}).group]);
        
        base_md_rec_md2_morph_rois_matched{1} = base_md_rec_md2_morph_rois_matched{1}(sorter{1});
        base_md_rec_md2_morph_rois_matched{2} = base_md_rec_md2_morph_rois_matched{2}(sorter{2});
        base_md_rec_md2_morph_rois_matched{3} = base_md_rec_md2_morph_rois_matched{3}(sorter{3});
        base_md_rec_md2_morph_rois_matched{4} = base_md_rec_md2_morph_rois_matched{4}(sorter{4});
        base_md_rec_md2_morph_rois_matched{5} = base_md_rec_md2_morph_rois_matched{5}(sorter{5});
        
        %         disp_groups(base_md_rec_md2_morph_rois_matched, [baseline baseline+1 baseline+recovery1  baseline+baseline2+1], [], mouse, exp ,roi_data, 25, datapath);
        base_md_rec_md2_resp_groups_matched_t = intersect([roi_data{baseline}.ROIs(base_md_rec_md2_resp_rois{baseline}).group], [roi_data{baseline+1}.ROIs(base_md_rec_md2_resp_rois{baseline+1}).group]);
        base_md_rec_md2_resp_groups_matched = intersect(base_md_rec_md2_resp_groups_matched_t, [roi_data{baseline+recovery1}.ROIs(base_md_rec_md2_resp_rois{baseline+recovery1}).group]);
        base_md_rec_md2_resp_groups_matched = intersect(base_md_rec_md2_resp_groups_matched, [roi_data{baseline+baseline2}.ROIs(base_md_rec_md2_resp_rois{baseline+baseline2}).group]);
        base_md_rec_md2_resp_groups_matched = intersect(base_md_rec_md2_resp_groups_matched, [roi_data{baseline+baseline2+1}.ROIs(base_md_rec_md2_resp_rois{baseline+baseline2+1}).group]);
        
        base_md_rec_md2_resp_rois_matched{1} = find(ismember([roi_data{baseline}.ROIs(:).group], base_md_rec_md2_resp_groups_matched)) ;
        base_md_rec_md2_resp_rois_matched{2} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_md_rec_md2_resp_groups_matched)) ;
        base_md_rec_md2_resp_rois_matched{3} = find(ismember([roi_data{baseline+recovery1}.ROIs(:).group], base_md_rec_md2_resp_groups_matched)) ;
        base_md_rec_md2_resp_rois_matched{4} = find(ismember([roi_data{baseline+baseline2}.ROIs(:).group], base_md_rec_md2_resp_groups_matched)) ;
        base_md_rec_md2_resp_rois_matched{5} = find(ismember([roi_data{baseline+baseline2+1}.ROIs(:).group], base_md_rec_md2_resp_groups_matched)) ;
        
        [~, sorter{1}] = sort([roi_data{baseline}.ROIs(base_md_rec_md2_resp_rois_matched{1}).group]);
        [~, sorter{2}] = sort([roi_data{baseline+1}.ROIs(base_md_rec_md2_resp_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline+recovery1}.ROIs(base_md_rec_md2_resp_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+baseline2}.ROIs(base_md_rec_md2_resp_rois_matched{4}).group]);
        [~, sorter{5}] = sort([roi_data{baseline+baseline2+1}.ROIs(base_md_rec_md2_resp_rois_matched{5}).group]);
        
        base_md_rec_md2_resp_rois_matched{1} = base_md_rec_md2_resp_rois_matched{1}(sorter{1});
        base_md_rec_md2_resp_rois_matched{2} = base_md_rec_md2_resp_rois_matched{2}(sorter{2});
        base_md_rec_md2_resp_rois_matched{3} = base_md_rec_md2_resp_rois_matched{3}(sorter{3});
        base_md_rec_md2_resp_rois_matched{4} = base_md_rec_md2_resp_rois_matched{4}(sorter{4});
        base_md_rec_md2_resp_rois_matched{5} = base_md_rec_md2_resp_rois_matched{5}(sorter{5});
        
        %         disp_groups(base_md_rec_md2_resp_rois_matched, [baseline baseline+1 baseline+recovery1  baseline+baseline2+1], [], mouse, exp ,roi_data, 25, datapath);
        
        
    end
    % end
    
    
    % disp_groups(base_rec_resp_rois_tuned_matched, [baseline baseline+recovery1], [], mouse, exp ,roi_data, 25, datapath);
    % make_sortmap_from_PSTH(ana_data([baseline baseline+recovery1]),base_rec_resp_rois_tuned_matched, 1);
    
    % [base_md_tuned_rois
    
    % disp_groups(base_md_resp_rois_matched, [baseline baseline+1], [], mouse, exp ,roi_data, 25, datapath);
    % base_md_morph_rois(1) = intersect([roi_data{1}.ROIs(selector{1}).group], [roi_data{2}.ROIs(selector{2}).group]);
    %% Data extraction: extract pre-MD post-MD amplitudes, ODIs, ORI mismatches etc
    for ko = 1:2;
        if ~fitamp
            contra_base_md_ODI(:,ko) = ( [ana_data{baseline+ko-1}.peaks(contra_base_md_rois_resp{ko}).maxAmpDelAve_contra] -  [ana_data{baseline+ko-1}.peaks(contra_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] ) ./ ([ana_data{baseline+ko-1}.peaks(contra_base_md_rois_resp{ko}).maxAmpDelAve_contra] +  [ana_data{baseline+ko-1}.peaks(contra_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] );
            ipsi_base_md_ODI(:,ko) = ( [ana_data{baseline+ko-1}.peaks(ipsi_base_md_rois_resp{ko}).maxAmpDelAve_contra] -  [ana_data{baseline+ko-1}.peaks(ipsi_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] ) ./ ([ana_data{baseline+ko-1}.peaks(ipsi_base_md_rois_resp{ko}).maxAmpDelAve_contra] +  [ana_data{baseline+ko-1}.peaks(ipsi_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] );
            bino_base_md_ODI(:,ko) = ( [ana_data{baseline+ko-1}.peaks(bino_base_md_rois_resp{ko}).maxAmpDelAve_contra] -  [ana_data{baseline+ko-1}.peaks(bino_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] ) ./ ([ana_data{baseline+ko-1}.peaks(bino_base_md_rois_resp{ko}).maxAmpDelAve_contra] +  [ana_data{baseline+ko-1}.peaks(bino_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] );
            all_base_md_ODI(:,ko) = ( [ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_contra] -  [ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] ) ./ ([ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_contra] +  [ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_ipsi] );
            
            %all ODIs all ROIS
            base_md_amp_contra{ko} =  [ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_contra];
            base_md_amp_ipsi{ko} =  [ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}).maxAmpDelAve_ipsi];
            
            %all ODIs all ROIS
            base_md_circvar_contra{ko} =  [ana_data{baseline+ko-1}.Fit(all_base_md_rois_resp{ko}).circvar_contra];
            base_md_circvar_ipsi{ko} =  [ana_data{baseline+ko-1}.Fit(all_base_md_rois_resp{ko}).circvar_ipsi];
            
            
            %         for i = 1:
            
            for ii = 1:length(all_base_md_rois_resp{ko})
                base_md_OSI_contra_resp(ko,ii) = TT_OrientationSelectivityIndex(ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}(ii)).deltapeaks_averagetrace_contra');
                %         cv_ctest(ko,ii) = TT_CircularVariance(ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}(ii)).deltapeaks_averagetrace_contra');
                base_md_OSI_ipsi_resp(ko,ii) = TT_OrientationSelectivityIndex(ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}(ii)).deltapeaks_averagetrace_ipsi');
                %         cv_itest(ko,ii) = TT_CircularVariance(ana_data{baseline+ko-1}.peaks(all_base_md_rois_resp{ko}(ii)).deltapeaks_averagetrace_ipsi');
            end
            base_md_amp_contra_morph{ko} =  [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_amp_ipsi_morph{ko} =  [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            
            
            base_md_trials_contra{:,ko} =  [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).deltatrace_trials_oris_contra];
            base_md_trials_ipsi{:,ko} =    [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).deltatrace_trials_oris_ipsi];
            
            all_base_md_ODI_morph(:,ko) = ( base_md_amp_contra_morph{ko} - base_md_amp_ipsi_morph{ko} ) ./ (base_md_amp_contra_morph{ko} + base_md_amp_ipsi_morph{ko} );
            all_base_md_ODI_morph_nan(:,ko) = all_base_md_ODI_morph(:,ko);
            
            all_base_md_ODI_morph_nan([ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).responder]'==0,ko) = NaN;
            % ODIs and deltaORIs of all tuned morph ROIs
            
            try
                base_md_tunedODI_morph(:,ko) = ( [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_tuned_matched{ko}).maxAmpDelAve_contra] -  [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_tuned_matched{ko}).maxAmpDelAve_ipsi] ) ./ ([ana_data{baseline+ko-1}.peaks(base_md_morph_rois_tuned_matched{ko}).maxAmpDelAve_contra] +  [ana_data{baseline+ko-1}.peaks(base_md_morph_rois_tuned_matched{ko}).maxAmpDelAve_ipsi] );
                
                FitIpsi = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_tuned_matched{ko}).ipsi];
                FitContra = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_tuned_matched{ko}).contra];
                
                plotarray_base_md_tuned_morph_PrefOriContra(:,ko) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                plotarray_base_md_tuned_morph_PrefOriIpsi(:,ko) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_md_tuned_morph_PrefOriContra(:,ko)*2, plotarray_base_md_tuned_morph_PrefOriIpsi(:,ko)*2);
                plotarray_base_md_tuned_morph_DelatOri(:,ko) =   MinAngDiff /2;
            end
            
            
            
        else
            %contra
            FitIpsi = [ana_data{baseline+ko-1}.Fit(contra_base_md_rois_resp{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+ko-1}.Fit(contra_base_md_rois_resp{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            contra_base_md_ODI(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            
            %ipsi
            FitIpsi = [ana_data{baseline+ko-1}.Fit(ipsi_base_md_rois_resp{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+ko-1}.Fit(ipsi_base_md_rois_resp{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            ipsi_base_md_ODI(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            
            %bino
            FitIpsi = [ana_data{baseline+ko-1}.Fit(bino_base_md_rois_resp{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+ko-1}.Fit(bino_base_md_rois_resp{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            bino_base_md_ODI(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            
            %all ODIs, pre-MD timepoint responsive
            FitIpsi = [ana_data{baseline+ko-1}.Fit(all_base_md_rois_resp{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+ko-1}.Fit(all_base_md_rois_resp{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_base_md_ODI(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            
            base_md_amp_contra{ko} = amp_contra;
            base_md_amp_ipsi{ko} = amp_ipsi;
            
            %all ODIs all morph ROIs
            FitIpsi = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_matched{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_matched{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            base_md_amp_contra_morph{ko} = amp_contra;
            base_md_amp_ipsi_morph{ko} = amp_ipsi;
            
            all_base_md_ODI_morph(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_base_md_ODI_morph_nan(:,ko) = all_base_md_ODI_morph(:,ko);
            
            all_base_md_ODI_morph_nan([ana_data{baseline+ko-1}.peaks(base_md_morph_rois_matched{ko}).responder]'==0,ko) = NaN;
            
            
            
            % ODIs and deltaORIs of all morph ROIs
            try
                FitIpsi = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_tuned_matched{ko}).ipsi];
                FitContra = [ana_data{baseline+ko-1}.Fit(base_md_morph_rois_tuned_matched{ko}).contra];
                
                amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
                amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
                
                base_md_tunedODI_morph(:,ko) =  ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
                
                plotarray_base_md_tuned_morph_PrefOriContra(:,ko) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                plotarray_base_md_tuned_morph_PrefOriIpsi(:,ko) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_md_tuned_morph_PrefOriContra(:,ko)*2, plotarray_base_md_tuned_morph_PrefOriIpsi(:,ko)*2);
                plotarray_base_md_tuned_morph_DelatOri(:,ko) =   MinAngDiff /2;
            end
            
        end
        
    end
    
    %% data extraction baseline MD and Recovery
    
    if ~fitamp
        ko = 1;
        base_md_rec_amp_contra_morph{ko} =  [ana_data{baseline}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_morph{ko} =  [ana_data{baseline}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        base_md_rec_amp_contra_resp{ko} =  [ana_data{baseline}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_resp{ko} =  [ana_data{baseline}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        ko = 2;
        base_md_rec_amp_contra_morph{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_morph{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        base_md_rec_amp_contra_resp{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_resp{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        ko=3;
        base_md_rec_amp_contra_morph{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_morph{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        base_md_rec_amp_contra_resp{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_contra];
        base_md_rec_amp_ipsi_resp{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        
    end
    
    %% data extraction pre baseline, baseline and MD
    
    if ~fitamp
        ko = 1;
        pre_base_md_amp_contra_morph{ko} =  [ana_data{baseline-1}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_morph{ko} =  [ana_data{baseline-1}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        pre_base_md_amp_contra_resp{ko} =  [ana_data{baseline-1}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_resp{ko} =  [ana_data{baseline-1}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        ko = 2;
        pre_base_md_amp_contra_morph{ko} =  [ana_data{baseline}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_morph{ko} =  [ana_data{baseline}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        pre_base_md_amp_contra_resp{ko} =  [ana_data{baseline}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_resp{ko} =  [ana_data{baseline}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        ko=3;
        pre_base_md_amp_contra_morph{ko} =  [ana_data{baseline+1}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_morph{ko} =  [ana_data{baseline+1}.peaks(pre_base_md_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        pre_base_md_amp_contra_resp{ko} =  [ana_data{baseline+1}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_contra];
        pre_base_md_amp_ipsi_resp{ko} =  [ana_data{baseline+1}.peaks(pre_base_md_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        
    end
    
    
    %% data extraction baseline DM and Recovery and MD2
    
    
    if baseline2;
        if ~fitamp
            ko = 1;
            base_md_rec_md2_amp_contra_morph{ko} =  [ana_data{baseline}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_morph{ko} =  [ana_data{baseline}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko = 2;
            base_md_rec_md2_amp_contra_morph{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_morph{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko = 3;
            base_md_rec_md2_amp_contra_morph{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_morph{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko =4;
            base_md_rec_md2_amp_contra_morph{ko} =  [ana_data{baseline+baseline2}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_morph{ko} =  [ana_data{baseline+baseline2}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko =5;
            base_md_rec_md2_amp_contra_morph{ko} =  [ana_data{baseline+baseline2+1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_morph{ko} =  [ana_data{baseline+baseline2+1}.peaks(base_md_rec_md2_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko = 1;
            base_md_rec_md2_amp_contra_resp{ko} =  [ana_data{baseline}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_resp{ko} =  [ana_data{baseline}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko = 2;
            base_md_rec_md2_amp_contra_resp{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_resp{ko} =  [ana_data{baseline+1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko = 3;
            base_md_rec_md2_amp_contra_resp{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_resp{ko} =  [ana_data{baseline+recovery1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko =4;
            base_md_rec_md2_amp_contra_resp{ko} =  [ana_data{baseline+baseline2}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_resp{ko} =  [ana_data{baseline+baseline2}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
            ko =5;
            base_md_rec_md2_amp_contra_resp{ko} =  [ana_data{baseline+baseline2+1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_contra];
            base_md_rec_md2_amp_ipsi_resp{ko} =  [ana_data{baseline+baseline2+1}.peaks(base_md_rec_md2_resp_rois_matched{ko}).maxAmpDelAve_ipsi];
        end
    end
end

% end
%basline n n+1
n_n1_morph_groups_matched = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_morph_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_morph_rois{baseline_pair(2)}).group]);
n_n1_morph_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n1_morph_groups_matched)) ;
n_n1_morph_rois_matched{2} = find(ismember(unique([roi_data{baseline_pair(2)}.ROIs(:).group]), n_n1_morph_groups_matched)) ;
[~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n1_morph_rois_matched{1}).group]);
[~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n1_morph_rois_matched{2}).group]);
n_n1_morph_rois_matched{1} = n_n1_morph_rois_matched{1}(sorter{1});
n_n1_morph_rois_matched{2} = n_n1_morph_rois_matched{2}(sorter{2});

% disp_groups(n_n1_morph_rois_matched, [baseline_pair(1) baseline_pair(2)], [], mouse, exp ,roi_data, 25, datapath);
% try
%basline n n+1 n+2
try
    n_n2_morph_groups_matched = intersect([roi_data{baseline-2}.ROIs(base_n3_morph_rois{baseline-2}).group], [roi_data{baseline-1}.ROIs(base_n3_morph_rois{baseline-1}).group]);
catch
    n_n2_morph_groups_matched = intersect([roi_data{baseline-1}.ROIs(base_n3_morph_rois{baseline-1}).group], [roi_data{baseline-1}.ROIs(base_n3_morph_rois{baseline-1}).group]);
end
n_n2_morph_groups_matched = intersect(n_n2_morph_groups_matched, [roi_data{baseline}.ROIs(base_n3_morph_rois{baseline}).group]);

try
    
    n_n2_morph_rois_matched{1} = find(ismember([roi_data{baseline-2}.ROIs(:).group], n_n2_morph_groups_matched)) ;
catch
    n_n2_morph_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], n_n2_morph_groups_matched)) ;
end
n_n2_morph_rois_matched{2} = find(ismember([roi_data{baseline-1}.ROIs(:).group], n_n2_morph_groups_matched)) ;
n_n2_morph_rois_matched{3} = find(ismember([roi_data{baseline}.ROIs(:).group], n_n2_morph_groups_matched)) ;
try
    [~, sorter{1}] = sort([roi_data{baseline-2}.ROIs(n_n2_morph_rois_matched{1}).group]);
catch
    [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(n_n2_morph_rois_matched{1}).group]);
end
[~, sorter{2}] = sort([roi_data{baseline-1}.ROIs(n_n2_morph_rois_matched{2}).group]);
[~, sorter{3}] = sort([roi_data{baseline}.ROIs(n_n2_morph_rois_matched{3}).group]);
n_n2_morph_rois_matched{1} = n_n2_morph_rois_matched{1}(sorter{1});
n_n2_morph_rois_matched{2} = n_n2_morph_rois_matched{2}(sorter{2});
n_n2_morph_rois_matched{3} = n_n2_morph_rois_matched{3}(sorter{3});

% make_sortmap_from_PSTH(ana_data(baseline:baseline+1),base_md_morph_rois_matched, 1);


% disp_groups(n_n2_morph_rois_matched, [baseline-2:baseline], [], mouse, exp ,roi_data, 25, datapath);

try
    %     make_sortmap_from_PSTH(ana_data(baseline-2:baseline),n_n2_morph_rois_matched, 1); set(gcf, 'name', 'n_n2_morph_rois_matched');
catch
    %     make_sortmap_from_PSTH(ana_data(baseline-1:baseline),n_n2_morph_rois_matched(2:3), 1); set(gcf, 'name', 'n_n2_morph_rois_matched');
end
%
try
    if md_incl %selectors only necessary for MD experiments!
        %3 point baseline responsive plus MD (no resp crit for MD)
        try
            base_n3_md_resp_groups_matched = intersect([roi_data{baseline-2}.ROIs(base_n3_md_resp_rois{baseline-2}).group], [roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group]);
        catch
            base_n3_md_resp_groups_matched = intersect([roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group], [roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group]);
        end
        base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline}.ROIs(base_n3_md_resp_rois{baseline}).group]);
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_resp_rois{baseline+1}).group]); disp('base_n3_md_resp_groups_matched setting MD point to responsive!')
        base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs.group]); %MD timepoint does not have to be responsive!
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_resp_rois{baseline+1}).group]);
        %
        try
            base_n3_md_resp_rois_matched{1} = find(ismember([roi_data{baseline-2}.ROIs(:).group], base_n3_md_resp_groups_matched)) ;
        catch
            base_n3_md_resp_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_resp_groups_matched)) ;
        end
        base_n3_md_resp_rois_matched{2} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_resp_groups_matched)) ;
        base_n3_md_resp_rois_matched{3} = find(ismember([roi_data{baseline}.ROIs(:).group], base_n3_md_resp_groups_matched)) ;
        base_n3_md_resp_rois_matched{4} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_n3_md_resp_groups_matched)) ;
        
        try
            [~, sorter{1}] = sort([roi_data{baseline-2}.ROIs(base_n3_md_resp_rois_matched{1}).group]);
        catch
            [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_resp_rois_matched{1}).group]);
        end
        [~, sorter{2}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_resp_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline}.ROIs(base_n3_md_resp_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+1}.ROIs(base_n3_md_resp_rois_matched{4}).group]);
        
        base_n3_md_resp_rois_matched{1} = base_n3_md_resp_rois_matched{1}(sorter{1});
        base_n3_md_resp_rois_matched{2} = base_n3_md_resp_rois_matched{2}(sorter{2});
        base_n3_md_resp_rois_matched{3} = base_n3_md_resp_rois_matched{3}(sorter{3});
        base_n3_md_resp_rois_matched{4} = base_n3_md_resp_rois_matched{4}(sorter{4});
        
        
        % disp_groups(base_n3_md_resp_rois_matched, [baseline-2:baseline+1], [], mouse, exp ,roi_data, 25, datapath);
        % make_sortmap_from_PSTH(ana_data(baseline-2:baseline+1),base_n3_md_resp_rois_matched, 3);
        
        
        %3 point baseline responsive plus MD RESPONSIVE
        try
            base_n3_resp_md_resp_groups_matched = intersect([roi_data{baseline-2}.ROIs(base_n3_md_resp_rois{baseline-2}).group], [roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group]);
        catch
            base_n3_resp_md_resp_groups_matched = intersect([roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group], [roi_data{baseline-1}.ROIs(base_n3_md_resp_rois{baseline-1}).group]);
        end
        base_n3_resp_md_resp_groups_matched = intersect(base_n3_resp_md_resp_groups_matched, [roi_data{baseline}.ROIs(base_n3_md_resp_rois{baseline}).group]);
        base_n3_resp_md_resp_groups_matched = intersect(base_n3_resp_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_resp_rois{baseline+1}).group]); disp('base_n3_md_resp_resp_groups_matched setting MD point to responsive!')
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs.group]); %MD timepoint does not have to be responsive!
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_resp_rois{baseline+1}).group]);
        %
        try
            base_n3_resp_md_resp_rois_matched{1} = find(ismember([roi_data{baseline-2}.ROIs(:).group], base_n3_resp_md_resp_groups_matched)) ;
        catch
            base_n3_resp_md_resp_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_resp_md_resp_groups_matched)) ;
        end
        base_n3_resp_md_resp_rois_matched{2} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_resp_md_resp_groups_matched)) ;
        base_n3_resp_md_resp_rois_matched{3} = find(ismember([roi_data{baseline}.ROIs(:).group], base_n3_resp_md_resp_groups_matched)) ;
        base_n3_resp_md_resp_rois_matched{4} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_n3_resp_md_resp_groups_matched)) ;
        
        try
            [~, sorter{1}] = sort([roi_data{baseline-2}.ROIs(base_n3_resp_md_resp_rois_matched{1}).group]);
        catch
            [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(base_n3_resp_md_resp_rois_matched{1}).group]);
        end
        [~, sorter{2}] = sort([roi_data{baseline-1}.ROIs(base_n3_resp_md_resp_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline}.ROIs(base_n3_resp_md_resp_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+1}.ROIs(base_n3_resp_md_resp_rois_matched{4}).group]);
        
        base_n3_resp_md_resp_rois_matched{1} = base_n3_resp_md_resp_rois_matched{1}(sorter{1});
        base_n3_resp_md_resp_rois_matched{2} = base_n3_resp_md_resp_rois_matched{2}(sorter{2});
        base_n3_resp_md_resp_rois_matched{3} = base_n3_resp_md_resp_rois_matched{3}(sorter{3});
        base_n3_resp_md_resp_rois_matched{4} = base_n3_resp_md_resp_rois_matched{4}(sorter{4});
        
        
        % disp_groups(base_n3_resp_md_resp_rois_matched, [baseline-2:baseline+1], [], mouse, exp ,roi_data, 25, datapath);
        % make_sortmap_from_PSTH(ana_data(baseline-2:baseline+1),base_n3_resp_md_resp_rois_matched, 3);
        
        %3 point baseline MORPH plus MD
        try
            base_n3_md_morph_groups_matched = intersect([roi_data{baseline-2}.ROIs(base_n3_md_morph_rois{baseline-2}).group], [roi_data{baseline-1}.ROIs(base_n3_md_morph_rois{baseline-1}).group]);
        catch
            base_n3_md_morph_groups_matched = intersect([roi_data{baseline-1}.ROIs(base_n3_md_morph_rois{baseline-1}).group], [roi_data{baseline-1}.ROIs(base_n3_md_morph_rois{baseline-1}).group]);
        end
        base_n3_md_morph_groups_matched = intersect(base_n3_md_morph_groups_matched, [roi_data{baseline}.ROIs(base_n3_md_morph_rois{baseline}).group]);
        base_n3_md_morph_groups_matched = intersect(base_n3_md_morph_groups_matched, [roi_data{baseline+1}.ROIs.group]); %MD timepoint does not have to be responsive!
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_morph_rois{baseline+1}).group]);
        %
        try
            base_n3_md_morph_rois_matched{1} = find(ismember([roi_data{baseline-2}.ROIs(:).group], base_n3_md_morph_groups_matched)) ;
        catch
            base_n3_md_morph_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_morph_groups_matched)) ;
        end
        base_n3_md_morph_rois_matched{2} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_morph_groups_matched)) ;
        base_n3_md_morph_rois_matched{3} = find(ismember([roi_data{baseline}.ROIs(:).group], base_n3_md_morph_groups_matched)) ;
        base_n3_md_morph_rois_matched{4} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_n3_md_morph_groups_matched)) ;
        
        try
            [~, sorter{1}] = sort([roi_data{baseline-2}.ROIs(base_n3_md_morph_rois_matched{1}).group]);
        catch
            [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_morph_rois_matched{1}).group]);
        end
        [~, sorter{2}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_morph_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline}.ROIs(base_n3_md_morph_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+1}.ROIs(base_n3_md_morph_rois_matched{4}).group]);
        
        base_n3_md_morph_rois_matched{1} = base_n3_md_morph_rois_matched{1}(sorter{1});
        base_n3_md_morph_rois_matched{2} = base_n3_md_morph_rois_matched{2}(sorter{2});
        base_n3_md_morph_rois_matched{3} = base_n3_md_morph_rois_matched{3}(sorter{3});
        base_n3_md_morph_rois_matched{4} = base_n3_md_morph_rois_matched{4}(sorter{4});
        
        
        %3 point baseline TUNED MORPH plus MD
        try
            base_n3_md_morph_tuned_groups_matched = intersect([roi_data{baseline-2}.ROIs(base_n3_md_tuned_rois{baseline-2}).group], [roi_data{baseline-1}.ROIs(base_n3_md_tuned_rois{baseline-1}).group]);
        catch
            base_n3_md_morph_tuned_groups_matched = intersect([roi_data{baseline-1}.ROIs(base_n3_md_tuned_rois{baseline-1}).group], [roi_data{baseline-1}.ROIs(base_n3_md_tuned_rois{baseline-1}).group]);
        end
        base_n3_md_morph_tuned_groups_matched = intersect(base_n3_md_morph_tuned_groups_matched, [roi_data{baseline}.ROIs(base_n3_md_tuned_rois{baseline}).group]);
        base_n3_md_morph_tuned_groups_matched = intersect(base_n3_md_morph_tuned_groups_matched, [roi_data{baseline+1}.ROIs.group]); %MD timepoint does not have to be responsive!
        % base_n3_md_resp_groups_matched = intersect(base_n3_md_resp_groups_matched, [roi_data{baseline+1}.ROIs(base_n3_md_morph_rois{baseline+1}).group]);
        %
        try
            base_n3_md_morph_tuned_rois_matched{1} = find(ismember([roi_data{baseline-2}.ROIs(:).group], base_n3_md_morph_tuned_groups_matched)) ;
        catch
            base_n3_md_morph_tuned_rois_matched{1} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_morph_tuned_groups_matched)) ;
        end
        base_n3_md_morph_tuned_rois_matched{2} = find(ismember([roi_data{baseline-1}.ROIs(:).group], base_n3_md_morph_tuned_groups_matched)) ;
        base_n3_md_morph_tuned_rois_matched{3} = find(ismember([roi_data{baseline}.ROIs(:).group], base_n3_md_morph_tuned_groups_matched)) ;
        base_n3_md_morph_tuned_rois_matched{4} = find(ismember([roi_data{baseline+1}.ROIs(:).group], base_n3_md_morph_tuned_groups_matched)) ;
        
        try
            [~, sorter{1}] = sort([roi_data{baseline-2}.ROIs(base_n3_md_morph_tuned_rois_matched{1}).group]);
        catch
            [~, sorter{1}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_morph_tuned_rois_matched{1}).group]);
        end
        [~, sorter{2}] = sort([roi_data{baseline-1}.ROIs(base_n3_md_morph_tuned_rois_matched{2}).group]);
        [~, sorter{3}] = sort([roi_data{baseline}.ROIs(base_n3_md_morph_tuned_rois_matched{3}).group]);
        [~, sorter{4}] = sort([roi_data{baseline+1}.ROIs(base_n3_md_morph_tuned_rois_matched{4}).group]);
        
        base_n3_md_morph_tuned_rois_matched{1} = base_n3_md_morph_tuned_rois_matched{1}(sorter{1});
        base_n3_md_morph_tuned_rois_matched{2} = base_n3_md_morph_tuned_rois_matched{2}(sorter{2});
        base_n3_md_morph_tuned_rois_matched{3} = base_n3_md_morph_tuned_rois_matched{3}(sorter{3});
        base_n3_md_morph_tuned_rois_matched{4} = base_n3_md_morph_tuned_rois_matched{4}(sorter{4});
        
    end
end

% end

% disp_groups(base_n3_md_morph_tuned_rois_matched, [baseline-2:baseline+1], [], mouse, exp ,roi_data, 25, datapath);
% make_sortmap_from_PSTH(ana_data(baseline-2:baseline+1),base_n3_md_morph_tuned_rois_matched, 3);


%basline n n+1 contra tuned
n_n1_tuned_contra_groups_matched = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_tuned_contra_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_tuned_contra_rois{baseline_pair(2)}).group]);
n_n1_tuned_contra_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n1_tuned_contra_groups_matched)) ;
n_n1_tuned_contra_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n1_tuned_contra_groups_matched)) ;
[~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n1_tuned_contra_rois_matched{1}).group]);
[~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n1_tuned_contra_rois_matched{2}).group]);
n_n1_tuned_contra_rois_matched{1} = n_n1_tuned_contra_rois_matched{1}(sorter{1});
n_n1_tuned_contra_rois_matched{2} = n_n1_tuned_contra_rois_matched{2}(sorter{2});

% disp_groups(n_n1_tuned_contra_rois_matched, [baseline_pair(1) baseline_pair(2)], [], mouse, exp ,roi_data, 25, datapath);

%basline n n+1 ipsi tuned
n_n1_tuned_ipsi_groups_matched = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_tuned_ipsi_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_tuned_ipsi_rois{baseline_pair(2)}).group]);
n_n1_tuned_ipsi_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n1_tuned_ipsi_groups_matched)) ;
n_n1_tuned_ipsi_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n1_tuned_ipsi_groups_matched)) ;
[~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n1_tuned_ipsi_rois_matched{1}).group]);
[~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n1_tuned_ipsi_rois_matched{2}).group]);
n_n1_tuned_ipsi_rois_matched{1} = n_n1_tuned_ipsi_rois_matched{1}(sorter{1});
n_n1_tuned_ipsi_rois_matched{2} = n_n1_tuned_ipsi_rois_matched{2}(sorter{2});
% disp_groups(n_n1_tuned_ipsi_rois_matched, [baseline_pair(1) baseline_pair(2)], [], mouse, exp ,roi_data, 25, datapath);

%basline n n+1 both tuned
n_n1_tuned_both_groups_matched = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_tuned_both_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_tuned_both_rois{baseline_pair(2)}).group]);
n_n1_tuned_both_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n1_tuned_both_groups_matched)) ;
n_n1_tuned_both_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n1_tuned_both_groups_matched)) ;
[~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n1_tuned_both_rois_matched{1}).group]);
[~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n1_tuned_both_rois_matched{2}).group]);
n_n1_tuned_both_rois_matched{1} = n_n1_tuned_both_rois_matched{1}(sorter{1});
n_n1_tuned_both_rois_matched{2} = n_n1_tuned_both_rois_matched{2}(sorter{2});
% disp_groups(n_n1_tuned_ipsi_rois_matched, [baseline_pair(1) baseline_pair(2)], [], mouse, exp ,roi_data, 25, datapath);
% make_sortmap_from_PSTH(ana_data([baseline_pair(1) baseline_pair(2)]),n_n1_tuned_both_rois_matched, 1);

if ~isempty(baseline_pair14)
    
    
    %basline n n+1
    n_n14_morph_groups_matched = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_morph_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_morph_rois{baseline_pair14(2)}).group]);
    n_n14_morph_rois_matched{1} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n14_morph_groups_matched)) ;
    n_n14_morph_rois_matched{2} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n14_morph_groups_matched)) ;
    [~, sorter{1}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n14_morph_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n14_morph_rois_matched{2}).group]);
    n_n14_morph_rois_matched{1} = n_n14_morph_rois_matched{1}(sorter{1});
    n_n14_morph_rois_matched{2} = n_n14_morph_rois_matched{2}(sorter{2});
    %     disp_groups(n_n14_morph_rois_matched, [baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    
    %basline n n+1 14d contra tuned
    n_n1_14d_tuned_contra_groups_matched = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_tuned_contra_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_tuned_contra_rois{baseline_pair14(2)}).group]);
    n_n1_14d_tuned_contra_rois_matched{1} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n1_14d_tuned_contra_groups_matched)) ;
    n_n1_14d_tuned_contra_rois_matched{2} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n1_14d_tuned_contra_groups_matched)) ;
    [~, sorter{1}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n1_14d_tuned_contra_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n1_14d_tuned_contra_rois_matched{2}).group]);
    n_n1_14d_tuned_contra_rois_matched{1} = n_n1_14d_tuned_contra_rois_matched{1}(sorter{1});
    n_n1_14d_tuned_contra_rois_matched{2} = n_n1_14d_tuned_contra_rois_matched{2}(sorter{2});
    % disp_groups(n_n1_tuned_contra_rois_matched, [baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    % make_sortmap_from_PSTH(ana_data([baseline_pair14(1) baseline_pair14(2)]),n_n1_14d_tuned_contra_rois_matched, 1);
    
    %basline n n+1 14d ipsi tuned
    n_n1_14d_tuned_ipsi_groups_matched = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_tuned_ipsi_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_tuned_ipsi_rois{baseline_pair14(2)}).group]);
    n_n1_14d_tuned_ipsi_rois_matched{1} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n1_14d_tuned_ipsi_groups_matched)) ;
    n_n1_14d_tuned_ipsi_rois_matched{2} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n1_14d_tuned_ipsi_groups_matched)) ;
    [~, sorter{1}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n1_14d_tuned_ipsi_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n1_14d_tuned_ipsi_rois_matched{2}).group]);
    n_n1_14d_tuned_ipsi_rois_matched{1} = n_n1_14d_tuned_ipsi_rois_matched{1}(sorter{1});
    n_n1_14d_tuned_ipsi_rois_matched{2} = n_n1_14d_tuned_ipsi_rois_matched{2}(sorter{2});
    % disp_groups(n_n1_14d_tuned_ipsi_rois_matched, [baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    
    % make_sortmap_from_PSTH(ana_data([baseline_pair14(1) baseline_pair14(2)]),n_n1_14d_tuned_ipsi_rois_matched, 1);
    
    %basline n n+1 both tuned
    n_n1_14d_tuned_both_groups_matched = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_tuned_both_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_tuned_both_rois{baseline_pair14(2)}).group]);
    n_n1_14d_tuned_both_rois_matched{1} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n1_14d_tuned_both_groups_matched)) ;
    n_n1_14d_tuned_both_rois_matched{2} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n1_14d_tuned_both_groups_matched)) ;
    [~, sorter{1}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n1_14d_tuned_both_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n1_14d_tuned_both_rois_matched{2}).group]);
    n_n1_14d_tuned_both_rois_matched{1} = n_n1_14d_tuned_both_rois_matched{1}(sorter{1});
    n_n1_14d_tuned_both_rois_matched{2} = n_n1_14d_tuned_both_rois_matched{2}(sorter{2});
    % disp_groups(n_n1_tuned_both_rois_matched, [baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    % make_sortmap_from_PSTH(ana_data([baseline_pair14(1) baseline_pair14(2)]),n_n1_14d_tuned_both_rois_matched, 1);
    
    
    %4 point n4 AND n14 matched MORPH
    
    n_n14_AND_n4_morph_groups_matched = intersect(n_n1_morph_groups_matched, n_n14_morph_groups_matched);
    n_n14_AND_n4_morph_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n14_AND_n4_morph_groups_matched)) ;
    n_n14_AND_n4_morph_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n14_AND_n4_morph_groups_matched)) ;
    n_n14_AND_n4_morph_rois_matched{3} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n14_AND_n4_morph_groups_matched)) ;
    n_n14_AND_n4_morph_rois_matched{4} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n14_AND_n4_morph_groups_matched)) ;
    
    
    [~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n14_AND_n4_morph_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n14_AND_n4_morph_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n14_AND_n4_morph_rois_matched{3}).group]);
    [~, sorter{4}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n14_AND_n4_morph_rois_matched{4}).group]);
    
    n_n14_AND_n4_morph_rois_matched{1} = n_n14_AND_n4_morph_rois_matched{1}(sorter{1});
    n_n14_AND_n4_morph_rois_matched{2} = n_n14_AND_n4_morph_rois_matched{2}(sorter{2});
    n_n14_AND_n4_morph_rois_matched{3} = n_n14_AND_n4_morph_rois_matched{3}(sorter{3});
    n_n14_AND_n4_morph_rois_matched{4} = n_n14_AND_n4_morph_rois_matched{4}(sorter{4});
    
    
    %     disp_groups(n_n14_AND_n4_morph_rois_matched, [baseline_pair(1) baseline_pair(2) baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    
    %4 point n4 AND n14 matched FUNC
    
    n_n14_AND_n4_resp_groups_matched_t = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_AND_n4_z_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_AND_n4_z_rois{baseline_pair14(2)}).group]);
    n_n14_AND_n4_resp_groups_matched_t2 = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_z_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_z_rois{baseline_pair(2)}).group]);
    
    n_n14_AND_n4_resp_groups_matched = intersect(n_n14_AND_n4_resp_groups_matched_t, n_n14_AND_n4_resp_groups_matched_t2);
    n_n14_AND_n4_resp_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n14_AND_n4_resp_groups_matched)) ;
    n_n14_AND_n4_resp_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n14_AND_n4_resp_groups_matched)) ;
    n_n14_AND_n4_resp_rois_matched{3} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n14_AND_n4_resp_groups_matched)) ;
    n_n14_AND_n4_resp_rois_matched{4} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n14_AND_n4_resp_groups_matched)) ;
    
    
    [~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n14_AND_n4_resp_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n14_AND_n4_resp_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n14_AND_n4_resp_rois_matched{3}).group]);
    [~, sorter{4}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n14_AND_n4_resp_rois_matched{4}).group]);
    
    n_n14_AND_n4_resp_rois_matched{1} = n_n14_AND_n4_resp_rois_matched{1}(sorter{1});
    n_n14_AND_n4_resp_rois_matched{2} = n_n14_AND_n4_resp_rois_matched{2}(sorter{2});
    n_n14_AND_n4_resp_rois_matched{3} = n_n14_AND_n4_resp_rois_matched{3}(sorter{3});
    n_n14_AND_n4_resp_rois_matched{4} = n_n14_AND_n4_resp_rois_matched{4}(sorter{4});
    
    %      disp_groups(n_n14_AND_n4_resp_rois_matched, [baseline_pair(1) baseline_pair(2) baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    
    
    %basline n n+1 4d AND 14d contra tuned
    n_n14_AND_n4_tuned_contra_groups_matched = intersect(n_n1_tuned_contra_groups_matched, n_n1_14d_tuned_contra_groups_matched);
    n_n14_AND_n4_tuned_contra_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n14_AND_n4_tuned_contra_groups_matched)) ;
    n_n14_AND_n4_tuned_contra_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n14_AND_n4_tuned_contra_groups_matched)) ;
    n_n14_AND_n4_tuned_contra_rois_matched{3} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n14_AND_n4_tuned_contra_groups_matched)) ;
    n_n14_AND_n4_tuned_contra_rois_matched{4} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n14_AND_n4_tuned_contra_groups_matched)) ;
    
    
    [~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n14_AND_n4_tuned_contra_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n14_AND_n4_tuned_contra_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n14_AND_n4_tuned_contra_rois_matched{3}).group]);
    [~, sorter{4}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n14_AND_n4_tuned_contra_rois_matched{4}).group]);
    
    n_n14_AND_n4_tuned_contra_rois_matched{1} = n_n14_AND_n4_tuned_contra_rois_matched{1}(sorter{1});
    n_n14_AND_n4_tuned_contra_rois_matched{2} = n_n14_AND_n4_tuned_contra_rois_matched{2}(sorter{2});
    n_n14_AND_n4_tuned_contra_rois_matched{3} = n_n14_AND_n4_tuned_contra_rois_matched{3}(sorter{3});
    n_n14_AND_n4_tuned_contra_rois_matched{4} = n_n14_AND_n4_tuned_contra_rois_matched{4}(sorter{4});
    
    %      disp_groups(n_n14_AND_n4_resp_rois_matched, [baseline_pair(1) baseline_pair(2) baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    
    %basline n n+1 4d AND 14d contra tuned
    n_n14_AND_n4_tuned_ipsi_groups_matched = intersect(n_n1_tuned_ipsi_groups_matched, n_n1_14d_tuned_ipsi_groups_matched);
    n_n14_AND_n4_tuned_ipsi_rois_matched{1} = find(ismember([roi_data{baseline_pair(1)}.ROIs(:).group], n_n14_AND_n4_tuned_ipsi_groups_matched)) ;
    n_n14_AND_n4_tuned_ipsi_rois_matched{2} = find(ismember([roi_data{baseline_pair(2)}.ROIs(:).group], n_n14_AND_n4_tuned_ipsi_groups_matched)) ;
    n_n14_AND_n4_tuned_ipsi_rois_matched{3} = find(ismember([roi_data{baseline_pair14(1)}.ROIs(:).group], n_n14_AND_n4_tuned_ipsi_groups_matched)) ;
    n_n14_AND_n4_tuned_ipsi_rois_matched{4} = find(ismember([roi_data{baseline_pair14(2)}.ROIs(:).group], n_n14_AND_n4_tuned_ipsi_groups_matched)) ;
    
    
    [~, sorter{1}] = sort([roi_data{baseline_pair(1)}.ROIs(n_n14_AND_n4_tuned_ipsi_rois_matched{1}).group]);
    [~, sorter{2}] = sort([roi_data{baseline_pair(2)}.ROIs(n_n14_AND_n4_tuned_ipsi_rois_matched{2}).group]);
    [~, sorter{3}] = sort([roi_data{baseline_pair14(1)}.ROIs(n_n14_AND_n4_tuned_ipsi_rois_matched{3}).group]);
    [~, sorter{4}] = sort([roi_data{baseline_pair14(2)}.ROIs(n_n14_AND_n4_tuned_ipsi_rois_matched{4}).group]);
    
    n_n14_AND_n4_tuned_ipsi_rois_matched{1} = n_n14_AND_n4_tuned_ipsi_rois_matched{1}(sorter{1});
    n_n14_AND_n4_tuned_ipsi_rois_matched{2} = n_n14_AND_n4_tuned_ipsi_rois_matched{2}(sorter{2});
    n_n14_AND_n4_tuned_ipsi_rois_matched{3} = n_n14_AND_n4_tuned_ipsi_rois_matched{3}(sorter{3});
    n_n14_AND_n4_tuned_ipsi_rois_matched{4} = n_n14_AND_n4_tuned_ipsi_rois_matched{4}(sorter{4});
    %      disp_groups(n_n14_AND_n4_tuned_ipsi_rois_matched, [baseline_pair(1) baseline_pair(2) baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    % disp_groups(n_n14_AND_n4_tuned_ipsi_rois_matched, [baseline_pair14(1) baseline_pair14(2)], [], mouse, exp ,roi_data, 25, datapath);
    % make_sortmap_from_PSTH(ana_data([baseline_pair14(1) baseline_pair14(2)]),n_n1_14d_tuned_contra_rois_matched, 1);
    
end
%% Data extraction: BASELINE N N+1 amplitudes, population crosssectional MORPH
for ko = 1:2;
    if ~fitamp
        %all ODIs all ROIS
        n_n1_amp_contra_morph{ko} =  [ana_data{baseline_pair(ko)}.peaks(n_n1_morph_rois_matched{ko}).maxAmpDelAve_contra];
        n_n1_amp_ipsi_morph{ko} =  [ana_data{baseline_pair(ko)}.peaks(n_n1_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
        
        all_n_n1_ODI_morph(:,ko) = ( n_n1_amp_contra_morph{ko} - n_n1_amp_ipsi_morph{ko} ) ./ (n_n1_amp_contra_morph{ko} + n_n1_amp_ipsi_morph{ko} );
        
    else
        %all ODIs all ROIS
        FitIpsi = [ana_data{baseline_pair(ko)}.Fit(n_n1_morph_rois_matched{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
        FitContra = [ana_data{baseline_pair(ko)}.Fit(n_n1_morph_rois_matched{ko}).contra];
        
        amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
        amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
        
        all_n_n1_ODI_morph(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
        
        n_n1_amp_contra_morph{ko} = amp_contra;
        n_n1_amp_ipsi_morph{ko} = amp_ipsi;
    end
    
end

%% Data extraction: BASELINE N N+1 amplitudes, population crosssectional 12-14d MORPH
if ~isempty(baseline_pair14)
    for ko = 1:2;
        if ~fitamp
            %all ODIs all ROIS
            n_n14_amp_contra_morph{ko} =  [ana_data{baseline_pair14(ko)}.peaks(n_n14_morph_rois_matched{ko}).maxAmpDelAve_contra];
            n_n14_amp_ipsi_morph{ko} =  [ana_data{baseline_pair14(ko)}.peaks(n_n14_morph_rois_matched{ko}).maxAmpDelAve_ipsi];
            
            all_n_n14_ODI_morph(:,ko) = ( n_n14_amp_contra_morph{ko} - n_n14_amp_ipsi_morph{ko} ) ./ (n_n14_amp_contra_morph{ko} + n_n14_amp_ipsi_morph{ko} );
            
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline_pair14(ko)}.Fit(n_n14_morph_rois_matched{ko}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline_pair14(ko)}.Fit(n_n14_morph_rois_matched{ko}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_n_n14_ODI_morph(:,ko) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            
            n_n14_amp_contra_morph{ko} = amp_contra;
            n_n14_amp_ipsi_morph{ko} = amp_ipsi;
        end
        
    end
end

%% Data extraction: 4 BASELINE  n4 and n14 pairs amplitudes, MORPH & Func & ORI tuned

if ~isempty(baseline_pair14)
    if ~fitamp
        n_n1_AND_n14_amp_contra_morph{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_morph_rois_matched{1}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_morph{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_morph_rois_matched{1}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_morph{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_morph_rois_matched{2}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_morph{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_morph_rois_matched{2}).maxAmpDelAve_ipsi];
        
        n_n1_AND_n14_amp_contra_morph{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_morph_rois_matched{3}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_morph{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_morph_rois_matched{3}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_morph{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_morph_rois_matched{4}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_morph{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_morph_rois_matched{4}).maxAmpDelAve_ipsi];
        
        
        n_n1_AND_n14_amp_contra_resp{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_resp_rois_matched{1}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_resp{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_resp_rois_matched{1}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_resp{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_resp_rois_matched{2}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_resp{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_resp_rois_matched{2}).maxAmpDelAve_ipsi];
        
        n_n1_AND_n14_amp_contra_resp{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_resp_rois_matched{3}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_resp{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_resp_rois_matched{3}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_resp{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_resp_rois_matched{4}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_resp{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_resp_rois_matched{4}).maxAmpDelAve_ipsi];
        
        
        n_n1_AND_n14_circvar_contra_resp{1} = [ana_data{baseline_pair(1)}.Fit(n_n14_AND_n4_resp_rois_matched{1}).circvar_contra];
        n_n1_AND_n14_circvar_ipsi_resp{1} = [ana_data{baseline_pair(1)}.Fit(n_n14_AND_n4_resp_rois_matched{1}).circvar_ipsi];
        n_n1_AND_n14_circvar_contra_resp{2} = [ana_data{baseline_pair(2)}.Fit(n_n14_AND_n4_resp_rois_matched{2}).circvar_contra];
        n_n1_AND_n14_circvar_ipsi_resp{2} = [ana_data{baseline_pair(2)}.Fit(n_n14_AND_n4_resp_rois_matched{2}).circvar_ipsi];
        
        n_n1_AND_n14_circvar_contra_resp{3} = [ana_data{baseline_pair14(1)}.Fit(n_n14_AND_n4_resp_rois_matched{3}).circvar_contra];
        n_n1_AND_n14_circvar_ipsi_resp{3} = [ana_data{baseline_pair14(1)}.Fit(n_n14_AND_n4_resp_rois_matched{3}).circvar_ipsi];
        n_n1_AND_n14_circvar_contra_resp{4} = [ana_data{baseline_pair14(2)}.Fit(n_n14_AND_n4_resp_rois_matched{4}).circvar_contra];
        n_n1_AND_n14_circvar_ipsi_resp{4} = [ana_data{baseline_pair14(2)}.Fit(n_n14_AND_n4_resp_rois_matched{4}).circvar_ipsi];
        
        
        n_n1_AND_n14_plotarray_Contra_resp{1}= [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_resp_rois_matched{1}).contra_responder];
        n_n1_AND_n14_plotarray_Ipsi_resp{1}  = [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_resp_rois_matched{1}).ipsi_responder];
        n_n1_AND_n14_plotarray_Contra_resp{2}= [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_resp_rois_matched{2}).contra_responder];
        n_n1_AND_n14_plotarray_Ipsi_resp{2}  = [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_resp_rois_matched{2}).ipsi_responder];
        
        n_n1_AND_n14_plotarray_Contra_resp{3}= [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_resp_rois_matched{3}).contra_responder];
        n_n1_AND_n14_plotarray_Ipsi_resp{3}  = [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_resp_rois_matched{3}).ipsi_responder];
        n_n1_AND_n14_plotarray_Contra_resp{4}= [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_resp_rois_matched{4}).contra_responder];
        n_n1_AND_n14_plotarray_Ipsi_resp{4}  = [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_resp_rois_matched{4}).ipsi_responder];
        
        for pop = 1:4;
            n_n1_AND_n14_ODI_resp{pop} = (n_n1_AND_n14_amp_contra_resp{pop} - n_n1_AND_n14_amp_ipsi_resp{pop})./(n_n1_AND_n14_amp_contra_resp{pop} + n_n1_AND_n14_amp_ipsi_resp{pop});
        end
        
        n_n1_AND_n14_amp_contra_tuned{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_tuned_contra_rois_matched{1}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_tuned{1} =  [ana_data{baseline_pair(1)}.peaks(n_n14_AND_n4_tuned_ipsi_rois_matched{1}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_tuned{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_tuned_contra_rois_matched{2}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_tuned{2} =  [ana_data{baseline_pair(2)}.peaks(n_n14_AND_n4_tuned_ipsi_rois_matched{2}).maxAmpDelAve_ipsi];
        
        n_n1_AND_n14_amp_contra_tuned{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_tuned_contra_rois_matched{3}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_tuned{3} =  [ana_data{baseline_pair14(1)}.peaks(n_n14_AND_n4_tuned_ipsi_rois_matched{3}).maxAmpDelAve_ipsi];
        n_n1_AND_n14_amp_contra_tuned{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_tuned_contra_rois_matched{4}).maxAmpDelAve_contra];
        n_n1_AND_n14_amp_ipsi_tuned{4} =  [ana_data{baseline_pair14(2)}.peaks(n_n14_AND_n4_tuned_ipsi_rois_matched{4}).maxAmpDelAve_ipsi];
        
    else
        %not implemented
    end
    try
        FitIpsi = [ana_data{baseline_pair(1)}.Fit(n_n14_AND_n4_tuned_ipsi_rois_matched{1}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
        FitContra = [ana_data{baseline_pair(1)}.Fit(n_n14_AND_n4_tuned_contra_rois_matched{1}).contra];
        plotarray_n_n1_AND_n14_PrefOriContra(:,1) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        plotarray_n_n1_AND_n14_PrefOriIpsi(:,1) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
        
        FitIpsi = [ana_data{baseline_pair(2)}.Fit(n_n14_AND_n4_tuned_ipsi_rois_matched{2}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
        FitContra = [ana_data{baseline_pair(2)}.Fit(n_n14_AND_n4_tuned_contra_rois_matched{2}).contra];
        plotarray_n_n1_AND_n14_PrefOriContra(:,2) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        plotarray_n_n1_AND_n14_PrefOriIpsi(:,2) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
        
        FitIpsi = [ana_data{baseline_pair14(1)}.Fit(n_n14_AND_n4_tuned_ipsi_rois_matched{3}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
        FitContra = [ana_data{baseline_pair14(1)}.Fit(n_n14_AND_n4_tuned_contra_rois_matched{3}).contra];
        plotarray_n_n1_AND_n14_PrefOriContra(:,3) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        plotarray_n_n1_AND_n14_PrefOriIpsi(:,3) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
        
        FitIpsi = [ana_data{baseline_pair14(2)}.Fit(n_n14_AND_n4_tuned_ipsi_rois_matched{4}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
        FitContra = [ana_data{baseline_pair14(2)}.Fit(n_n14_AND_n4_tuned_contra_rois_matched{4}).contra];
        plotarray_n_n1_AND_n14_PrefOriContra(:,4) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        plotarray_n_n1_AND_n14_PrefOriIpsi(:,4) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
    end
end

%% Data extraction: 3 BASELINE  amplitudes, population crosssectional
try
    pp=1;
    for ko = 2:-1:0;
        
        if ~fitamp
            %all ODIs all ROIS
            n_n2_amp_contra_morph{pp} =  [ana_data{baseline-ko}.peaks(n_n2_morph_rois_matched{pp}).maxAmpDelAve_contra];
            n_n2_amp_ipsi_morph{pp} =  [ana_data{baseline-ko}.peaks(n_n2_morph_rois_matched{pp}).maxAmpDelAve_ipsi];
            all_n_n2_ODI_morph(:,pp) = ( n_n2_amp_contra_morph{pp} - n_n2_amp_ipsi_morph{pp} ) ./ (n_n2_amp_contra_morph{pp} + n_n2_amp_ipsi_morph{pp} );
            all_n_n2_ODI_morph_nan(:,pp)  = all_n_n2_ODI_morph(:,pp) ;
            all_n_n2_ODI_morph_nan([ana_data{baseline-ko}.peaks(n_n2_morph_rois_matched{pp}).responder]'==0,pp) = NaN;
            pp = pp +1;
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline-ko}.Fit(n_n2_morph_rois_matched{pp}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline-ko}.Fit(n_n2_morph_rois_matched{pp}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_n_n2_ODI_morph(:,pp) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_n_n2_ODI_morph_nan(:,pp)  = all_n_n2_ODI_morph(:,pp) ;
            all_n_n2_ODI_morph_nan([ana_data{baseline-ko}.peaks(n_n2_morph_rois_matched{pp}).responder]'==0,pp) = NaN;
            
            n_n2_amp_contra_morph{pp} = amp_contra;
            n_n2_amp_ipsi_morph{pp} = amp_ipsi;
            pp = pp +1;
        end
        
    end
end
% try

if md_incl %selectors only necessary for MD experiments!
    %% Data extraction: 3 BASELINE  RESPONSIVE MD no resp crit amplitudes, population crosssectional
    pp=1;
    for ko = 3:-1:0;
        
        if ~fitamp
            %all ODIs all ROIS
            base_n3_md_amp_contra{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}).maxAmpDelAve_contra];
            base_n3_md_amp_ipsi{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}).maxAmpDelAve_ipsi];
            all_base_n3_md_ODI(:,pp) = ( base_n3_md_amp_contra{pp} - base_n3_md_amp_ipsi{pp} ) ./ (base_n3_md_amp_contra{pp} + base_n3_md_amp_ipsi{pp} );
            all_base_n3_md_ODI_nan(:,pp)  = all_base_n3_md_ODI(:,pp) ;
            all_base_n3_md_ODI_nan([ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}).responder]'==0,pp) = NaN;
            for ii = 1:length(base_n3_md_resp_rois_matched{1})
                all_base_n3_md_OSI_contra(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}(ii)).deltapeaks_averagetrace_contra');
                all_base_n3_md_OSI_ipsi(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}(ii)).deltapeaks_averagetrace_ipsi');
            end
            
            all_base_n3_md_circvar_contra(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_md_resp_rois_matched{pp}).circvar_contra];
            all_base_n3_md_circvar_ipsi(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_md_resp_rois_matched{pp}).circvar_ipsi];
            
            
            all_base_n3_md_trials_contra{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}).deltatrace_trials_oris_contra];
            all_base_n3_md_trials_ipsi{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_resp_rois_matched{pp}).deltatrace_trials_oris_ipsi];
            
            
            pp = pp +1;
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline+1-ko}.Fit(base_n3_md_resp_rois_matched{pp}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+1-ko}.Fit(base_n3_md_resp_rois_matched{pp}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_base_n3_md_ODI(:,pp) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_base_n3_md_ODI_nan(:,pp)  = all_n_n2_ODI_morph(:,pp) ;
            all_base_n3_md_ODI_nan([ana_data{baseline-ko}.peaks(base_n3_md_resp_rois_matched{pp}).responder]'==0,pp) = NaN;
            
            base_n3_md_amp_contra{pp} = amp_contra;
            base_n3_md_amp_ipsi{pp} = amp_ipsi;
            pp = pp +1;
        end
        
    end
    %% Data extraction: 3 BASELINE AND MD  RESPONSIVE amplitudes, population crosssectional
    
    pp=1;
    for ko = 3:-1:0;
        
        if ~fitamp
            %all ODIs all ROIS
            base_n3_resp_md_amp_contra{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).maxAmpDelAve_contra];
            base_n3_resp_md_amp_ipsi{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).maxAmpDelAve_ipsi];
            all_base_n3_resp_md_ODI(:,pp) = ( base_n3_resp_md_amp_contra{pp} - base_n3_resp_md_amp_ipsi{pp} ) ./ (base_n3_resp_md_amp_contra{pp} + base_n3_resp_md_amp_ipsi{pp} );
            all_base_n3_resp_resp_md_ODI_nan(:,pp)  = all_base_n3_resp_md_ODI(:,pp) ;
            all_base_n3_resp_resp_md_ODI_nan([ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).responder]'==0,pp) = NaN;
            for ii = 1:length(base_n3_resp_md_resp_rois_matched{1})
                all_base_n3_resp_md_OSI_contra(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}(ii)).deltapeaks_averagetrace_contra');
                all_base_n3_resp_md_OSI_ipsi(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}(ii)).deltapeaks_averagetrace_ipsi');
            end
            
            all_base_n3_resp_md_circvar_contra(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_resp_md_resp_rois_matched{pp}).circvar_contra];
            all_base_n3_resp_md_circvar_ipsi(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_resp_md_resp_rois_matched{pp}).circvar_ipsi];
            
            all_base_n3_resp_md_trials_contra{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).deltatrace_trials_oris_contra];
            all_base_n3_resp_md_trials_ipsi{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).deltatrace_trials_oris_ipsi];
            
            
            all_base_n3_resp_md_Contra_resp{:,pp}= [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).contra_responder];
            all_base_n3_resp_md_Ipsi_resp{:,pp}  = [ana_data{baseline+1-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).ipsi_responder];
            
            
            
            pp = pp +1;
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline+1-ko}.Fit(base_n3_resp_md_resp_rois_matched{pp}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+1-ko}.Fit(base_n3_resp_md_resp_rois_matched{pp}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_base_n3_resp_md_ODI(:,pp) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_base_n3_resp_resp_md_ODI_nan(:,pp)  = all_n_n2_ODI_morph(:,pp) ;
            all_base_n3_resp_resp_md_ODI_nan([ana_data{baseline-ko}.peaks(base_n3_resp_md_resp_rois_matched{pp}).responder]'==0,pp) = NaN;
            
            base_n3_resp_md_amp_contra{pp} = amp_contra;
            base_n3_resp_md_amp_ipsi{pp} = amp_ipsi;
            pp = pp +1;
        end
    end
    
    
    %% Data extraction: 3 BASELINE AND MD  MORPH amplitudes, population crosssectional
    pp=1;
    for ko = 3:-1:0;
        
        if ~fitamp
            %all ODIs all ROIS
            base_n3_md_amp_contra_morph{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).maxAmpDelAve_contra];
            base_n3_md_amp_ipsi_morph{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).maxAmpDelAve_ipsi];
            all_base_n3_md_ODI_morph(:,pp) = ( base_n3_md_amp_contra_morph{pp} - base_n3_md_amp_ipsi_morph{pp} ) ./ (base_n3_md_amp_contra_morph{pp} + base_n3_md_amp_ipsi_morph{pp} );
            all_base_n3_md_ODI_nan_morph(:,pp)  = all_base_n3_md_ODI_morph(:,pp) ;
            all_base_n3_md_ODI_nan_morph([ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).responder]'==0,pp) = NaN;
            base_n3_md_morph_respidx{pp} = [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).responder];
            
            all_base_n3_md_trials_contra_morph{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).deltatrace_trials_oris_contra];
            all_base_n3_md_trials_ipsi_morph{:,pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_rois_matched{pp}).deltatrace_trials_oris_ipsi];
            
            
            pp = pp +1;
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_rois_matched{pp}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_rois_matched{pp}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_base_n3_md_ODI_morph(:,pp) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_base_n3_md_ODI_nan_morph(:,pp)  = all_base_n3_md_ODI_morph(:,pp) ;
            all_base_n3_md_ODI_nan_morph([ana_data{baseline-ko}.peaks(base_n3_md_morph_rois_matched{pp}).responder]'==0,pp) = NaN;
            
            base_n3_md_amp_contra_morph{pp} = amp_contra;
            base_n3_md_amp_ipsi_morph{pp} = amp_ipsi;
            pp = pp +1;
        end
        
    end
    
    %% Data extraction: 3 BASELINE AND MD  TUNED amplitudes, population crosssectional
    pp=1;
    for ko = 3:-1:0;
        
        if ~fitamp
            %all ODIs all ROIS
            if ~isempty(base_n3_md_morph_tuned_rois_matched{pp})
                base_n3_md_amp_contra_tuned{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}).maxAmpDelAve_contra];
                base_n3_md_amp_ipsi_tuned{pp} =  [ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}).maxAmpDelAve_ipsi];
                all_base_n3_md_ODI_tuned(:,pp) = ( base_n3_md_amp_contra_tuned{pp} - base_n3_md_amp_ipsi_tuned{pp} ) ./ (base_n3_md_amp_contra_tuned{pp} + base_n3_md_amp_ipsi_tuned{pp} );
                all_base_n3_md_ODI_nan_tuned(:,pp)  = all_base_n3_md_ODI_tuned(:,pp) ;
                all_base_n3_md_ODI_nan_tuned([ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}).responder]'==0,pp) = NaN;
                for ii = 1:length(base_n3_md_morph_tuned_rois_matched{1})
                    all_base_n3_md_OSI_contra_tuned(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}(ii)).deltapeaks_averagetrace_contra');
                    all_base_n3_md_OSI_ipsi_tuned(ii,pp) =      TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}(ii)).deltapeaks_averagetrace_ipsi');
                    %             all_base_n3_md_deltaOri_ipsi_tuned(ii,pp) = TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}(ii)).deltapeaks_averagetrace_ipsi');
                    %             all_base_n3_md_OSI_ipsi_tuned(ii,pp) =      TT_OrientationSelectivityIndex(ana_data{baseline+1-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}(ii)).deltapeaks_averagetrace_ipsi');
                    
                end
                
                all_base_n3_md_circvar_contra_tuned(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).circvar_contra];
                all_base_n3_md_circvar_ipsi_tuned(:,pp) = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).circvar_ipsi];
                
                FitIpsi = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).ipsi];
                FitContra =[ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).contra];
                
                base_n3_md_morph_tuned_rois_matched_PrefOriContra(:,pp) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                base_n3_md_morph_tuned_rois_matched_PrefOriIpsi(:,pp) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                
                [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(base_n3_md_morph_tuned_rois_matched_PrefOriContra(:,pp)*2, base_n3_md_morph_tuned_rois_matched_PrefOriIpsi(:,pp)*2) ;
                base_n3_md_morph_tuned_rois_matched_PrefDeltaOri(:,pp) = MinAngDiff /2;
                base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi(:,pp) = AngDiff /2;
                %         base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset(:,pp) = median(base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi(:,pp));
                try
                    base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset(:,pp) = tuned_both_PrefDeltaOri_bidi_offset{baseline+1-ko}; % use population offset
                catch
                    base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset(:,pp) = [];
                end
                base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset_co(:,pp) = mod([base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi(:,pp) + 90] -  base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset(:,pp),180)-90;
            else
                base_n3_md_amp_contra_tuned{pp} = [];
                base_n3_md_amp_ipsi_tuned{pp} = [];
                %             base_n3_md_morph_tuned_rois_matched = [];
                all_base_n3_md_ODI_tuned =[];
                base_n3_md_morph_tuned_rois_matched_PrefDeltaOri_bidi_offset_co = [];
            end
            pp = pp +1;
        else
            %all ODIs all ROIS
            FitIpsi = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).ipsi]; %MAKE A FUNCTION OUT OF THIS !!!!!!
            FitContra = [ana_data{baseline+1-ko}.Fit(base_n3_md_morph_tuned_rois_matched{pp}).contra];
            
            amp_contra = ([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]);
            amp_ipsi =   ([FitIpsi.PrefRsp]   + [FitIpsi.BaselineRsp])   + ([FitIpsi.OppResp]   + [FitIpsi.BaselineRsp]);
            
            all_base_n3_md_ODI(:,pp) = ( amp_contra - amp_ipsi) ./ (amp_contra + amp_ipsi);
            all_base_n3_md_ODI_nan(:,pp)  = all_n_n2_ODI_morph(:,pp) ;
            all_base_n3_md_ODI_nan([ana_data{baseline-ko}.peaks(base_n3_md_morph_tuned_rois_matched{pp}).responder]'==0,pp) = NaN;
            
            base_n3_md_amp_contra_tuned{pp} = amp_contra;
            base_n3_md_amp_ipsi_tuned{pp} = amp_ipsi;
            pp = pp +1;
        end
    end
end
% end

% try
if md_incl %selectors only necessary for MD experiments!
    
    
    %% Data extraction: binary classification of OD types over 1st MD III - Shift classes and types
    shiftmag = 0.25; %CONVERT TO BASELINE OD - BASED MAgnitude -> shift out of class!
    
    strongshifters_contra = length(find(diff(contra_base_md_ODI')<-shiftmag)) / total_cells;
    strongshifters_ipsi_MD = length(find(diff(ipsi_base_md_ODI')<-shiftmag))/ total_cells;
    
    
    noshifters_contra = (contra_md - length(find(diff(contra_base_md_ODI')<-shiftmag))) / total_cells;
    noshifters_ipsi_MD = (ipsi_MD_base - length(find(diff(ipsi_base_md_ODI')<-shiftmag)))/ total_cells;
    
    contra_lost_perc = contra_lost / total_cells;
    ipsi_gained_perc_MD = ipsi_MD_gained / total_cells;
    
    a = [contra_base_frac strongshifters_contra  contra_MD_md_frac strongshifters_ipsi_MD]';
    b = [bino_base_frac noshifters_contra  bino_MD_md_frac noshifters_ipsi_MD]';
    c = [ipsi_base_frac contra_lost_perc  ipsi_MD_md_frac ipsi_gained_perc_MD]';
    
    %% Figure 2: Bar graph of shift classes and types
    figure(24626236)
    bar_h=bar(1:4, [a b c], 0.5, 'stack', 'FaceColor' ,'w');
    bar_child=get(bar_h,'Children')
    set(bar_child{1}, 'FaceColor' ,'b');
    set(bar_child{2}, 'FaceColor' ,'w');
    set(bar_child{3}, 'FaceColor' ,'r');
    ylabel('responsive (&re fraction')
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'pre-md';  'contra pre-md'; 'post-md' ;'ipsi post-md'})
    title(mouse, 'Interpreter', 'none')
    
    %% Figure 3: Line graph ipsi/contro survivors/recruits
    figure(2352616)
    subplot(1,2,1), plot(ipsi_base_md_ODI', 'k')
    ylabel('ODI')
    set(gca, 'XTick', [1 2 ], 'XTickLabel', {'pre-MD'; 'post-MD'})
    title(mouse, 'Interpreter', 'none')
    subplot(1,2,2), plot(contra_base_md_ODI', 'k')
    ylabel('ODI')
    set(gca, 'XTick', [1 2 ], 'XTickLabel', {'pre-MD'; 'post-MD'})
    title(mouse, 'Interpreter', 'none')
    ylim([-1 1])
    xlim([0.75 2.25])
    
    %% Figure 4.1: plot amplitude changes N - N+1(aka hebbian / non-hebbian cross)
    
    figure(1992188);
    
    %generate core variables
    delta_amp_contra_n_n1_morph = [n_n1_amp_contra_morph{2} - n_n1_amp_contra_morph{1} ]';
    delta_amp_ipsi_n_n1_morph = [n_n1_amp_ipsi_morph{2} - n_n1_amp_ipsi_morph{1} ]';
    
    
    thresh = [std(delta_amp_contra_n_n1_morph) std(delta_amp_ipsi_n_n1_morph)];
    
    %z-score!
    zscore = 0;
    zthresh = 1;
    if zscore
        delta_amp_contra_n_n1_morph = (delta_amp_contra_n_n1_morph - mean(delta_amp_contra_n_n1_morph)) ./ thresh(1);
        delta_amp_ipsi_n_n1_morph = (delta_amp_ipsi_n_n1_morph - mean(delta_amp_ipsi_n_n1_morph)) ./ thresh(2);
        crit = abs(delta_amp_contra_n_n1_morph)>zthresh | abs(delta_amp_ipsi_n_n1_morph)>zthresh
    else
        crit = abs(delta_amp_contra_n_n1_morph)>thresh(1) | abs(delta_amp_ipsi_n_n1_morph)>thresh(2)
    end
    
    [theta, rho] = cart2pol(delta_amp_contra_n_n1_morph(crit),delta_amp_ipsi_n_n1_morph(crit));
    
    colorlevels = 11;
    coc = cbrewer('div', 'RdBu', colorlevels); %toned down bipolar maps... use for print
    basehdl1 = subplot(2,2,1), scatter(delta_amp_contra_n_n1_morph(crit), delta_amp_ipsi_n_n1_morph(crit), [], all_n_n1_ODI_morph(crit,1),  'fill', 'MarkerEdgeColor','k'); colormap(coc);hold on;
    subplot(2,2,1), scatter(delta_amp_contra_n_n1_morph(~crit), delta_amp_ipsi_n_n1_morph(~crit), 5, all_n_n1_ODI_morph(~crit,1)); colormap(coc)
    
    xl = xlim; yl = ylim; mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    xlim(xl);ylim(yl);hline(0); vline(0);
    title('plot hebbian / non-hebbian cross')
    
    basehdl2 =subplot(2,2,2), scatter(delta_amp_contra_n_n1_morph(crit), delta_amp_ipsi_n_n1_morph(crit), [], all_n_n1_ODI_morph(crit,2),  'fill', 'MarkerEdgeColor','k'); colormap(coc); hold on;
    subplot(2,2,2), scatter(delta_amp_contra_n_n1_morph(~crit), delta_amp_ipsi_n_n1_morph(~crit), 5, all_n_n1_ODI_morph(~crit,2)); colormap(coc)
    
    xl = xlim; yl = ylim; mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    xlim(xl);ylim(yl);hline(0); vline(0);
    title('plot hebbian / non-hebbian cross')
    
    subplot(2,2,3), rose2(theta,8);line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    subplot(2,2,4), polar(theta,rho, 'ok');line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    
    % scatter(delta_amp_contra_morph, delta_amp_ipsi_morph,abs(all_base_md_ODI_morph(:,2)-all_base_md_ODI_morph(:,1))*50,all_base_md_ODI_morph(:,1),  'fill', 'MarkerEdgeColor','k')
    % scatter(delta_amp_contra_morph, delta_amp_ipsi_morph,abs(all_base_md_ODI_morph(:,2)-all_base_md_ODI_morph(:,1))*50,all_base_md_ODI_morph(:,2),  'fill', 'MarkerEdgeColor','k')
    %
    % plot(delta_amp_contra, delta_amp_ipsi, 'ok', 'MarkerFaceColor','k'); hold on;
    % plot(mean(delta_amp_contra), mean(delta_amp_ipsi), 'or', 'MarkerFaceColor','r', 'MarkerSize', 10)
    %
    % plot(delta_amp_contra_morph, delta_amp_ipsi_morph, '*k', 'MarkerFaceColor','b'); hold on;
    % plot(mean(delta_amp_contra_morph), mean(delta_amp_ipsi_morph), '*r', 'MarkerFaceColor','b', 'MarkerSize', 10)
    
    xl = xlim;
    yl = ylim;
    mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    
    % line([-1*mxlim mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    
    xlim(xl);
    ylim(yl);
    
    hline(0); vline(0);
    title('plot hebbian / non-hebbian cross')
    %% Figure 4.2: plot amplitude changes base - MD (aka hebbian / non-hebbian cross)
    figure(199288);
    
    %generate core variables
    delta_amp_contra = [base_md_amp_contra{2} - base_md_amp_contra{1} ]';
    delta_amp_ipsi = [base_md_amp_ipsi{2} - base_md_amp_ipsi{1} ]';
    
    delta_amp_contra_morph = [base_md_amp_contra_morph{2} - base_md_amp_contra_morph{1} ]';
    delta_amp_ipsi_morph = [base_md_amp_ipsi_morph{2} - base_md_amp_ipsi_morph{1} ]';
    
    %z-score!
    if zscore
        delta_amp_contra_morph = (delta_amp_contra_morph - mean(delta_amp_contra_n_n1_morph)) ./ thresh(1);
        delta_amp_ipsi_morph = (delta_amp_ipsi_morph - mean(delta_amp_ipsi_n_n1_morph)) ./ thresh(2);
        ampcrit = abs(delta_amp_contra_morph)>zthresh | abs(delta_amp_ipsi_morph)>zthresh;
        ampcrit_resp = abs(delta_amp_contra)>zthresh | abs(delta_amp_ipsi)>zthresh;
    else
        ampcrit = abs(delta_amp_contra_morph)>thresh(1) | abs(delta_amp_ipsi_morph)>thresh(2);
        ampcrit_resp = abs(delta_amp_contra)>thresh(1) | abs(delta_amp_ipsi)>thresh(2);
    end
    % threshold is taken from n n+1 changes!
    [amptheta, amprho] = cart2pol(delta_amp_contra_morph(ampcrit),delta_amp_ipsi_morph(ampcrit));
    OD_premd_crit = all_base_md_ODI_morph(ampcrit,1);
    OD_post_crit = all_base_md_ODI_morph(ampcrit,2);
    % [~, rhoOD] = cart2pol(delta_amp_contra_morph(crit),all_base_md_ODI_morph(crit,2));
    
    colorlevels = 11;
    coc = cbrewer('div', 'RdBu', colorlevels); %toned down bipolar maps... use for print
    mdhdl1 =  subplot(2,2,1), scatter(delta_amp_contra_morph(ampcrit), delta_amp_ipsi_morph(ampcrit), [], all_base_md_ODI_morph(ampcrit,1),  'fill', 'MarkerEdgeColor','k'); colormap(coc);hold on;
    subplot(2,2,1), scatter(delta_amp_contra_morph(~ampcrit), delta_amp_ipsi_morph(~ampcrit), 5, all_base_md_ODI_morph(~ampcrit,1)); colormap(coc)
    xl = xlim; yl = ylim; mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    xlim(xl);ylim(yl);hline(0); vline(0);
    title('\Delta[Ipsi, Contra]  pre- post MD')
    
    mdhdl2 =  subplot(2,2,2), scatter(delta_amp_contra_morph(ampcrit), delta_amp_ipsi_morph(ampcrit), [], all_base_md_ODI_morph(ampcrit,2),  'fill', 'MarkerEdgeColor','k'); colormap(coc); hold on;
    subplot(2,2,2), scatter(delta_amp_contra_morph(~ampcrit), delta_amp_ipsi_morph(~ampcrit), 5, all_base_md_ODI_morph(~ampcrit,2)); colormap(coc)
    xl = xlim; yl = ylim; mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    xlim(xl);ylim(yl);hline(0); vline(0);
    title('\Delta[Ipsi, Contra] Vector Response Magnitude Distribution')
    
    xl = xlim;
    yl = ylim;
    mxlim = max([abs(max([xl yl])) abs(min([xl yl]))]);
    
    % line([-1*mxlim mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    
    xlim(xl);
    ylim(yl);
    
    hline(0); vline(0);
    title('plot hebbian / non-hebbian cross')
    
    
    subplot(2,2,3), rose2(amptheta,8);line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    % subplot(2,2,4), polar(theta,rho, 'ok');line([mxlim -1*mxlim],[mxlim -1*mxlim],'LineStyle', '--', 'Color', 'k')
    
    
    % change sectors
    ipsi_strngth_q =[135 45];
    contra_strngth_q =[45 -45];
    ipsi_depr_q = [-45 -135];
    contra_depr_q = [-135 135];
    
    ipsi_strngth_idx = find(amptheta<ipsi_strngth_q(1)*pi/180 & amptheta>ipsi_strngth_q(2)*pi/180);
    contra_strngth_idx = find(amptheta<contra_strngth_q(1)*pi/180 & amptheta>contra_strngth_q(2)*pi/180);
    ipsi_depr_idx = find(amptheta<ipsi_depr_q(1)*pi/180 & amptheta>ipsi_depr_q(2)*pi/180);
    contra_depr_idx= find(amptheta<contra_depr_q(1)*pi/180 & amptheta>contra_depr_q(2)*pi/180);
    contra_depr_idx= [find(amptheta>contra_depr_q(2)*pi/180); find(amptheta<contra_depr_q(1)*pi/180)]
    
    subplot(2,2,4);
    hold on;
    polar(amptheta(contra_strngth_idx), amprho(contra_strngth_idx), 'ob');
    polar(amptheta(ipsi_strngth_idx), amprho(ipsi_strngth_idx), 'or');
    polar(amptheta(ipsi_depr_idx), amprho(ipsi_depr_idx), '+r');
    polar(amptheta(contra_depr_idx), amprho(contra_depr_idx), '*b');hold on;
    hold off;
    
    labelmatrix = {'ipsi_strngth_q'; 'contra_strngth_q'; 'ipsi_depr_q'; 'contra_depr_q'}
    
    figure(33312)
    subplot(2,1,1), bar([1:4], [mean(OD_premd_crit(ipsi_strngth_idx)) mean(OD_premd_crit(contra_strngth_idx)) mean(OD_premd_crit(ipsi_depr_idx)) mean(OD_premd_crit(contra_depr_idx))])
    ylabel('Pre-MD ODI')
    set(gca,'XTickLabel',labelmatrix)
    subplot(2,1,2), bar([1:4], [mean(OD_post_crit(ipsi_strngth_idx)) mean(OD_post_crit(contra_strngth_idx)) mean(OD_post_crit(ipsi_depr_idx)) mean(OD_post_crit(contra_depr_idx))])
    ylabel('Post-MD ODI')
    set(gca,'XTickLabel',labelmatrix)
    
    [figure_handle,count,speeds,directions,Table] = WindRose((180.*amptheta)./pi,amprho, 'centeredin0', 1, 'ndirections',12, 'ODplot', 1,...
        'labels',{'Ipsi Strengthening','Ipsi Depression','Contra Strengthening','Contra Depression'}, 'legendvariable', '||\DeltaR/R_{0}||', ...
        'lablegend', 'Magnitude of Response Vector Change','freqlabelangle', [45], 'nspeeds', [], 'titlestring', ...
        {['ODI Rose ' mouse{1}]; '\Delta[Ipsi, Contra] Vector Response Magnitude Distribution';' '}, 'legendtype', 2 );
    
    % TODO! Do the same for BASELINE - BEST WAY OF ASSESSING VARIABILITY!
    
    
    %% Figure 5:  - population ODI over ALL timepoints
    % subtract timestamp of pre-MD recording
    TB=TB-TB(baseline);
    
    % % dump some of the core variables to workspace
    % assignin('base', 'pop_amp_ODI',pop_amp_ODI);
    % assignin('base', 'pop_amp_ODI_med',pop_amp_ODI_med);
    % assignin('base', 'TB',TB);
    
    if ~noplot
        figure(4234)
        plot(TB, pop_amp_ODI, '-ok'); hold on;
        plot(TB, pop_amp_ODI_allresp, ':ok'); hold on;
        % plot(TB, pop_amp_ODI_med, '-or'); hold on;
        % plot(TB, pop_amp_ODI_med_allresp, ':or'); hold on;
        plot(TB, pop_amp_ODI_all_full, '-or'); hold on;
        
        legend('full group responders (session wise - not matched for responsiveness throughout) [fullidx_z]', 'session wise crossectional responders [respidx_z]',  'all non-rejects [fullidx]');
        xlabel('time relative to MD onset')
        ylabel('ODI')
        title('population ODI - all ungrouped z-responders. (>8 Z-SCORE IN 50% OF THE TRIALS)')
        %
        %     if reanalyze
        %         close(h);
        %     end
        
    end
    
    %% - - - - DATA EXTRACTION: grouping of morph. and function - - - -
    
    % choose responsiveness criterion (uncomment if needed)
    
    selector = fullidx_z_rois;
    % selector = respidx_z;
    
    
    % disp('Group selection criterion: Significant response in atleast one of ALL timepoints')
    disp('Group selection criterion: Significant response (z>8)in ALL timepoints: fullidx_z_rois')
    
    combined = [];
    combined_all = [];
    
    for i = 1:length(roi_data)
        % combine all selected (e.g. signigficantly responsive) cells of each
        % session in one array
        combined = [combined roi_data{i}.ROIs(selector{i}).group];
        combined_all = [combined_all roi_data{i}.ROIs(fullidx_rois{i}).group];
    end
    
    %% Loop through groups and intersect the overlaps
    % Just take the cells that are significantly responsive either troughout or to specific timepoints.
    
    % first step: do this for the first two timepoints
    cum_int_groups = intersect([roi_data{1}.ROIs(selector{1}).group], [roi_data{2}.ROIs(selector{2}).group]);
    cum_int_base_groups = intersect([roi_data{1}.ROIs(selector{1}).group], [roi_data{2}.ROIs(selector{2}).group]);
    if baseline >3; cum_int_base_groups = intersect([roi_data{2}.ROIs(selector{2}).group], [roi_data{3}.ROIs(selector{3}).group]); end
    cum_int_all_groups = intersect([roi_data{1}.ROIs(fullidx_rois{1}).group], [roi_data{2}.ROIs(fullidx_rois{2}).group]);
    
    
    
    
    % Baseline n and n14
    cum_int_n_groups = intersect([roi_data{baseline_pair(1)}.ROIs(base_n_z_rois{baseline_pair(1)}).group], [roi_data{baseline_pair(2)}.ROIs(base_n_z_rois{baseline_pair(2)}).group]);
    if ~isempty(baseline_pair14)
        cum_int_n14_groups = intersect([roi_data{baseline_pair14(1)}.ROIs(base_n14_z_rois{baseline_pair14(1)}).group], [roi_data{baseline_pair14(2)}.ROIs(base_n14_z_rois{baseline_pair14(2)}).group]);
    end
    % just MD (Baseline and md) group
    cum_int_MD_groups = intersect([roi_data{baseline}.ROIs(selector{baseline}).group], [roi_data{baseline+1}.ROIs(selector{baseline+1}).group]);
    
    % Dual MD / Recovery (responsive at last TP efore the Ms. Good for
    % amplitude- based comparisons] NOT IMPLEMENTED
    try
        cum_int_MD2_groups = intersect([roi_data{baseline}.ROIs(selector{baseline}).group], [roi_data{baseline+recovery1}.ROIs(selector{baseline+recovery1}).group]);
    end
    
    % then loop through the rest of the data and only keep overlap GROUPS
    k = 1;
    for i = 2:length(roi_data)-1
        cum_int_groups = intersect(cum_int_groups, [roi_data{i+1}.ROIs(selector{i+1}).group]);
        cum_int_all_groups = intersect(cum_int_all_groups, [roi_data{i+1}.ROIs(fullidx_rois{i+1}).group]);
        
        if  i<= baseline
            %baseline responders
            cum_int_base_groups = intersect(cum_int_base_groups, [roi_data{i+1}.ROIs(selector{i+1}).group]);
        end
    end
    
    %% Data extraction: initialize plotarrays and select master selector
    % We have four major ways of selecting cells. Either
    % 1) responisive in at least one timepoint:
    %     groups_selected = sort(unique(combined));
    % or 2) responsive (i.e. selected) in all timepoints
    groups_selected = cum_int_groups; %COMMENT / UNCOMMENT
    % or 3) responsive during baseline
    %     groups_selected = cum_int_base;
    % or 4) only during first MD.
    %     groups_selected = cum_int_MD;
    
    %  values to be plotted - load to array - initialize as NaN to allow for
    %  dropout-means
    
    plotarray_OD = NaN(length(groups_selected),length(roi_data));
    plotarray_OD_mean = NaN(length(groups_selected),length(roi_data));
    plotarray_AmpContra = NaN(length(groups_selected),length(roi_data));
    plotarray_AmpIpsi = NaN(length(groups_selected),length(roi_data));
    plotarray_CircvarContra = NaN(length(groups_selected),length(roi_data));
    plotarray_CircvarIpsi = NaN(length(groups_selected),length(roi_data));
    plotarray_CircvarContra_ori = NaN(length(groups_selected),length(roi_data));
    plotarray_CircvarIpsi_ori = NaN(length(groups_selected),length(roi_data));
    plotarray_peak_zscore = NaN(length(groups_selected),length(roi_data));
    plotarray_mean_zscore = NaN(length(groups_selected),length(roi_data));
    plotarray_AnovIpsi = NaN(length(groups_selected),length(roi_data));
    plotarray_AnovContra = NaN(length(groups_selected),length(roi_data));
    plotarray_AmpContra_all = NaN(length(cum_int_all_groups),length(roi_data));
    plotarray_AmpIpsi_all = NaN(length(cum_int_all_groups),length(roi_data));
    
    plotarray_Contra_resp= NaN(length(groups_selected),length(roi_data));
    plotarray_Ipsi_resp= NaN(length(groups_selected),length(roi_data));
    
    plotarray_Contra_low_resp= NaN(length(groups_selected),length(roi_data));
    plotarray_Ipsi_low_resp= NaN(length(groups_selected),length(roi_data));
    
    plotarray_CircvarContra_all = NaN(length(cum_int_all_groups),length(roi_data));
    plotarray_CircvarIpsi_all = NaN(length(cum_int_all_groups),length(roi_data));
    
    plotarray_OD_pop_all =  NaN(length(cum_int_all_groups),length(roi_data));
    plotarray_OD_amp_contra_all =  NaN(length(cum_int_all_groups),length(roi_data));
    plotarray_OD_amp_ipsi_all = NaN(length(cum_int_all_groups),length(roi_data));
    
    %% Data extraction: get pixelmaps
    pixelborder = 25; %size of ROI images (include ROI at some point)
    
    if pixelmaps;
        [ROIimgs_red , ROIimgs_green ,ROIimgs_act, ROIimgs_odi, ROIimgs_ovl, ROIimgs_ori_ctr, ROIimgs_ori_ips, ROIimgs_ori_bino  majorline_r minorline_r majorline_g minorline_g line_array_g line_array_r] = cell_group_display([],mouse,exp,roi_data,pixelborder, datapath);
    end
    
    %% [+] - - - - DATA EXTRACTION: consolidate data into arrays - loop through sessions again and match groups with overlap matched cum-groups
    
    % THE MASTERLOOP STARTS HERE!
    kl =1; kl2 = 1; kl3 = 1; kl4 = 1;
    for i = 1:length(roi_data); %recordings
        
        %% Data extraction: consolitate the overlap ROI arrays
        %index into ROINr that includes all ROIs that are imaged in all
        %sessions (a) responsive (b) no responsiveness criterion
        %a)
        [~, ovlpidx{i}] = intersect([roi_data{i}.ROIs(:).group],groups_selected);
        resp_ovlpidx_rois{i} = ovlpidx{i};
        %b)
        [~, ovlpidx_all{i}] = intersect([roi_data{i}.ROIs(:).group],cum_int_all_groups);
        resp_ovlpidx_all_rois{i} = ovlpidx_all{i};
        
        [~, ovlpidx_base_all{i}] = intersect([roi_data{i}.ROIs(:).group],cum_int_base_groups);
        resp_ovlpidx_base_all_rois{i} = ovlpidx_base_all{i};
        
        
        [~, ovlpidx_base_n{i}] = intersect([roi_data{i}.ROIs(:).group],cum_int_n_groups);
        resp_ovlpidx_base_n_rois{i} = ovlpidx_base_n{i};
        
        if ~isempty(baseline_pair14)
            [~, ovlpidx_base_n14{i}] = intersect([roi_data{i}.ROIs(:).group],cum_int_n14_groups);
            resp_ovlpidx_base_n14_rois{i} = ovlpidx_base_n14{i};
        end
        %      [~, resp_ovlpidx{i}] = intersect(ovlpidx{i}', selector{i})
        
        %% Data extraction: consolidate pixelmap cell chunks in array
        if pixelmaps;
            if responly
                for kk = 1:length(groups_selected);
                    cell_array_ovl{kk,i} = ROIimgs_ovl{i}.ROI{resp_ovlpidx_rois{i}(kk)};
                    cell_array_odi{kk,i} = ROIimgs_odi{i}.ROI{resp_ovlpidx_rois{i}(kk)};
                    cell_array_ori_ctr{kk,i} = ROIimgs_ori_ctr{i}.ROI{resp_ovlpidx_rois{i}(kk)};
                    cell_array_ori_ips{kk,i} = ROIimgs_ori_ips{i}.ROI{resp_ovlpidx_rois{i}(kk)};
                    
                    if ~isempty(ROIimgs_ori_bino)
                        cell_array_ori_bino{kk,i} = ROIimgs_ori_bino{i}.ROI{resp_ovlpidx_rois{i}(kk)};
                    end
                end
            elseif baselinerois
                for kk = 1:length(cum_int_base_groups);
                    cell_array_ovl{kk,i} = ROIimgs_ovl{i}.ROI{resp_ovlpidx_base_all_rois{i}(kk)};
                    cell_array_odi{kk,i} = ROIimgs_odi{i}.ROI{resp_ovlpidx_base_all_rois{i}(kk)};
                    cell_array_ori_ctr{kk,i} = ROIimgs_ori_ctr{i}.ROI{resp_ovlpidx_base_all_rois{i}(kk)};
                    cell_array_ori_ips{kk,i} = ROIimgs_ori_ips{i}.ROI{resp_ovlpidx_base_all_rois{i}(kk)};
                    if ~isempty(ROIimgs_ori_bino)
                        cell_array_ori_bino{kk,i} = ROIimgs_ori_bino{i}.ROI{resp_ovlpidx_base_all_rois{i}(kk)};
                    end
                end
            else
                for kk = 1:length(cum_int_all_groups);
                    cell_array_ovl{kk,i} = ROIimgs_ovl{i}.ROI{resp_ovlpidx_all_rois{i}(kk)};
                    cell_array_odi{kk,i} = ROIimgs_odi{i}.ROI{resp_ovlpidx_all_rois{i}(kk)};
                    cell_array_ori_ctr{kk,i} = ROIimgs_ori_ctr{i}.ROI{resp_ovlpidx_all_rois{i}(kk)};
                    cell_array_ori_ips{kk,i} = ROIimgs_ori_ips{i}.ROI{resp_ovlpidx_all_rois{i}(kk)};
                    if ~isempty(ROIimgs_ori_bino)
                        cell_array_ori_bino{kk,i} = ROIimgs_ori_bino{i}.ROI{resp_ovlpidx_all_rois{i}(kk)};
                    end
                end
            end
        end
        
        
        
        %% Data extraction: Data arrays - z-responders all timepoints refound, all timepoints responsive
        
        
        
        %             resp_ovlpidx_responder = [ana_data{i}.peaks(resp_ovlpidx{i}).responder];
        if fitamp
            FitIpsi = [ana_data{i}.Fit(resp_ovlpidx_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(resp_ovlpidx_rois{i}).contra];
            plotarray_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
            plotarray_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
        else
            plotarray_AmpContra(:,i) = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).maxAmpDelAve_contra]';
            plotarray_AmpIpsi(:,i) = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).maxAmpDelAve_ipsi]';
            plotarray_AmpContra_trials{:,i} = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).deltatrace_trials_oris_contra];
            plotarray_AmpIpsi_trials{:,i} = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).deltatrace_trials_oris_ipsi];
        end
        for rk =1:length(resp_ovlpidx_rois{i})
            [plotarray_CircvarContra(rk,i), ~, plotarray_CircvarContra_ori(rk,i)] = TT_CircularVariance_ORI(ana_data{i}.peaks(resp_ovlpidx_rois{i}(rk)).deltapeaks_averagetrace_contra');
            [plotarray_CircvarIpsi(rk,i), ~, plotarray_CircvarIpsi_ori(rk,i)] = TT_CircularVariance_ORI(ana_data{i}.peaks(resp_ovlpidx_rois{i}(rk)).deltapeaks_averagetrace_ipsi');
        end
        % Get rid of crazy amplitudes
        if delcrazies
            crazie_idx{i} = find(plotarray_AmpContra(:,i)>cra | plotarray_AmpIpsi(:,i)>cra)
            %                 resp_ovlpidx{i}(crazie_idx{i}) = [];
        end
        plotarray_OD(:,i) =  (plotarray_AmpContra(:,i) - plotarray_AmpIpsi(:,i)) ./ (plotarray_AmpContra(:,i) + plotarray_AmpIpsi(:,i));
        
        plotarray_peak_zscore(:,i) = max(max([ana_data{i}.peaks(resp_ovlpidx_rois{i}).zscore_peaks_contra ana_data{i}.peaks(resp_ovlpidx_rois{i}).zscore_peaks_ipsi]));
        plotarray_mean_zscore(:,i) = mean(mean([ana_data{i}.peaks(resp_ovlpidx_rois{i}).zscore_peaks_contra ana_data{i}.peaks(resp_ovlpidx_rois{i}).zscore_peaks_ipsi]));
        
        plotarray_AnovIpsi(:,i) =    [ana_data{i}.peaks(resp_ovlpidx_rois{i}).Tune_Anova_maxAmpDelAve_ipsi];
        plotarray_AnovContra(:,i) =  [ana_data{i}.peaks(resp_ovlpidx_rois{i}).Tune_Anova_maxAmpDelAve_contra];
        
        plotarray_Contra_resp(:,i)= [ana_data{i}.peaks(resp_ovlpidx_rois{i}).contra_responder];
        plotarray_Ipsi_resp(:,i)  = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).ipsi_responder];
        
        plotarray_Contra_low_resp(:,i)= [ana_data{i}.peaks(resp_ovlpidx_rois{i}).contra_low_responder];
        plotarray_Ipsi_low_resp(:,i)  = [ana_data{i}.peaks(resp_ovlpidx_rois{i}).ipsi_low_responder];
        %
        %         master{1}.ana_data{1}.peaks(1).Tune_Anova_maxAmpDelAve_contra
        %% Data extraction:  Data arrays - non-matched population of responders, cross-sectional responders regardless of group size
        
        
        if fitamp
            FitIpsi = [ana_data{i}.Fit(respidx_z_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(respidx_z_rois{i}).contra];
            plotarray_OD_amp_contra{i} = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
            plotarray_OD_amp_ipsi{i} = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])'  + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])';
        else
            plotarray_OD_amp_contra{i} = [ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_contra];
            plotarray_OD_amp_ipsi{i} = [ana_data{i}.peaks(respidx_z_rois{i}).maxAmpDelAve_ipsi];
        end
        % Get rid of crazy amplitudes
        if delcrazies
            plotarray_OD_amp_contra{i}(find([plotarray_OD_amp_contra{:,i}]>cra)) = NaN;
            plotarray_OD_amp_ipsi{i}(find([plotarray_OD_amp_ipsi{:,i}]>cra)) = NaN;
            plotarray_OD_amp_contra{i}(find([plotarray_OD_amp_contra{:,i}]<0)) = 0;
            plotarray_OD_amp_ipsi{i}(find([plotarray_OD_amp_ipsi{:,i}]<0)) = 0;
        end
        plotarray_OD_pop{i} =  ([plotarray_OD_amp_contra{i}] - [plotarray_OD_amp_ipsi{i}]) ./ ([plotarray_OD_amp_contra{i}] + [plotarray_OD_amp_ipsi{i}]);
        
        %% Data extraction:  Data arrays - all full groups (no responsiveness crit) - use for test of fractional responsiveness
        
        if fitamp
            
            FitIpsi = [ana_data{i}.Fit(resp_ovlpidx_all_rois{i}).ipsi];
            FitContra = [ana_data{i}.Fit(resp_ovlpidx_all_rois{i}).contra];
            plotarray_OD_amp_contra_all(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
            plotarray_OD_amp_ipsi_all(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])';
        else
            plotarray_OD_amp_contra_all(:,i) = [ana_data{i}.peaks(resp_ovlpidx_all_rois{i}).maxAmpDelAve_contra]';
            plotarray_OD_amp_ipsi_all(:,i) = [ana_data{i}.peaks(resp_ovlpidx_all_rois{i}).maxAmpDelAve_ipsi]';
            plotarray_OD_amp_contra_all_trials{:,i} = [ana_data{i}.peaks(resp_ovlpidx_all_rois{i}).deltatrace_trials_oris_contra];
            plotarray_OD_amp_ipsi_all_trials{:,i} = [ana_data{i}.peaks(resp_ovlpidx_all_rois{i}).deltatrace_trials_oris_ipsi];
            
        end
        for rk =1:length(resp_ovlpidx_all_rois{i})
            plotarray_CircvarContra_all(rk,i) = TT_CircularVariance(ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(rk)).deltapeaks_averagetrace_contra');
            
            plotarray_CircvarIpsi_all(rk,i) = TT_CircularVariance(ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(rk)).deltapeaks_averagetrace_ipsi');
        end
        % Get rid of crazy amplitudes
        if delcrazies
            crazie_idx_all{i} = find(plotarray_OD_amp_contra_all(:,i)>cra | plotarray_OD_amp_ipsi_all(:,i)>cra)
            %                 resp_ovlpidx_all{i}(crazie_idx_all{i}) = [];
        end
        plotarray_OD_pop_all(:,i) = (plotarray_OD_amp_contra_all(:,i) - plotarray_OD_amp_ipsi_all(:,i)) ./ (plotarray_OD_amp_contra_all(:,i) + plotarray_OD_amp_ipsi_all(:,i));
        
        %% Data extraction:  Data arrays - z-responders BASELINE n n+1, tuned ipsi or contra or both
        if i == baseline_pair(1) | i == baseline_pair(2)
            
            
            % first amplitude and OD. Responsive throughout
            
            if fitamp
                FitIpsi = [ana_data{i}.Fit(resp_ovlpidx_base_n_rois{i}).ipsi];
                FitContra = [ana_data{i}.Fit(resp_ovlpidx_base_n_rois{i}).contra];
                plotarray_base_n_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
                plotarray_base_n_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
            else
                plotarray_base_n_AmpContra(:,i) =  [ana_data{i}.peaks(resp_ovlpidx_base_n_rois{i}).maxAmpDelAve_contra]';
                plotarray_base_n_AmpIpsi(:,i) = [ana_data{i}.peaks(resp_ovlpidx_base_n_rois{i}).maxAmpDelAve_ipsi]';
                
                plotarray_base_n_trials_contra{:,i} =  [ana_data{i}.peaks(resp_ovlpidx_base_n_rois{i}).deltatrace_trials_oris_contra];
                plotarray_base_n_trials_ipsi{:,i} =    [ana_data{i}.peaks(resp_ovlpidx_base_n_rois{i}).deltatrace_trials_oris_ipsi];
            end
            % Get rid of crazy amplitudes
            if delcrazies
                crazie_n_base_idx{i} = find(plotarray_base_n_AmpContra(:,i)>cra | plotarray_base_n_AmpIpsi(:,i)>cra)
                %                 resp_ovlpidx{i}(crazie_idx{i}) = [];
            end
            
            plotarray_base_n_OD(:,i) =  (plotarray_base_n_AmpContra(:,i) - plotarray_base_n_AmpIpsi(:,i)) ./ (plotarray_base_n_AmpContra(:,i) + plotarray_base_n_AmpIpsi(:,i));
            
            
            
            %Orientation/Direction tuning parameters BASELINE N N+1. TUNING
            %SELECTOR EITHER IPSI OR CONTRA
            try
                FitIpsi = [ana_data{i}.Fit(n_n1_tuned_ipsi_rois_matched{kl}).ipsi];
                FitContra = [ana_data{i}.Fit(n_n1_tuned_contra_rois_matched{kl}).contra];
                
                plotarray_base_n_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                plotarray_base_n_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                
                plotarray_base_n_PrefDirContra(:,i) = [FitContra.PrefDir];
                plotarray_base_n_PrefDirIpsi(:,i) = [FitIpsi.PrefDir];
                
                plotarray_base_n_OpPrefDirContra(:,i) = mod(round([FitContra.PrefDir]+360/2)-1, 360)+1;
                plotarray_base_n_OpPrefDirIpsi(:,i) = mod(round([FitIpsi.PrefDir]+360/2)-1, 360)+1;
                
                %                 plotarray_base_n_PrefDeltaOri(:,i) = abs(plotarray_base_n_PrefOriContra(:,i) - plotarray_base_n_PrefOriIpsi(:,i)) ;
                plotarray_base_n_TWContra(:,i) = [FitContra.Sigma]';
                plotarray_base_n_TWIspi(:,i) = [FitIpsi.Sigma]';
            end
            
            %Orientation/Direction tuning parameters TUNING MISMATCH DIFFERENT
            %SELECTOR: baseline_pair14(2) TUNED!
            try
                FitIpsi = [ana_data{i}.Fit(n_n1_tuned_both_rois_matched{kl}).ipsi];
                FitContra = [ana_data{i}.Fit(n_n1_tuned_both_rois_matched{kl}).contra];
                
                plotarray_base_n_both_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                plotarray_base_n_both_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                
                [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_n_both_PrefOriContra(:,i)*2, plotarray_base_n_both_PrefOriIpsi(:,i)*2) ;
                plotarray_base_n_PrefDeltaOri(:,i) = MinAngDiff /2;
                plotarray_base_n_PrefDeltaOri_bidi(:,i) = AngDiff /2;
                plotarray_base_n_PrefDeltaOri_bidi_offset(:,i) = median(plotarray_base_n_PrefDeltaOri_bidi(:,i));
                %             offset_array = repmat(plotarray_base_n_PrefDeltaOri_bidi_offset(:,i), length(plotarray_base_n_PrefDeltaOri_bidi(:,i)),1);
                %
                plotarray_base_n_PrefDeltaOri_bidi_offset_corr(:,i) = mod([plotarray_base_n_PrefDeltaOri_bidi(:,i) + 90] -  plotarray_base_n_PrefDeltaOri_bidi_offset(:,i),180)-90;
                
                figure(22492);hold all
                %             phdl3 = cdfplot(abs(plotarray_base_n_PrefDeltaOri_bidi(:,i)- plotarray_base_n_PrefDeltaOri_bidi_offset(:,i)));
                phdl34 = cdfplot(abs(plotarray_base_n_PrefDeltaOri_bidi(:,i)));
                phdl35 = cdfplot(abs(plotarray_base_n_PrefDeltaOri_bidi_offset_corr(:,i)));
                set(phdl35, 'LineStyle', ':')
                set([phdl34 phdl35], 'Color', colors(kl,:))
                vline(median(abs(plotarray_base_n_PrefDeltaOri_bidi(:,i))));
                title('\Delta Ori n n+1. ')
                xlabel('abs. \Delta Ori [^{\circ}]')
                % consolidate traces and tuningcurves
                for kk = 1:length(n_n1_tuned_contra_rois_matched{kl});
                    plotarray_base_n_OSIContra(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(n_n1_tuned_contra_rois_matched{kl}(kk)).contra.FittedData);
                    plotarray_base_n_DSIContra(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(n_n1_tuned_contra_rois_matched{kl}(kk)).contra.FittedData);
                end
                for kk = 1:length(n_n1_tuned_ipsi_rois_matched{kl});
                    plotarray_base_n_OSIIpsi(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(n_n1_tuned_ipsi_rois_matched{kl}(kk)).ipsi.FittedData);
                    plotarray_base_n_DSIIpsi(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(n_n1_tuned_ipsi_rois_matched{kl}(kk)).ipsi.FittedData);
                end
            end
            kl=kl+1;
        end
        %% Data extraction:  Data arrays - z-responders BASELINE n n+1 14d  tuned ipsi or contra or both
        if ~isempty(baseline_pair14)
            if i == baseline_pair14(1) | i == baseline_pair14(2)
                
                % first amplitude and OD. Responsive throughout
                
                if fitamp
                    FitIpsi = [ana_data{i}.Fit(resp_ovlpidx_base_n14_rois{i}).ipsi];
                    FitContra = [ana_data{i}.Fit(resp_ovlpidx_base_n14_rois{i}).contra];
                    plotarray_base_n14_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
                    plotarray_base_n14_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
                    plotarray_base_n14_trials_contra{:,i} =  [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).deltatrace_trials_oris_contra];
                    plotarray_base_n14_trials_ipsi{:,i} =    [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).deltatrace_trials_oris_ipsi];
                else
                    plotarray_base_n14_AmpContra(:,i) =  [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).maxAmpDelAve_contra]';
                    plotarray_base_n14_AmpIpsi(:,i) = [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).maxAmpDelAve_ipsi]';
                    
                    plotarray_base_n14_trials_contra{:,i} =  [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).deltatrace_trials_oris_contra];
                    plotarray_base_n14_trials_ipsi{:,i} =    [ana_data{i}.peaks(resp_ovlpidx_base_n14_rois{i}).deltatrace_trials_oris_ipsi];
                    
                    
                end
                % Get rid of crazy amplitudes
                if delcrazies
                    crazie_n14_base_idx{i} = find(plotarray_base_n14_AmpContra(:,i)>cra | plotarray_base_n14_AmpIpsi(:,i)>cra)
                    %                 resp_ovlpidx{i}(crazie_idx{i}) = [];
                end
                
                plotarray_base_n14_OD(:,i) =  (plotarray_base_n14_AmpContra(:,i) - plotarray_base_n14_AmpIpsi(:,i)) ./ (plotarray_base_n14_AmpContra(:,i) + plotarray_base_n14_AmpIpsi(:,i));
                
                
                
                %Orientation/Direction tuning parameters BASELINE N N+1. TUNING
                %SELECTOR EITHER IPSI OR CONTRA
                try
                    FitIpsi = [ana_data{i}.Fit(n_n1_14d_tuned_ipsi_rois_matched{kl2}).ipsi];
                    FitContra = [ana_data{i}.Fit(n_n1_14d_tuned_contra_rois_matched{kl2}).contra];
                    
                    plotarray_base_n14_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                    plotarray_base_n14_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                    
                    plotarray_base_n14_PrefDirContra(:,i) = [FitContra.PrefDir];
                    plotarray_base_n14_PrefDirIpsi(:,i) = [FitIpsi.PrefDir];
                    
                    plotarray_base_n14_OpPrefDirContra(:,i) = mod(round([FitContra.PrefDir]+360/2)-1, 360)+1;
                    plotarray_base_n14_OpPrefDirIpsi(:,i) = mod(round([FitIpsi.PrefDir]+360/2)-1, 360)+1;
                    
                    %                 plotarray_base_n14_PrefDeltaOri(:,i) = abs(plotarray_base_n14_PrefOriContra(:,i) - plotarray_base_n14_PrefOriIpsi(:,i)) ;
                    plotarray_base_n14_TWContra(:,i) = [FitContra.Sigma]';
                    plotarray_base_n14_TWIspi(:,i) = [FitIpsi.Sigma]';
                end
                
                %Orientation/Direction tuning parameters TUNING MISMATCH DIFFERENT
                %SELECTOR: baseline_pair14(2) TUNED!
                try
                    FitIpsi = [ana_data{i}.Fit(n_n1_14d_tuned_both_rois_matched{kl2}).ipsi];
                    FitContra = [ana_data{i}.Fit(n_n1_14d_tuned_both_rois_matched{kl2}).contra];
                    
                    plotarray_base_n14_both_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                    plotarray_base_n14_both_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                    
                    [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_n14_both_PrefOriContra(:,i)*2, plotarray_base_n14_both_PrefOriIpsi(:,i)*2) ;
                    plotarray_base_n14_PrefDeltaOri(:,i) = MinAngDiff /2;
                    plotarray_base_n14_PrefDeltaOri_bidi(:,i) = AngDiff /2;
                    plotarray_base_n14_PrefDeltaOri_bidi_offset(:,i) = median(plotarray_base_n14_PrefDeltaOri_bidi(:,i));
                    
                    %                 offset_array = repmat(plotarray_base_n_PrefDeltaOri_bidi_offset(:,i), length(plotarray_base_n_PrefDeltaOri_bidi(:,i)),1);
                    plotarray_base_n14_PrefDeltaOri_bidi_offset_corr(:,i) = mod([plotarray_base_n14_PrefDeltaOri_bidi(:,i) + 90] - plotarray_base_n14_PrefDeltaOri_bidi_offset(:,i),180)-90;
                    
                    figure(224923);hold all
                    %                 plhdl1 = cdfplot(abs(plotarray_base_n14_PrefDeltaOri_bidi(:,i)- plotarray_base_n14_PrefDeltaOri_bidi_offset(:,i)));
                    plhdl1 = cdfplot(abs(plotarray_base_n14_PrefDeltaOri_bidi(:,i)));
                    plhdl11 = cdfplot(abs(plotarray_base_n14_PrefDeltaOri_bidi_offset_corr(:,i)));
                    set(plhdl11, 'LineStyle', ':')
                    set([plhdl1 plhdl11], 'Color', colors(kl2,:))
                    
                    vline(median(abs(plotarray_base_n14_PrefDeltaOri_bidi(:,i))));
                    title('\Delta Ori n n+1 14d. Corrected for mean DeltaOri')
                    xlabel('abs. \Delta Ori [^{\circ}]')
                    % consolidate traces and tuningcurves
                    for kk = 1:length(n_n1_14d_tuned_contra_rois_matched{kl2});
                        plotarray_base_n14_OSIContra(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(n_n1_14d_tuned_contra_rois_matched{kl2}(kk)).contra.FittedData);
                        plotarray_base_n14_DSIContra(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(n_n1_14d_tuned_contra_rois_matched{kl2}(kk)).contra.FittedData);
                    end
                    for kk = 1:length(n_n1_14d_tuned_ipsi_rois_matched{kl2});
                        plotarray_base_n14_OSIIpsi(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(n_n1_14d_tuned_ipsi_rois_matched{kl2}(kk)).ipsi.FittedData);
                        plotarray_base_n14_DSIIpsi(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(n_n1_14d_tuned_ipsi_rois_matched{kl2}(kk)).ipsi.FittedData);
                    end
                end
                kl2=kl2+1;
            end
        end
        %% Data extraction:  Data arrays - z-responders BASELINE MD  tuned ipsi or contra or both
        if ~isempty(baseline)
            try
                if i == baseline | i == baseline +1
                    
                    % first amplitude and OD. Responsive throughout
                    
                    if fitamp
                        FitIpsi = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).contra];
                        plotarray_base_md_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
                        plotarray_base_md_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
                    else
                        plotarray_base_md_AmpContra(:,i) =  [ana_data{i}.peaks(base_md_morph_rois_tuned_matched{kl4}).maxAmpDelAve_contra]';
                        plotarray_base_md_AmpIpsi(:,i) = [ana_data{i}.peaks(base_md_morph_rois_tuned_matched{kl4}).maxAmpDelAve_ipsi]';
                    end
                    % Get rid of crazy amplitudes
                    if delcrazies
                        crazie_rec1_base_idx{i} = find(plotarray_base_md_AmpContra(:,i)>cra | plotarray_base_md_AmpIpsi(:,i)>cra)
                        %                 resp_ovlpidx{i}(crazie_idx{i}) = [];
                    end
                    
                    plotarray_base_md_tuned_OD(:,i) =  (plotarray_base_md_AmpContra(:,i) - plotarray_base_md_AmpIpsi(:,i)) ./ (plotarray_base_md_AmpContra(:,i) + plotarray_base_md_AmpIpsi(:,i));
                    
                    % first amplitude and OD. Only tuned at last baseline point
                    
                    
                    
                    if fitamp
                        FitIpsi = [ana_data{i}.Fit(base_md_morph_rois_base_tuned_matched{kl4}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_md_morph_rois_base_tuned_matched{kl4}).contra];
                        plotarray_base_md_base_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
                        plotarray_base_md_base_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
                    else
                        plotarray_base_md_base_AmpContra(:,i) =  [ana_data{i}.peaks(base_md_morph_rois_base_tuned_matched{kl4}).maxAmpDelAve_contra]';
                        plotarray_base_md_base_AmpIpsi(:,i) = [ana_data{i}.peaks(base_md_morph_rois_base_tuned_matched{kl4}).maxAmpDelAve_ipsi]';
                    end
                    
                    
                    %Orientation/Direction tuning parameters BASELINE and MD TUNING
                    
                    try
                        FitIpsi = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).contra];
                        
                        plotarray_base_md_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                        plotarray_base_md_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                        
                        plotarray_base_md_PrefDirContra(:,i) = [FitContra.PrefDir];
                        plotarray_base_md_PrefDirIpsi(:,i) = [FitIpsi.PrefDir];
                        
                        plotarray_base_md_OpPrefDirContra(:,i) = mod(round([FitContra.PrefDir]+360/2)-1, 360)+1;
                        plotarray_base_md_OpPrefDirIpsi(:,i) = mod(round([FitIpsi.PrefDir]+360/2)-1, 360)+1;
                        
                        %                 plotarray_base_md_PrefDeltaOri(:,i) = abs(plotarray_base_md_PrefOriContra(:,i) - plotarray_base_md_PrefOriIpsi(:,i)) ;
                        plotarray_base_md_TWContra(:,i) = [FitContra.Sigma]';
                        plotarray_base_md_TWIspi(:,i) = [FitIpsi.Sigma]';
                    end
                    %Orientation/Direction tuning parameters TUNING MISMATCH DIFFERENT
                    %SELECTOR: baseline only tuned
                    try
                        FitIpsi = [ana_data{i}.Fit(base_md_morph_rois_base_tuned_matched{kl4}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_md_morph_rois_base_tuned_matched{kl4}).contra];
                        
                        plotarray_base_md_base_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                        plotarray_base_md_base_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                        
                        [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_md_base_PrefOriContra(:,i)*2, plotarray_base_md_base_PrefOriIpsi(:,i)*2) ;
                        plotarray_base_md_base_PrefDeltaOri(:,i) = MinAngDiff /2;
                        plotarray_base_md_base_PrefDeltaOri_bidi(:,i) = AngDiff /2;
                        plotarray_base_md_base_PrefDeltaOri_bidi_offset(:,i) = median(plotarray_base_md_base_PrefDeltaOri_bidi(:,i));
                        plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr(:,i) = mod([plotarray_base_md_base_PrefDeltaOri_bidi(:,i) + 90] - plotarray_base_md_base_PrefDeltaOri_bidi_offset(:,i),180)-90;
                    end
                    %                 plotarray_base_md_base_PrefDeltaOri_offset_corr(:,i)  = abs(plotarray_base_md_base_PrefDeltaOri_bidi(:,i)- plotarray_base_md_base_PrefDeltaOri_bidi_offset(:,i));
                    
                    
                    
                    
                    %Orientation/Direction tuning parameters TUNING MISMATCH
                    %DIFFERENT (REDUNDANT! SEE ABOVE)
                    %SELECTOR: baseline and mD both tuned!
                    try
                        FitIpsi = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}).contra];
                        
                        plotarray_base_md_both_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                        plotarray_base_md_both_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                        
                        [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_md_both_PrefOriContra(:,i)*2, plotarray_base_md_both_PrefOriIpsi(:,i)*2) ;
                        plotarray_base_md_PrefDeltaOri(:,i) = MinAngDiff /2;
                        plotarray_base_md_PrefDeltaOri_bidi(:,i) = AngDiff /2;
                        plotarray_base_md_PrefDeltaOri_bidi_offset(:,i) = median(plotarray_base_md_PrefDeltaOri_bidi(:,i));
                        plotarray_base_md_PrefDeltaOri_bidi_offset_corr(:,i) = mod([plotarray_base_md_PrefDeltaOri_bidi(:,i) + 90] - plotarray_base_md_PrefDeltaOri_bidi_offset(:,i),180)-90;
                    end
                    
                    figure(249236);hold all
                    %                 plhdl2 =  cdfplot(abs(plotarray_base_md_PrefDeltaOri_bidi(:,i)- plotarray_base_md_PrefDeltaOri_bidi_offset(:,i)));
                    plhdl2 =  cdfplot(abs(plotarray_base_md_PrefDeltaOri_bidi(:,i)));
                    plhdl21 =  cdfplot(abs(plotarray_base_md_PrefDeltaOri_bidi_offset_corr(:,i)));
                    set(plhdl21, 'LineStyle', ':')
                    set([plhdl2 plhdl21], 'Color', colors(kl4,:)) ;
                    vline(median(abs(plotarray_base_md_PrefDeltaOri_bidi(:,i))));
                    title('\Delta Ori base MD Corrected for mean DeltaOri')
                    xlabel('abs. \Delta Ori [^{\circ}]')
                    % consolidate traces and tuningcurves
                    for kk = 1:length(base_md_morph_rois_tuned_matched{kl4});
                        plotarray_base_md_OSIContra(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}(kk)).contra.FittedData);
                        plotarray_base_md_DSIContra(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}(kk)).contra.FittedData);
                    end
                    for kk = 1:length(base_md_morph_rois_tuned_matched{kl4});
                        plotarray_base_md_OSIIpsi(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}(kk)).ipsi.FittedData);
                        plotarray_base_md_DSIIpsi(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(base_md_morph_rois_tuned_matched{kl4}(kk)).ipsi.FittedData);
                    end
                    kl4=kl4+1;
                end
            end
        end
        %% Data extraction:  Data arrays - z-responders BASELINE Recovery  tuned ipsi or contra or both
        if ~isempty(recovery1)
            if i == baseline | i == baseline +recovery1
                if ~isempty(base_rec_morph_rois_tuned_matched{kl3})
                    % first amplitude and OD. Responsive throughout
                    
                    
                    if fitamp
                        FitIpsi = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).ipsi];
                        FitContra = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).contra];
                        plotarray_base_rec1_AmpContra(:,i) = ([FitContra.PrefRsp] + [FitContra.BaselineRsp])' + ([FitContra.OppResp] + [FitContra.BaselineRsp])';
                        plotarray_base_rec1_AmpIpsi(:,i) = ([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp])' + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])'
                    else
                        plotarray_base_rec1_AmpContra(:,i) =  [ana_data{i}.peaks(base_rec_morph_rois_tuned_matched{kl3}).maxAmpDelAve_contra]';
                        plotarray_base_rec1_AmpIpsi(:,i) = [ana_data{i}.peaks(base_rec_morph_rois_tuned_matched{kl3}).maxAmpDelAve_ipsi]';
                    end
                    % Get rid of crazy amplitudes
                    if delcrazies
                        crazie_rec1_base_idx{i} = find(plotarray_base_rec1_AmpContra(:,i)>cra | plotarray_base_rec1_AmpIpsi(:,i)>cra)
                        %                 resp_ovlpidx{i}(crazie_idx{i}) = [];
                    end
                    
                    plotarray_base_rec1_OD(:,i) =  (plotarray_base_rec1_AmpContra(:,i) - plotarray_base_rec1_AmpIpsi(:,i)) ./ (plotarray_base_rec1_AmpContra(:,i) + plotarray_base_rec1_AmpIpsi(:,i));
                    
                    
                    
                    %Orientation/Direction tuning parameters BASELINE N N+1. TUNING
                    %SELECTOR EITHER IPSI OR CONTRA
                    
                    FitIpsi = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).ipsi];
                    FitContra = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).contra];
                    
                    plotarray_base_rec1_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                    plotarray_base_rec1_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                    
                    plotarray_base_rec1_PrefDirContra(:,i) = [FitContra.PrefDir];
                    plotarray_base_rec1_PrefDirIpsi(:,i) = [FitIpsi.PrefDir];
                    
                    plotarray_base_rec1_OpPrefDirContra(:,i) = mod(round([FitContra.PrefDir]+360/2)-1, 360)+1;
                    plotarray_base_rec1_OpPrefDirIpsi(:,i) = mod(round([FitIpsi.PrefDir]+360/2)-1, 360)+1;
                    
                    %                 plotarray_base_rec1_PrefDeltaOri(:,i) = abs(plotarray_base_rec1_PrefOriContra(:,i) - plotarray_base_rec1_PrefOriIpsi(:,i)) ;
                    plotarray_base_rec1_TWContra(:,i) = [FitContra.Sigma]';
                    plotarray_base_rec1_TWIspi(:,i) = [FitIpsi.Sigma]';
                    
                    
                    %Orientation/Direction tuning parameters TUNING MISMATCH DIFFERENT
                    %SELECTOR: baseline_pair14(2) TUNED!
                    
                    FitIpsi = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).ipsi];
                    FitContra = [ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}).contra];
                    
                    plotarray_base_rec1_both_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
                    plotarray_base_rec1_both_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
                    
                    [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_base_rec1_both_PrefOriContra(:,i)*2, plotarray_base_rec1_both_PrefOriIpsi(:,i)*2) ;
                    plotarray_base_rec1_PrefDeltaOri(:,i) = MinAngDiff /2;
                    plotarray_base_rec1_PrefDeltaOri_bidi(:,i) = AngDiff /2;
                    plotarray_base_rec1_PrefDeltaOri_bidi_offset(:,i) = median(plotarray_base_rec1_PrefDeltaOri_bidi(:,i));
                    plotarray_base_rec1_PrefDeltaOri_bidi_offset_corr(:,i) = mod([plotarray_base_rec1_PrefDeltaOri_bidi(:,i) + 90] - plotarray_base_rec1_PrefDeltaOri_bidi_offset(:,i),180)-90;
                    
                    
                    
                    
                    figure(2249236);hold all
                    cdfplot(abs(plotarray_base_rec1_PrefDeltaOri_bidi(:,i)- plotarray_base_rec1_PrefDeltaOri_bidi_offset(:,i)))
                    vline(median(abs(plotarray_base_rec1_PrefDeltaOri_bidi(:,i)- plotarray_base_rec1_PrefDeltaOri_bidi_offset(:,i))))
                    title('\Delta Ori base rec. Corrected for mean DeltaOri')
                    xlabel('abs. \Delta Ori [^{\circ}]')
                    % consolidate traces and tuningcurves
                    for kk = 1:length(base_rec_morph_rois_tuned_matched{kl3});
                        plotarray_base_rec1_OSIContra(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}(kk)).contra.FittedData);
                        plotarray_base_rec1_DSIContra(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}(kk)).contra.FittedData);
                    end
                    for kk = 1:length(base_rec_morph_rois_tuned_matched{kl3});
                        plotarray_base_rec1_OSIIpsi(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}(kk)).ipsi.FittedData);
                        plotarray_base_rec1_DSIIpsi(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(base_rec_morph_rois_tuned_matched{kl3}(kk)).ipsi.FittedData);
                    end
                    kl3=kl3+1;
                end
            end
        end
        %% Orientation/Direction tuning parameters ALL TIMEPOINTS Responsive and refound
        
        FitIpsi = [ana_data{i}.Fit(resp_ovlpidx_rois{i}).ipsi];
        FitContra = [ana_data{i}.Fit(resp_ovlpidx_rois{i}).contra];
        
        plotarray_PrefOriContra(:,i) = [mod([FitContra.PrefDir]+180-1, 180)+1]';
        plotarray_PrefOriIpsi(:,i) = [mod([FitIpsi.PrefDir]+180-1, 180)+1]';
        
        plotarray_PrefDirContra(:,i) = [FitContra.PrefDir];
        plotarray_PrefDirIpsi(:,i) = [FitIpsi.PrefDir];
        
        plotarray_OpPrefDirContra(:,i) = mod(round([FitContra.PrefDir]+360/2)-1, 360)+1;
        plotarray_OpPrefDirIpsi(:,i) = mod(round([FitIpsi.PrefDir]+360/2)-1, 360)+1;
        
        [MinAngDiff, AngDiff, ~, ~] = TT_AngularDifference(plotarray_PrefOriContra(:,i)*2, plotarray_PrefOriIpsi(:,i)*2);
        plotarray_PrefDeltaOri(:,i) = MinAngDiff /2;
        plotarray_PrefDeltaOri_bidi(:,i) = AngDiff / 2;
        plotarray_PrefDeltaOri_bidi_offset(:,i) = mean(plotarray_PrefDeltaOri_bidi(:,i));
        plotarray_PrefDeltaOri_bidi(:,i) =plotarray_PrefDeltaOri_bidi(:,i) - plotarray_PrefDeltaOri_bidi_offset(:,i);
        
        plotarray_TWContra(:,i) = [FitContra.Sigma]';
        plotarray_TWIspi(:,i) = [FitIpsi.Sigma]';
        
        %             plotarray_DSIContra(:,i) = (([FitContra.PrefRsp] + [FitContra.BaselineRsp]) -  ([FitContra.OppResp] + [FitContra.BaselineRsp])) ./ (([FitContra.PrefRsp] + [FitContra.BaselineRsp]) + ([FitContra.OppResp] + [FitContra.BaselineRsp]));
        %             plotarray_DSIIpsi(:,i) = (([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) -  ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp])) ./ (([FitIpsi.PrefRsp] + [FitIpsi.BaselineRsp]) + ([FitIpsi.OppResp] + [FitIpsi.BaselineRsp]));
        
        
        % consolidate traces and tuningcurves
        maxsamp = 94
        for kk = 1:length(groups_selected);
            [~,midx_c] =  max(ana_data{i}.peaks(resp_ovlpidx_rois{i}(kk)).deltapeaks_averagetrace_contra);
            [~,midx_i] =  max(ana_data{i}.peaks(resp_ovlpidx_rois{i}(kk)).deltapeaks_averagetrace_ipsi);
            
            plotarray_PrefOriContra_PSTH{kk,i} = ana_data{i}.peaks(resp_ovlpidx_rois{i}(kk)).deltatrace_averagetrace_oris_contra{midx_c}(1:maxsamp);
            plotarray_PrefOriIpsi_PSTH{kk,i} = ana_data{i}.peaks(resp_ovlpidx_rois{i}(kk)).deltatrace_averagetrace_oris_ipsi{midx_i}(1:maxsamp);
            
            plotarray_360Dir_ctr_tuning{kk,i} = ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).contra.FittedData;
            plotarray_360Dir_ips_tuning{kk,i} = ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).ipsi.FittedData;
            
            plotarray_180Ori_ctr_tuning{kk,i} = TT_FoldedTuningCurve(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).contra.FittedData);
            plotarray_180Ori_ips_tuning{kk,i} = TT_FoldedTuningCurve(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).ipsi.FittedData);
            
            plotarray_OSIContra(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).contra.FittedData);
            plotarray_OSIIpsi(kk,i) = TT_OrientationSelectivityIndex(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).ipsi.FittedData);
            
            plotarray_DSIContra(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).contra.FittedData);
            plotarray_DSIIpsi(kk,i) =  TT_DirectionIndex(ana_data{i}.Fit(resp_ovlpidx_rois{i}(kk)).ipsi.FittedData);
        end
        for kk = 1:length(resp_ovlpidx_all_rois{1});
            [~,midxa_c] = max(ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(kk)).deltapeaks_averagetrace_contra);
            [~,midxa_i] = max(ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(kk)).deltapeaks_averagetrace_ipsi);
            plotarray_PrefOriContra_PSTH_a{kk,i} = ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(kk)).deltatrace_averagetrace_oris_contra{midxa_c}(1:maxsamp);
            plotarray_PrefOriIpsi_PSTH_a{kk,i} = ana_data{i}.peaks(resp_ovlpidx_all_rois{i}(kk)).deltatrace_averagetrace_oris_ipsi{midxa_i}(1:maxsamp);
        end
        % THE MASTERLOOP ENDS HERE!
    end
    % try
    %% Figure XXX: Delta ORI vs. abs ODI and shift magnitude
    
    
    try
        figure(141820)
        deltaODI_base_md_tuned = plotarray_base_md_tuned_OD(:,baseline+1)-plotarray_base_md_tuned_OD(:,baseline);
        plot(abs(plotarray_base_md_PrefDeltaOri_bidi(:,baseline)), abs(deltaODI_base_md_tuned), '+k');lsline; hold all
        plot(abs(plotarray_base_md_PrefDeltaOri_bidi(:,baseline+1)), abs(deltaODI_base_md_tuned), '+r');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}]'); ylabel('\Delta ODI (red: post-mD, black pre-MD)');
        
        figure(141849)
        plot(abs(plotarray_base_md_PrefDeltaOri_bidi_offset_corr(:,baseline)), abs(deltaODI_base_md_tuned), '+k');lsline; hold all
        plot(abs(plotarray_base_md_PrefDeltaOri_bidi_offset_corr(:,baseline+1)), abs(deltaODI_base_md_tuned), '+r');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}] - median'); ylabel('\Delta ODI (red: post-mD, black pre-MD)');
        
        base_md_base_tuned_delta_amp_contra =  plotarray_base_md_base_AmpContra(:,baseline+1) - plotarray_base_md_base_AmpContra(:,baseline);
        base_md_base_tuned_delta_amp_ipsi =  plotarray_base_md_base_AmpIpsi(:,baseline+1) - plotarray_base_md_base_AmpIpsi(:,baseline);
        [shift_dir, shift_mag_tuned] = cart2pol(base_md_base_tuned_delta_amp_contra,base_md_base_tuned_delta_amp_ipsi);
        
        
        try
            figure(1422220)
            deltaODI_base_rec_tuned = plotarray_base_rec1_OD(:,baseline)-plotarray_base_rec1_OD(:,baseline+recovery1);
            plot(abs(plotarray_base_rec1_PrefDeltaOri_bidi(:,baseline)), (deltaODI_base_rec_tuned), '+k');lsline; hold all
            plot(abs(plotarray_base_rec1_PrefDeltaOri_bidi(:,baseline+recovery1)), (deltaODI_base_rec_tuned), '+r');lsline; hold all
            xlabel('abs. \Delta Ori [^{\circ}]'); ylabel('\Delta ODI recover vs. baseline');
        end
        
        figure(1421825)
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi(:,baseline)), base_md_base_tuned_delta_amp_contra, '+b');lsline; hold all
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi(:,baseline)), base_md_base_tuned_delta_amp_ipsi, '+r');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}]'); ylabel('\Delta R/R_0');
        
        figure(1421429)
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr(:,baseline)), base_md_base_tuned_delta_amp_contra, '+b');lsline; hold all
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr(:,baseline)), base_md_base_tuned_delta_amp_ipsi, '+r');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}] - median'); ylabel('\Delta R/R_0');
        
        figure(1221823)
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi(:,baseline)), shift_dir, '+k');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}]'); ylabel('\DeltaMD Vector Direction ||Resp.|| [\DeltaR/R_0]');
        
        
        figure(1421820)
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi(:,baseline)), shift_mag_tuned, '+k');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}]'); ylabel('\DeltaMD Vector Magnitude ||Resp.|| [\DeltaR/R_0]');
        
        figure(1421829)
        plot(abs(plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr(:,baseline)), shift_mag_tuned, '+k');lsline; hold all
        xlabel('abs. \Delta Ori [^{\circ}] - median'); ylabel('\DeltaMD Vector Magnitude ||Resp.|| [\DeltaR/R_0]');
        
        
        figure(1221821)
        subplot(2,1,1), hist(plotarray_base_md_base_PrefDeltaOri_bidi(:,baseline))
        vline(plotarray_base_md_base_PrefDeltaOri_bidi_offset(:,baseline)); xlabel('\Delta Ori (base)'); ylabel('n');
        
        subplot(2,1,2), hist(plotarray_base_md_PrefDeltaOri_bidi(:,baseline))
        vline(plotarray_base_md_PrefDeltaOri_bidi_offset(:,baseline)); xlabel('\Delta Ori(base & MD)'); ylabel('n');
        
        
        figure(1221825)
        subplot(2,1,1), hist(plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr(:,baseline))
        vline(plotarray_base_md_base_PrefDeltaOri_bidi_offset(:,baseline)); xlabel('\Delta Ori (base) - median'); ylabel('n');
        
        subplot(2,1,2), hist(plotarray_base_md_PrefDeltaOri_bidi_offset_corr(:,baseline))
        vline(plotarray_base_md_PrefDeltaOri_bidi_offset(:,baseline)); xlabel('\Delta Ori(base & MD) - median'); ylabel('n');
        
    end
    %     plotarray_base_md_base_PrefDeltaOri_bidi_offset_corr
    % end
    %     disp_groups(resp_ovlpidx_base_n, [baseline_pair(1) baseline_pair(2)], [], mouse, exp ,roi_data, 25);
    
    %% Data extraction: remove outlier groups (crazies...)
    if delcrazies
        % fullidx_z throughout (full groups, always responsive)
        killgroups = unique(cell2mat(crazie_idx'));
        plotarray_OD(killgroups,:) = [];
        plotarray_AmpContra(killgroups,:) = [];
        plotarray_AmpIpsi(killgroups,:) = [];
        plotarray_peak_zscore(killgroups,:) = [];
        plotarray_mean_zscore(killgroups,:) = [];
        
        plotarray_PrefOriContra(killgroups,:) = [];
        plotarray_PrefOriIpsi(killgroups,:) = [];
        
        plotarray_PrefDirContra(killgroups,:) = [];
        plotarray_PrefDirIpsi(killgroups,:) = [];
        
        plotarray_OpPrefDirContra(killgroups,:) = [];
        plotarray_OpPrefDirIpsi(killgroups,:) = [];
        
        plotarray_PrefDeltaOri(killgroups,:) = [];
        plotarray_PrefDeltaOri_bidi(killgroups,:) = [];
        plotarray_TWContra(killgroups,:) = [];
        plotarray_TWIspi(killgroups,:) = [];
        
        plotarray_PrefOriContra_PSTH(killgroups,:) = [];
        plotarray_PrefOriIpsi_PSTH(killgroups,:) = [];
        plotarray_360Dir_ctr_tuning(killgroups,:) = [];
        plotarray_360Dir_ips_tuning(killgroups,:) = [];
        
        plotarray_180Ori_ctr_tuning(killgroups,:) = [];
        plotarray_180Ori_ips_tuning(killgroups,:) = [];
        
        plotarray_OSIContra(killgroups,:) = [];
        plotarray_OSIIpsi(killgroups,:) = [];
        
        plotarray_DSIContra(killgroups,:) =  [];
        plotarray_DSIIpsi(killgroups,:) = [];
        
        % fullidx (full groups, no responsiveness crit)
        killgroups_all = unique(cell2mat(crazie_idx_all'));
        plotarray_OD_pop_all(killgroups_all,:) = [];
        plotarray_OD_amp_contra_all(killgroups_all,:) = [];
        plotarray_OD_amp_ipsi_all(killgroups_all,:) = [];
        plotarray_PrefOriContra_PSTH_a(killgroups_all,:) = [];
        plotarray_PrefOriIpsi_PSTH_a(killgroups_all,:) = [];
        
        
        killgroups_base = unique(cell2mat(crazie_n_base_idx'));
        plotarray_base_n_OD(killgroups_base,:) = [];
        plotarray_base_n_AmpContra(killgroups_base,:) = [];
        plotarray_base_n_AmpIpsi(killgroups_base,:) = [];
        
        
        plotarray_base_n_PrefOriContra(killgroups_base,:) = [];
        plotarray_base_n_PrefOriIpsi(killgroups_base,:) = [];
        
        plotarray_base_n_PrefDirContra(killgroups_base,:) = [];
        plotarray_base_n_PrefDirIpsi(killgroups_base,:) = [];
        
        plotarray_base_n_OpPrefDirContra(killgroups_base,:) = [];
        plotarray_base_n_OpPrefDirIpsi(killgroups_base,:) = [];
        
        plotarray_base_n_PrefDeltaOri(killgroups_base,:) = [];
        plotarray_base_n_TWContra(killgroups_base,:) = [];
        plotarray_base_n_TWIspi(killgroups_base,:) = [];
        plotarray_base_n_PrefDeltaOri_bidi(killgroups_base,:) = [];
        %         plotarray_base_n_PrefDeltaOri_bidi_offset(killgroups_base,:) =[];
        
        
        plotarray_base_n_OSIContra(killgroups_base,:) = [];
        plotarray_base_n_OSIIpsi(killgroups_base,:) = [];
        
        plotarray_base_n_DSIContra(killgroups_base,:) =  [];
        plotarray_base_n_DSIIpsi(killgroups_base,:) = [];
        
        
        for kk = 1:length(roi_data);
            resp_ovlpidx_rois{kk}(killgroups) =  [];
            resp_ovlpidx_all_rois{kk}(killgroups_all) =  [];
            try % this is a noncomplete group index.
                resp_ovlpidx_base_n_rois{kk}(killgroups_base) =  [];
            end
            
            %to allow reindexing based on resp_ovlpidx_all (the main
            %amplitude index w/o threshold)
            [~, resp_ovlpidx_all_responder{kk}] = intersect(resp_ovlpidx_all_rois{kk}, resp_ovlpidx_rois{kk});
            [~, resp_ovlpidx_base_n_responder{kk}] = intersect(resp_ovlpidx_all_rois{kk}, resp_ovlpidx_base_n_rois{kk});
            [~, resp_ovlpidx_cross_z_responder{kk}] = intersect(resp_ovlpidx_all_rois{kk}, respidx_z_rois{kk});
        end
        
    end
    
    %% Figure 6: grouped pixel map ROI images
    if pixelmaps
        if sortmap && ~responly
            sortOD = mean(plotarray_OD_pop_all(:,1:baseline),2) .* (mean(plotarray_OD_amp_contra_all(:,1:baseline),2) + mean(plotarray_OD_amp_ipsi_all(:,1:baseline),2)) ;
            [~, srtidx] =sort(sortOD);
        elseif  sortmap && responly
            sortOD = mean(plotarray_OD(:,baseline:baseline),2); %last baseline response for sorting
            %         sortOD = mean(plotarray_OD(:,1:baseline),2); %mean of last  baseline responses for sorting
            [~, srtidx] =sort(sortOD, 'descend'); %descend: contra ->ipsi
            
        else
            srtidx = 1:size(cell_array_ovl,1);
        end
        
        
        if sigrespmap
            sm = 1;
            sigshift = abs(plotarray_OD(:,baseline+1)-plotarray_OD(:,baseline))>sm * std(plotarray_OD(:,baseline-2:baseline),[],2);
            srtidx = srtidx(sigshift);
        end
        
        ovlp_img = cell2mat(cell_array_ovl(srtidx,:));
        odi_img = cell2mat(cell_array_odi(srtidx,:));
        ori_ctr_img = cell2mat(cell_array_ori_ctr(srtidx,:));
        ori_ips_img = cell2mat(cell_array_ori_ips(srtidx,:));
        
        if ~isempty(ROIimgs_ori_bino)
            ori_bino_img = cell2mat(cell_array_ori_bino(srtidx,:));
        end
        
        %         adjust_contrast = 1;
        splitmd = 0;
        blck =0;
        ctr_ovl = stretchlim(im2uint8(ovlp_img),[0.001 .999]);
        ctr_odi = stretchlim(odi_img,[0.001 .999]);
        ctr_ori_ctr = stretchlim(ori_ctr_img,[0.001 .999]);
        ctr_ori_ips = stretchlim(ori_ips_img,[0.001 .999]);
        if ~isempty(ROIimgs_ori_bino)
            ctr_ori_bino = stretchlim(ori_bino_img,[0.001 .999]);
        end
        
        
        ovlp_img = imadjust(im2uint8(ovlp_img), ctr_ovl);
        odi_img = imadjust(odi_img, ctr_odi);
        ori_ctr_img = imadjust(ori_ctr_img, ctr_ori_ctr);
        ori_ips_img = imadjust(ori_ips_img, ctr_ori_ips);
        if ~isempty(ROIimgs_ori_bino)
            ori_ips_bino = imadjust(ori_bino_img, ctr_ori_bino);
        end
        %
        %     view_tiff(ovlp_img) ;
        %     view_tiff(odi_img);
        
        if splitmd
            splitfact = 1.5;
            splitfact2 = 1.5;
            splitcol = 2^8/3;
            splitcol2 = 0;
            
            splitsize =  [ round(( size(ovlp_img,2) / length(roi_data)  ) / splitfact)  size(cell_array_ovl{1},1)];
            splitsize2 =  [ round(( size(ovlp_img,2) / length(roi_data)  ) / splitfact2)  size(cell_array_ovl{1},1)  ];
            
            grayimg = repmat(ones(splitsize(2), splitsize(1)),length(srtidx),1,3) * splitcol;
            grayimg2 = repmat(ones(splitsize2(2), splitsize2(1)),length(srtidx),1,3) * splitcol;
            blackimg = repmat(zeros(splitsize2(2), splitsize2(1)),length(srtidx),1,3)* splitcol2 ;
            
            splitpoint = size(ovlp_img,2) / length(roi_data) * baseline
            
            if blck
                if baseline2
                    splitpoint2 = size(ovlp_img,2) / length(roi_data) * (baseline + baseline2)
                    ovlp_img = [ovlp_img(:,1:splitpoint,:) grayimg ovlp_img(:,splitpoint+1:splitpoint2,:) grayimg ovlp_img(:,splitpoint2+1:end,:) blackimg];
                    odi_img = [odi_img(:,1:splitpoint,:) grayimg odi_img(:,splitpoint+1:splitpoint2,:) grayimg odi_img(:,splitpoint2+1:end,:) blackimg];
                    ori_ctr_img = [ori_ctr_img(:,1:splitpoint,:) grayimg ori_ctr_img(:,splitpoint+1:splitpoint2,:) grayimg ori_ctr_img(:,splitpoint2+1:end,:) blackimg];
                    ori_ips_img = [ori_ips_img(:,1:splitpoint,:) grayimg ori_ips_img(:,splitpoint+1:splitpoint2,:) grayimg ori_ips_img(:,splitpoint2+1:end,:) blackimg];
                    if ~isempty(ROIimgs_ori_bino)
                        ori_bino_img = [ori_bino_img(:,1:splitpoint,:) grayimg ori_bino_img(:,splitpoint+1:splitpoint2,:) grayimg ori_bino_img(:,splitpoint2+1:end,:) blackimg];
                    end
                else
                    splitpoint2 = size(ovlp_img,2) / length(roi_data) * (baseline + recovery1)
                    ovlp_img = [ovlp_img(:,1:splitpoint,:) grayimg ovlp_img(:,splitpoint+1:splitpoint2,:)  blackimg];
                    odi_img = [odi_img(:,1:splitpoint,:) grayimg odi_img(:,splitpoint+1:splitpoint2,:) blackimg];
                    ori_ctr_img = [ori_ctr_img(:,1:splitpoint,:) grayimg ori_ctr_img(:,splitpoint+1:splitpoint2,:) blackimg];
                    ori_ips_img = [ori_ips_img(:,1:splitpoint,:) grayimg ori_ips_img(:,splitpoint+1:splitpoint2,:) blackimg];
                    if ~isempty(ROIimgs_ori_bino)
                        ori_bino_img = [ori_bino_img(:,1:splitpoint,:) grayimg ori_bino_img(:,splitpoint+1:splitpoint2,:) blackimg];
                    end
                end
            else
                if baseline2
                    splitpoint2 = size(ovlp_img,2) / length(roi_data) * (baseline + baseline2);
                    ovlp_img = [ovlp_img(:,1:splitpoint,:) grayimg ovlp_img(:,splitpoint+1:splitpoint2,:) grayimg2 ovlp_img(:,splitpoint2+1:end,:) ];
                    odi_img = [odi_img(:,1:splitpoint,:) grayimg odi_img(:,splitpoint+1:splitpoint2,:) grayimg2 odi_img(:,splitpoint2+1:end,:) ];
                    ori_ctr_img = [ori_ctr_img(:,1:splitpoint,:) grayimg ori_ctr_img(:,splitpoint+1:splitpoint2,:) grayimg2 ori_ctr_img(:,splitpoint2+1:end,:) ];
                    ori_ips_img = [ori_ips_img(:,1:splitpoint,:) grayimg ori_ips_img(:,splitpoint+1:splitpoint2,:) grayimg2 ori_ips_img(:,splitpoint2+1:end,:) ];
                    if ~isempty(ROIimgs_ori_bino)
                        ori_bino_img = [ori_bino_img(:,1:splitpoint,:) grayimg ori_bino_img(:,splitpoint+1:splitpoint2,:) grayimg2 ori_bino_img(:,splitpoint2+1:end,:) ];
                    end
                else
                    splitpoint2 = size(ovlp_img,2) / length(roi_data) * (baseline + recovery1);
                    ovlp_img = [ovlp_img(:,1:splitpoint,:) grayimg ovlp_img(:,splitpoint+1:splitpoint2,:)  ];
                    odi_img = [odi_img(:,1:splitpoint,:) grayimg odi_img(:,splitpoint+1:splitpoint2,:) ];
                    ori_ctr_img = [ori_ctr_img(:,1:splitpoint,:) grayimg ori_ctr_img(:,splitpoint+1:splitpoint2,:) ];
                    ori_ips_img = [ori_ips_img(:,1:splitpoint,:) grayimg ori_ips_img(:,splitpoint+1:splitpoint2,:) ];
                    if ~isempty(ROIimgs_ori_bino)
                        ori_bino_img = [ori_bino_img(:,1:splitpoint,:) grayimg ori_bino_img(:,splitpoint+1:splitpoint2,:) ];
                    end
                    
                end
            end
        end
        if isempty(ROIimgs_ori_bino)
            view_tiff([ovlp_img ori_ctr_img ori_ips_img odi_img]);
        else
            view_tiff([ovlp_img ori_ctr_img ori_ips_img ori_bino_img odi_img]);
        end
        view_tiff([ovlp_img   ]);
        view_tiff([ ori_ctr_img  ]);
        view_tiff([  ori_ips_img ]);
        if ~isempty(ROIimgs_ori_bino)
            view_tiff([  ori_bino_img ]);
        end
        view_tiff([   odi_img]);
        
        svnm = fullfile(dirparts{end}(1), dirparts{end}(2), dirparts{end}(3), dirparts{end}(4), [ 'ovlp_img' num2str(exp(1)) '.tif']);
        imwrite(ovlp_img, svnm{1})
        svnm = fullfile(dirparts{end}(1), dirparts{end}(2), dirparts{end}(3), dirparts{end}(4), [ 'odi_img' num2str(exp(1)) '.tif']);
        imwrite(odi_img, svnm{1})
        svnm = fullfile(dirparts{end}(1), dirparts{end}(2), dirparts{end}(3), dirparts{end}(4), [ 'oric_img' num2str(exp(1)) '.tif']);
        imwrite(ori_ctr_img, svnm{1})
        svnm = fullfile(dirparts{end}(1), dirparts{end}(2), dirparts{end}(3), dirparts{end}(4), [ 'orii_img' num2str(exp(1)) '.tif']);
        imwrite(ori_ips_img, svnm{1})
    end
    
    %% Figure 7: grouped PSTH plots over time
    if PSTHs
        % group PSTHs
        cellimg = 1;
        
        %     rec_groups = [resp_ovlpidx_rois{:}]; type = '_resp_ovlpidx_rois'; groups_selected_disp = groups_selected% allways responder
        rec_groups = [resp_ovlpidx_base_all_rois{:}]; type = '_resp_ovlpidx_base_all_rois'; groups_selected_disp = cum_int_base_groups% base responder
        %             rec_groups = [resp_ovlpidx_all_rois{:}]; type = '_resp_ovlpidx_all_rois';groups_selected_disp = cum_int_all_groups %all rois
        if cellimg && ~exist('ROIimgs_red')
            [ROIimgs_red , ROIimgs_green ,ROIimgs_act, ROIimgs_odi, ROIimgs_ovl, ROIimgs_ori_ctr, ROIimgs_ori_ips, ROIimgs_ori_bino] = cell_group_display([],mouse,exp,roi_data,pixelborder, datapath);
        end
        ROIimgs_ovl = ROIimgs_ori_ctr;
        ROIimgs_odi = ROIimgs_ori_ips;
        
        
        for kk = 1:length(groups_selected_disp);
            
            ROId = rec_groups(kk,:)
            tic
            if ~cellimg
                [DFpeak F02  SigF0 delta_max delta_meanF0 delta_sigmaF0 mean_plotdata ]  = PSTH_plot_chronic(roi_data, ana_data, ids_new,ROId,ratio,0,1, neuropilfct, eyenum, baseline, runexclude );
            else
                [DFpeak F02  SigF0 delta_max delta_meanF0 delta_sigmaF0 mean_plotdata ]  = PSTH_plot_chronic(roi_data, ana_data, ids_new,ROId,ratio,0,1, neuropilfct, eyenum, baseline, runexclude, ROIimgs_ovl, ROIimgs_odi);
            end
            %             set(hdl, 'Color', 'w', 'position', [5 49 1272 1481]);
            %         set(gcf, 'position', [20 49 1272 1481]);
            
            
            % format figure
            % Defaults for Cell press. 1 col: 85mm, 1.5 col: 114mm, 2col:174mm
            % Defaults for Nature press. 1 col: 89mm, 1.5 col: 136mm, 2col:183mm
            width = 11.4;                  % Width in cm
            height = width;     % Height in cm (golden ratio default  1/1.618)
            xLeft = (21-width)/2; yTop = (30-height)/2;
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'Color', 'w');
            set(gcf,'PaperPosition',[xLeft yTop width height])
            pos = get(gcf, 'Position');
            set(gcf, 'Position', [pos(1) pos(2)-width*100 width*100, height*100]); %<- Set size
            
            tightfig;
            
            %                 try
            %                     dirname = fullfile(dirparts{end}(1), dirparts{end}(2), dirparts{end}(3), dirparts{end}(4), dirparts{end}(5), 'tuningcurves', num2str(exp(1)));
            %                 catch
            dirname = fullfile(dirparts2{end}(1), dirparts2{end}(2), dirparts2{end}(3), dirparts2{end}(4), dirparts2{end}(5), 'tuningcurves', num2str(exp(1)));
            %                 end
            
            if ~runexclude
                if ~ratio
                    filename = [mouse{1} '_ROIgroup' num2str(groups_selected(kk)) '_site' num2str(exp(1)) type '.pdf'];
                    %                         filename = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '.pdf'];
                    filename2 = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '.png'];
                else
                    filename = [mouse{1} '_ROIgroup' num2str(groups_selected(kk)) '_site' num2str(exp(1))  type '_ratio.pdf'];
                    %                         filename = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_ratio.pdf'];
                    filename2 = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1))  type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_ratio.png'];
                end
            else
                if ~ratio
                    filename = [mouse{1} '_ROIgroup' num2str(groups_selected(kk)) '_site' num2str(exp(1)) type '.pdf'];
                    %                         filename = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_runexclude.pdf'];
                    filename2 = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_runexclude.png'];
                else
                    filename = [mouse{1} '_ROIgroup' num2str(groups_selected(kk)) '_site' num2str(exp(1))  type '_ratio.pdf'];
                    %                         filename = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1)) type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_runexclude_ratio.pdf'];
                    filename2 = [mouse{1} '_ROIgroup' num2str(groups_selected_disp(kk)) '_site' num2str(exp(1))  type '_zth' num2str(z_thresh) '_zf' num2str(z_thresh_fraction) '_runexclude_ratio.png'];
                end
            end
            
            ax=axes('Units','Normal','Position',[.095 .095 .85 .85],'Visible','off');
            
            title(filename2);
            set(get(ax,'Title'),'Visible','on', 'FontSize', 20, 'Interpreter', 'none')
            
            mkdir(dirname{1});
            
            %                 plot2svg(fullfile(dirname{1},filename),gcf); %yeah - vector +
            export_fig(gcf,fullfile(dirname{1},filename), '-q101');
            try
                %                     export_fig(gcf,fullfile(dirname{1},filename), '-q101');
                
                export_fig(gcf,fullfile(dirname{1},filename2));
            end
            
            close(gcf);
            
        end
    end
    
    %% Figure 8: preferred ORI trace PSTH
    % circshift this for all ORIs! Like this we get a nice average tuning
    % curce PSTH display. Cut off afet 2 sec post-peak
    if ~noplot
        figure(23522342)
    end
    for xx = 1: length(roi_data)
        PSTH_cont(:,:,xx) = cell2mat(plotarray_PrefOriContra_PSTH(:,xx))';
        PSTH_ips(:,:,xx) = cell2mat(plotarray_PrefOriIpsi_PSTH(:,xx))';
        
        PSTH_cont_a(:,:,xx) = cell2mat(plotarray_PrefOriContra_PSTH_a(:,xx))';
        PSTH_ips_a(:,:,xx) = cell2mat(plotarray_PrefOriIpsi_PSTH_a(:,xx))';
        
        x =linspace(0,size(PSTH_cont,1)./SamplingFreq1(xx), size(PSTH_cont,1));
        
        if ~noplot && size(PSTH_cont,2) >1
            subplot(2,length(roi_data),xx)
            shadedErrorBar(x,nanmean(PSTH_cont(:,:,xx)'),nanstd(PSTH_cont(:,:,xx)')./sqrt(size(PSTH_cont,2)),'b',1); hold on
            shadedErrorBar(x,nanmean(PSTH_ips(:,:,xx)'),nanstd(PSTH_ips(:,:,xx)')./sqrt(size(PSTH_cont,2)),'r',1);
            
            xlabel('responders only (fulidx_z, matched) [s]')
            title(['Session ' num2str(xx)])
            ylim([0 1000])
            xlim([x(1) x(end)])
            
            
            subplot(2,length(roi_data),xx+length(roi_data))
            shadedErrorBar(x,nanmean(PSTH_cont_a(:,:,xx)'),nanstd(PSTH_cont_a(:,:,xx)')./sqrt(size(PSTH_cont_a,2)),'b',1); hold on
            shadedErrorBar(x,nanmean(PSTH_ips_a(:,:,xx)'),nanstd(PSTH_ips_a(:,:,xx)')./sqrt(size(PSTH_cont_a,2)),'r',1);
            xlabel('all (fullidx, matched) [s]')
            title(['Session ' num2str(xx)])
            ylim([0 500])
            xlim([x(1) x(end)])
        end
    end
    
    try
        if ~noplot
            tightfig
        end
    end
    
    %% Figure 9: preferred ORI centered FIT (I split this and the previous for better overview)
    if ~noplot
        figure(235222)
    end
    try
        for xx = 1: length(roi_data)
            for cc = 1:size(plotarray_PrefDirContra,1);
                xdc = circshift(cell2mat(plotarray_360Dir_ctr_tuning(cc,xx)), [0 -round(plotarray_PrefDirContra(cc,xx))+180]);
                xdi = circshift(cell2mat(plotarray_360Dir_ips_tuning(cc,xx)), [0 -round(plotarray_PrefDirIpsi(cc,xx))+180]);
                Tune_cont(cc,:,xx) = xdc;
                Tune_ips(cc,:,xx) = xdi;
            end
            x = 1:size(Tune_cont,2); x = x-180;
            if ~noplot
                subplot(3,length(roi_data),xx);
                shadedErrorBar(x,nanmean(Tune_cont(:,:,xx)),nanstd(Tune_cont(:,:,xx))./sqrt(size(Tune_cont,1)),'b',1); hold on
                shadedErrorBar(x,nanmean(Tune_ips(:,:,xx)),nanstd(Tune_ips(:,:,xx))./sqrt(size(Tune_ips,1)),'r',1);
                ylim([-20 500]);
                xlabel('Preferred direction [in {\circ}]');
                title(['Session ' num2str(xx)]);
            end
        end
        
        for xx = 1: length(roi_data)
            for cc = 1:size(plotarray_PrefDirContra,1)
                
                xdc = circshift(cell2mat(plotarray_360Dir_ctr_tuning(cc,xx)), [0 -round(plotarray_OpPrefDirContra(cc,xx))+180]);
                xdi = circshift(cell2mat(plotarray_360Dir_ips_tuning(cc,xx)), [0 -round(plotarray_OpPrefDirIpsi(cc,xx))+180]);
                Tune_cont(cc,:,xx) = xdc;
                Tune_ips(cc,:,xx) = xdi;
            end
            x = 1:size(Tune_cont,2);x = x-180;
            if ~noplot
                subplot(3,length(roi_data),xx+length(roi_data))
                shadedErrorBar(x,nanmean(Tune_cont(:,:,xx)),nanstd(Tune_cont(:,:,xx))./sqrt(size(Tune_cont,1)),'b',1); hold on
                shadedErrorBar(x,nanmean(Tune_ips(:,:,xx)),nanstd(Tune_ips(:,:,xx))./sqrt(size(Tune_ips,1)),'r',1);
                ylim([-20 500])
                xlabel('Opposite direction [in {\circ}]')
                %      title('Average tuning curves - centered on opposite direction')
                title(['Session ' num2str(xx)])
            end
        end
        
        for xx = 1: length(roi_data)
            for cc = 1:size(plotarray_PrefDirContra,1)
                
                xdc = circshift(cell2mat(plotarray_180Ori_ctr_tuning(cc,xx)), [0 -round(plotarray_PrefOriContra(cc,xx))+90]);
                xdi = circshift(cell2mat(plotarray_180Ori_ips_tuning(cc,xx)), [0 -round(plotarray_PrefOriIpsi(cc,xx))+90]);
                Tune_contO(cc,:,xx) = xdc;
                Tune_ipsO(cc,:,xx) = xdi;
            end
            x = 1:size(Tune_contO,2); x= x-90;
            if ~noplot
                subplot(3,length(roi_data),xx+length(roi_data)*2)
                shadedErrorBar(x,nanmean(Tune_contO(:,:,xx)),nanstd(Tune_contO(:,:,xx))./sqrt(size(Tune_contO,1)),'b',1); hold on
                shadedErrorBar(x,nanmean(Tune_ipsO(:,:,xx)),nanstd(Tune_ipsO(:,:,xx))./sqrt(size(Tune_ipsO,1)),'r',1);
                ylim([-20 1000])
                xlabel('Preferred orientation [in {\circ}]')
                %      title('Average tuning curves - centered on preferred orientation')
                title(['Session ' num2str(xx)])
            end
        end
        if ~noplot
            tightfig
        end
        title('Direction Tuning-Curve averages. Both centered on preferred direction')
    end
    %% Figure 10: ROIbased colormaps
    if roimaps
        xreso = 1024;
        yreso = 1024;
        
        if collapse
            xreso =xreso/2;
            yreso = yreso/2;
        end
        
        img = zeros(xreso,yreso);
        
        colorlevels = 11;
        
        coc = cbrewer('div', 'RdBu', colorlevels); %toned down bipolar maps... use for print
        
        %     coc = makeColorMap([1 0 0], [1 1 1], [0 0 1], colorlevels) % full
        % red/green RGB range. Use for screen?
        
        intmap = ColorMapInt(coc, colorlevels); % generate linear scaled intensity map
        ODI_HSV = ones( xreso, yreso, 3, length(roi_data) ).*.5;
        ODI_HSVns = ones( xreso, yreso, 3, length(roi_data) ).*.5;
        temp = zeros( xreso, yreso);
        
        % image border
        border = 0;
        boldroi = 2;
        
        ODI_HSV(1:yreso,1:border,:) = 0;
        ODI_HSV(1:border ,1:xreso,:) = 0;
        
        ODI_HSVns(1:yreso,1:border,:) = 0;
        ODI_HSVns(1:border ,1:xreso,:) = 0;
        
        %rescale ODIs to 0->1 for display purposes (indexing)
        if plotarray_OD_map
            ODIs2 = (plotarray_OD +1)/2;
            ODIs2(isnan(ODIs2)) = 0;
            ODIs2(ODIs2>1) = 1;
            ODIs2(ODIs2<-1) = -1;
            Ampli = (plotarray_AmpContra + plotarray_AmpIpsi)./2;
            Ampli(isnan(Ampli)) = 0;
        elseif baselinerois_all_map
            ODIs2 =  (all_n_n2_ODI_morph_nan +1)/2;
            %         ODIs2(isnan(ODIs2)) = 0;
            ODIs2(ODIs2>1) = 1;
            ODIs2(ODIs2<-1) = -1;
            Ampli = (cell2mat(n_n2_amp_contra_morph')' + cell2mat(n_n2_amp_ipsi_morph')')./2;
            Ampli(isnan(Ampli)) = 0;
        elseif baselinerois_map
            ODIs2 = (plotarray_base_n_OD +1)/2;
            ODIs2(isnan(ODIs2)) = 0;
            ODIs2(ODIs2>1) = 1;
            ODIs2(ODIs2<-1) = -1;
            Ampli = (plotarray_base_n_AmpContra + plotarray_base_n_AmpIpsi)./2;
            Ampli(isnan(Ampli)) = 0;
        elseif baselinerois_md_all_map;
            %         all_base_md_ODI_morph_nan
            ODIs2 =  (all_base_md_ODI_morph_nan +1)/2;
            %         ODIs2(isnan(ODIs2)) = 0;
            ODIs2(ODIs2>1) = 1;
            ODIs2(ODIs2<-1) = -1;
            Ampli = (cell2mat(base_md_amp_contra_morph')' + cell2mat(base_md_amp_ipsi_morph')')./2;
            Ampli(isnan(Ampli)) = 0;
        end
        
        %set the Ampli cutoffs
        mp = prctile(Ampli(:),90);
        minp = prctile(Ampli(:),10);
        Ampli = Ampli-minp;
        Ampli = Ampli/(mp-minp);
        
        Ampli(Ampli>1)=1;
        Ampli(Ampli<0)=0;
        
        if collapse
            xreso = xreso*2;
            yreso = yreso*2;
        end
        
        kkp = 0;
        for i = 1:length(roi_data);
            if plotarray_OD_map
                selidx = resp_ovlpidx_rois{i};
                kkp = kkp + 1;
            elseif baselinerois_map
                if i==baseline_pair(1) | i==baseline_pair(2)
                    selidx = resp_ovlpidx_base_n_rois{i};
                    kkp = kkp + 1;
                else
                    continue
                end
            elseif baselinerois_all_map
                if i==baseline-2 | i== baseline- 1 | i==baseline
                    kkp = kkp + 1;
                    selidx = n_n2_morph_rois_matched{kkp};
                    
                else
                    continue
                end
                
            elseif baselinerois_md_all_map
                if i==baseline | i== baseline+1
                    kkp = kkp + 1;
                    selidx = base_md_morph_rois_matched{kkp};
                    
                else
                    continue
                end
            end
            
            for k = 1:length(roi_data{i}.ROIs(selidx));
                if ~isnan(ODIs2(k,kkp))
                    [ROIx ROIy]= ind2sub([xreso yreso],roi_data{i}.ROIs(selidx(k)).indices);
                else
                    
                    temp(roi_data{i}.ROIs(selidx(k)).indices) = 1;
                    bw_img_perim = bwperim(temp);
                    se = strel('square',boldroi);
                    bw_img_perim = imdilate(bw_img_perim,se);
                    perim_indices = find(bw_img_perim==1);
                    [ROIx ROIy]= ind2sub([xreso yreso],perim_indices);
                    temp = zeros( xreso, yreso);
                end
                if collapse
                    if min(ROIx)>xreso/2
                        ROIx = ROIx-xreso/2;
                        ROIx(ROIx<1)=1;
                    end
                    if min(ROIy)>yreso/2
                        ROIy = ROIy-yreso/2;
                        ROIy(ROIy<1)=1;
                    end
                end
                for p = 1:length(ROIx)
                    if ~isnan(ODIs2(k,kkp))
                        try
                            ODI_HSV(ROIx(p), ROIy(p), :,kkp) = intmap{1+round(Ampli(k,kkp)*(colorlevels-1))}(1+round(ODIs2(k,kkp)*(colorlevels-1)),:);
                            ODI_HSVns(ROIx(p), ROIy(p), :,kkp)  = intmap{colorlevels-1}(1+round(ODIs2(k,kkp)*(colorlevels-1)),:);
                        catch
                            disp('WTF1?');
                        end
                    else
                        try
                            ODI_HSV(ROIx(p), ROIy(p), :,kkp) = 1;
                            ODI_HSVns(ROIx(p), ROIy(p), :,kkp)  = 0;
                        catch
                            disp('WTF2?');
                        end
                    end
                end
            end
            
            
            try
                odi_txt_img = ~text2im([mouse{1} '_exp' num2str(exp(i))]).*1;
                %             odi_txt_img = ~text2im([mouse{1} ' d ' num2str(round(TB(i)))]).*1;%
                %                  odi_txt_img_b(:,find(odi_txt_img(:,:) == 0)) = 0.5;
            catch
                odi_txt_img = ~text2im([mouse '_exp' num2str(exp(i))]).*1;
                %             odi_txt_img = ~text2im([mouse '_' yearmonthday(i,:) '_exp' num2str(exp(i))]).*1;
            end
            odi_txt_img=imresize(odi_txt_img,1.1);
            %             odi_txt_img_b=imresize(odi_txt_img_b,1.1);
            
            %         ODI_HSV(:,:,:,kkp) = implace(ODI_HSV(:,:,:,kkp),odi_txt_img,10,10);
            %         ODI_HSVns(:,:,:,kkp) = implace(ODI_HSVns(:,:,:,kkp),odi_txt_img,10,10);
        end
        
        
        try
            %         montageTR(ODI_HSV);
            %         if ~baselinerois_map
            %             view_tiff(ODI_HSV(:,:,:,baseline));
            %             view_tiff(ODI_HSV(:,:,:,baseline+1));
            %             view_tiff(ODI_HSV(:,:,:,baseline+recovery1));
            %         else
            %             view_tiff(ODI_HSV(:,:,:, baseline_pair(1)));
            %             view_tiff(ODI_HSV(:,:,:, baseline_pair(2)));
            %         end
            
        end
        try
            montageTR(ODI_HSVns);
            if plotarray_OD_map
                view_tiff(ODI_HSVns(:,:,:,baseline));
                view_tiff(ODI_HSVns(:,:,:,baseline+1));
                view_tiff(ODI_HSVns(:,:,:,baseline+recovery1));
            elseif baselinerois_map
                view_tiff(ODI_HSVns(:,:,:, baseline_pair(1)));
                view_tiff(ODI_HSVns(:,:,:, baseline_pair(2)));
            elseif baselinerois_all_map
                view_tiff(ODI_HSVns(:,:,:,1));
                view_tiff(ODI_HSVns(:,:,:,2));
                view_tiff(ODI_HSVns(:,:,:,3));
            elseif baselinerois_md_all_map
                view_tiff(ODI_HSVns(:,:,:,1));
                view_tiff(ODI_HSVns(:,:,:,2));
            end
        end
        
    end
    
    MouseID = mouse;
    mouse = MouseID;
    
    % - - - -- - - - - - - - - - - - - - - - - - - - -- - - -- --- - -
    
    %% [+] Figure 11: plot single-cell ODI over days plotarray_OD
    if ~noplot
        figure(2342)
        plot(TB,plotarray_OD', 'k'); hold on
        plot(TB,nanmean(plotarray_OD',2), 'r', 'LineWidth', 3)
        plot(TB,nanmedian(plotarray_OD',2), ':r', 'LineWidth', 3)
        
        hleg1 = legend('mean (plotarray_OD )','median (plotarray_OD)');
        
        a = title([MouseID ' :  ODscore (groupsize: ' num2str(max([roi_data{i}.ROIs(:).groupsize])) 'TPs). ' num2str(size(groups_selected)) ' cells in all TPs > 8z in baseline (resp_ovlpidx)']);
        set(a, 'Interpreter', 'none')
        xlabel('Days after MD onset')
        ylabel(' ODI (Fit): ( Contra-Ipsi ) / ( Contra + Ipsi )');
        
        
        if ~isnan(baseline2)
            colorarray = [repmat([0 0 0],baseline,1) ; [0.5 0 0] ; repmat([0 0 0],recovery1,1) ; [1 0 0] ; repmat([0 0 0],baseline2,1)]
        else
            colorarray = [repmat([0 0 0],baseline,1) ; [0.5 0 0] ; repmat([0 0 0], recovery1,1)]
        end
        
        colorarray = repmat([0 0 0],length(plotarray_OD),1)
        
        %% Figure 12: plot single-cell ODI distributions over days both population-based ad single cell based
        colorarray = num2cell(colorarray',1);
        
        figure(461771)
        
        subplot(3,1,1);
        distributionPlot([plotarray_OD],'histOpt',1, 'addSpread',0, 'showMM',5, 'color', colorarray, 'histOri' , 'right');  hold off
        hline(0);hline(1);hline(-1);ylim([-1 1]);
        xlabel('Refound responder single cell ODI (plotarray_OD - resp_ovlpidx)')
        
        subplot(3,1,2);
        distributionPlot([plotarray_OD_pop],'histOpt',1, 'addSpread',0, 'showMM',5, 'color', colorarray,  'histOri' , 'right');  hold off
        hline(0);hline(1);hline(-1);ylim([-1 1]);
        xlabel('All responder single cell ODI (plotarray_OD_pop - fullidx_z  ')
        % distributionPlot([plotarray_OD_amp_contra],'histOpt',1, 'addSpread',1, 'showMM',6, 'colormap', gray);  hold off
        % distributionPlot([plotarray_OD_amp_ipsi],'histOpt',1, 'addSpread',1, 'showMM',6, 'colormap', gray);  hold off
    end
    for i = 1:length(ana_data)
        plotarray_amp_pop(i) = ( nanmean(plotarray_OD_amp_contra{i}) - nanmean(plotarray_OD_amp_ipsi{i}) ) ./ ( nanmean(plotarray_OD_amp_contra{i}) + nanmean(plotarray_OD_amp_ipsi{i}) );
        plotarray_amp_pop_ref(i) = ( nanmean(plotarray_AmpContra(:,i)) - nanmean(plotarray_AmpIpsi(:,i)) ) ./ ( nanmean(plotarray_AmpContra(:,i)) + nanmean(plotarray_AmpIpsi(:,i)) );
        plotarray_amp_pop_all(i) = ( nanmean(plotarray_OD_amp_contra_all(:,i)) - nanmean(plotarray_OD_amp_ipsi_all(:,i)) ) ./ ( nanmean(plotarray_OD_amp_contra_all(:,i)) + nanmean(plotarray_OD_amp_ipsi_all(:,i)) );
    end
    if ~noplot
        subplot(3,1,3);
        
        % distributionPlot([plotarray_OD_amp_pop],'histOpt',1, 'addSpread',0, 'showMM',5, 'color', colorarray, 'histOri' , 'right');  hold off
        plot(TB,plotarray_amp_pop,'-ok'); hold on
        plot(TB,nanmean(plotarray_OD',2), '-or');
        plot(TB,nanmean(plotarray_amp_pop_ref',2), ':or');
        
        
        % plot(nanmean([plotarray_OD_pop{:}]',2), '-ob');
        
        legend('non-matched cross-sectional responder pop. amplitude average OD (ampli ave respidx_z)','Refound responder single cell ODI (ODI  ave resp_ovlpidx)', 'Refound responder single cell ODI (ampli ave resp_ovlpidx)')
        hline(0); vline(baseline)
        
        
        figure(4234)
        hold on
        plot(TB,plotarray_amp_pop,'-*k'); hold on
        plot(TB,nanmean(plotarray_OD',2), '-*r');
        plot(TB,nanmean(plotarray_amp_pop_all',2), '-*b');
        
        % hline(0);hline(1);hline(-1);ylim([-1 1]);
        %% Figure 13: plot amplitude over days
        figure(23242)
        try
            shadedErrorBar(TB,nanmean(plotarray_AmpContra',2), nanstd(plotarray_AmpContra)'./sqrt(size(plotarray_AmpContra,1)),'b',1); hold on
            %      plot(TB,nanmedian(plotarray_AmpContra',2), '--b', 'LineWidth', 3); hold on
            shadedErrorBar(TB,nanmean(plotarray_AmpIpsi',2), nanstd(plotarray_AmpIpsi)'./sqrt(size(plotarray_AmpIpsi,1)),'r',1);
            %      plot(TB,nanmedian(plotarray_AmpIpsi',2), '--r', 'LineWidth', 3)
            
            vline([TB(baseline)  TB(baseline+1) TB(baseline+recovery1)]);
            
            a = title([MouseID 'Amplitude refound responders (groupsize: ' num2str(max([roi_data{i}.ROIs(:).groupsize])) 'TPs). ' num2str(size(groups_selected)) ' cells in all TPs > 8z in baseline ']);
            set(a, 'Interpreter', 'none')
            xlabel('Days after MD onset')
            ylabel('\DeltaR/R');
        end
    end
    
    %     %% plotoverall z-score over time
    %     figure(231242)
    %     % plot(TB,plotarray_AmpContra', ':g'); hold on
    %     % % plot(TB,plotarray_AmpIpsi', ':r'); hold on
    %     %
    %     % crit1 = mean(plotarray_OD(:,1:baseline),2)<0.3;
    %     % crit2 = mean(plotarray_OD(:,1:baseline),2)>0.3;
    %
    %     plot(TB,plotarray_peak_zscore, '-k', 'LineWidth', 3);
    %     vline([TB(baseline) TB(baseline+recovery1)]);
    %
    %     a = title([MouseID ' :  peak z-score (ipsi or contra) over simulus baseline' ]);
    %     set(a, 'Interpreter', 'none')
    %     xlabel('Days after MD onset')
    %     ylabel('peak z-score over stimulus baseline');
    
    %% Data extraction - correlations - shared variables
    try
        baselineOD = nanmean(plotarray_OD(:,baseline-2:baseline),2);
        baselineAMPc = nanmean(plotarray_AmpContra(:,baseline-2:baseline),2);
        baselineAMPi = nanmean(plotarray_AmpIpsi(:,baseline-2:baseline),2);
    catch
        baselineOD = nanmean(plotarray_OD(:,baseline-1:baseline),2);
        baselineAMPc = nanmean(plotarray_AmpContra(:,baseline-2:baseline),2);
        baselineAMPi = nanmean(plotarray_AmpIpsi(:,baseline-2:baseline),2);
    end
    
    baselineOD_last = plotarray_OD(:,baseline);
    MDOD = plotarray_OD(:,baseline+1);
    deltaMD1 = plotarray_OD(:,baseline+1) - baselineOD;
    RecOD = plotarray_OD(:,baseline+recovery1);
    
    baselineAMPc_last = plotarray_AmpContra(:,baseline);
    MDAMPc = plotarray_AmpContra(:,baseline+1);
    deltaMD1_AMPc = plotarray_AmpContra(:,baseline+1) - baselineAMPc;
    RecAMPc = plotarray_AmpContra(:,baseline+recovery1);
    
    baselineAMPi_last = plotarray_AmpIpsi(:,baseline);
    MDAMPi = plotarray_AmpIpsi(:,baseline+1);
    deltaMD1_AMPi = plotarray_AmpIpsi(:,baseline+1) - baselineAMPi;
    RecAMPi = plotarray_AmpIpsi(:,baseline+recovery1);
    
    
    try
        STDbaselineOD = nanstd(plotarray_OD(:,baseline-2:baseline)')';
        STDbaselineAMPc = nanstd(plotarray_AmpContra(:,baseline-2:baseline)')';
        STDbaselineAMPi = nanstd(plotarray_AmpIpsi(:,baseline-2:baseline)')';
    catch
        STDbaselineOD = nanstd(plotarray_OD(:,baseline-1:baseline)')';
        STDbaselineAMPc = nanstd(plotarray_AmpContra(:,baseline-1:baseline)')';
        STDbaselineAMPi = nanstd(plotarray_AmpIpsi(:,baseline-1:baseline)')';
    end
    try
        n4baselineOD = plotarray_OD(:,[baseline_pair(1) baseline_pair(2)]);
        n4baselineOD2 = plotarray_OD(:,[baseline-1 baseline]);
        n4baselineOD3 = plotarray_OD(:,[baseline-2 baseline]);
    catch
        n4baselineOD = [];
        n4baselineOD2 = [];
        n4baselineOD3 = [];
    end
    try
        n14baselineOD = plotarray_OD(:,[baseline_pair14(1) baseline_pair14(2)]);
    catch
        n14baselineOD = [];
    end
    %% Figure 14: mean baseline ODI vs. MD1 change
    % deltaMDampIpsiMD1 = plotarray_AmpIpsi(:,baseline) - plotarray_AmpIpsi(:,baseline+1);
    % deltaMDampContraMD1 = plotarray_AmpContra(:,baseline) - plotarray_AmpContra(:,baseline+1);
    
    % crit = plotarray_OD(:,baseline)<0;
    % crit = logical(ones(size(baselineOD)));
    crit = deltaMD1>1*STDbaselineOD;
    % crit = mean(plotarray_OD(:,1:baseline),2)>0;
    % crit = mean(plotarray_OD(:,1:baseline),2)>0 & mean(plotarray_OD(:,1:baseline),2)<0.9;
    
    % crit = sum(plotarray_AmpContra(:,1:3)<prctile(plotarray_AmpContra(:),50),2)>1;
    % crit = sum(plotarray_AmpIpsi(:,1:3)>prctile(plotarray_AmpIpsi(:),80),2)>1;
    try
        if ~noplot
            figure(353);
            scatter(baselineOD(crit), deltaMD1(crit), 'ok', 'filled'); hold on;
            xlabel('mean baseline ODI');
            ylabel('\DeltaMD1');
            [p r] = corrcoef(baselineOD(crit), deltaMD1(crit));
            text(0.5,0.5, ['r = ' num2str(p(2)) ', p = ' num2str(r(2))])
            
            line([-1 1],[1 -1],'LineStyle', '--', 'Color', 'k')
            
            vline(0, ':k');
            hline(0, ':k');
            
            
            title(' Mean baseline ODI vs. DeltaMD1 ');
        end
    end
    % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
    % set(k, 'Interpreter', 'none');
    
    %% Figure 15:  plot corellations: SD baseline ODI vs. MD1 change
    
    % figure(2426257)
    
    %
    % deltaMDampIpsiMD1 = plotarray_AmpIpsi(:,baseline) - plotarray_AmpIpsi(:,baseline+1);
    % deltaMDampContraMD1 = plotarray_AmpContra(:,baseline) - plotarray_AmpContra(:,baseline+1);
    
    % crit = plotarray_OD(:,baseline)<0;
    crit = logical(ones(size(baselineOD)));
    % crit = mean(plotarray_OD(:,1:baseline),2)>0;
    % crit = mean(plotarray_OD(:,1:baseline),2)>0 & mean(plotarray_OD(:,1:baseline),2)<0.9;
    
    % crit = sum(plotarray_AmpContra(:,1:3)<prctile(plotarray_AmpContra(:),50),2)>1;
    % crit = sum(plotarray_AmpIpsi(:,1:3)>prctile(plotarray_AmpIpsi(:),80),2)>1;
    if ~noplot
        figure(35253);
        scatter(STDbaselineOD(crit), deltaMD1(crit), 'ok', 'filled'); hold on;
        xlabel('SD baseline ODI');
        ylabel('\DeltaMD1');
        [p r] = corrcoef(STDbaselineOD(crit), deltaMD1(crit));
        try
            text(0.5,0.5, ['r = ' num2str(p(2)) ', p = ' num2str(r(2))])
        end
        % line([0 1],[0 1],'LineStyle', '--', 'Color', 'k')
        %
        % vline(0, ':k');
        % hline(0, ':k');
        
        
        title(' SD baseline ODI vs. DeltaMD1 ');
        hold off
        % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
        % set(k, 'Interpreter', 'none');
        
        %% plot corellations: baseline ODI vs. full recovery ODI
        
        % figure(2426257)
        %     baselineOD = mean(plotarray_OD(:,1:baseline),2);
        % baselineOD = mean(plotarray_OD(:,baseline:baseline),2);
        
        % STDbaselineOD = std(plotarray_OD(:,1:baseline)')';
        %
        % deltaMDampIpsiMD1 = plotarray_AmpIpsi(:,baseline) - plotarray_AmpIpsi(:,baseline+1);
        % deltaMDampContraMD1 = plotarray_AmpContra(:,baseline) - plotarray_AmpContra(:,baseline+1);
        
        % crit = plotarray_OD(:,baseline)<0;
        %         crit = logical(ones(size(baselineOD)));
        % crit = mean(plotarray_OD(:,1:baseline),2)>0;
        % crit = mean(plotarray_OD(:,1:baseline),2)>0 & mean(plotarray_OD(:,1:baseline),2)<0.9;
        
        % crit = sum(plotarray_AmpContra(:,1:3)<prctile(plotarray_AmpContra(:),50),2)>1;
        % crit = sum(plotarray_AmpIpsi(:,1:3)>prctile(plotarray_AmpIpsi(:),80),2)>1;
    end
    try
        if ~noplot
            figure(15253);
            crit = abs(deltaMD1)>1*STDbaselineOD;
            scatter(baselineOD(crit), RecOD(crit), 'ok', 'filled'); hold on;
            xlabel('mean baseline ODI');
            ylabel('Recovery ODI');
            [p r] = corrcoef(baselineOD(crit), RecOD(crit));
            try
                text(0.5,0.5, ['r = ' num2str(p(2)) ', p = ' num2str(r(2))])
            end
            line([1 -1],[1 -1],'LineStyle', '--', 'Color', 'k')
            %
            vline(0, ':k');
            hline(0, ':k');
            
            
            title('baseline ODI vs. full recovery ODI');
            hold off
            % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
            % set(k, 'Interpreter', 'none');
        end
    end
    if ~noplot
        %% baseline DeltaMD1 vs.DeltaMD2
        if baseline2
            
            %
            %         deltaMD1_2 =  plotarray_OD(:,baseline+1) - plotarray_OD(:,baseline);
            deltaMD2 = plotarray_OD(:,baseline+baseline2+1) - plotarray_OD(:,baseline+baseline2);
            
            deltaMDampIpsiMD1 =  plotarray_AmpIpsi(:,baseline+1) - plotarray_AmpIpsi(:,baseline);
            deltaMDampContraMD1 =  plotarray_AmpContra(:,baseline+1) - plotarray_AmpContra(:,baseline);
            deltaMDampIpsiMD2 =  plotarray_AmpIpsi(:,baseline+baseline2+1) - plotarray_AmpIpsi(:,baseline+baseline2);
            deltaMDampContraMD2 =  plotarray_AmpContra(:,baseline+baseline2+1) - plotarray_AmpContra(:,baseline+baseline2);
            
            % crit = plotarray_OD(:,baseline)<0;
            %          crit = mean(plotarray_OD(:,1:baseline),2)>0;
            %                 crit = abs(deltaMD1_2)>abs(STDbaselineOD);
            crit = abs(deltaMD1)>abs(STDbaselineOD);
            %                 crit = logical(ones(size(deltaMD1_2)));
            
            % crit = mean(plotarray_OD(:,1:baseline),2)>0.1 & mean(plotarray_OD(:,1:baseline),2)<0.6;
            
            %        crit = sum(plotarray_AmpContra(:,1:3)<prctile(plotarray_AmpContra(:),50),2)>1;
            % crit = sum(plotarray_AmpIpsi(:,1:3)>prctile(plotarray_AmpIpsi(:),80),2)>1;
            
            figure(352563);
            %         scatter(deltaMD1_2(crit), deltaMD2(crit), 'ok', 'filled')
            scatter(deltaMD1(crit), deltaMD2(crit), 'ok', 'filled')
            
            xlabel('\DeltaMD1');
            ylabel('\DeltaMD2');
            
            %         [p r] = corrcoef(deltaMD1_2(crit), deltaMD2(crit));
            [p r] = corrcoef(deltaMD1(crit), deltaMD2(crit));
            try
                text(0.5,0.5, ['r = ' num2str(p(2)) ', p = ' num2str(r(2))])
            end
            % line([1 -1],[1 -1],'LineStyle', '--', 'Color', 'k')
            
            line([-1 1],[-1 1],'LineStyle', '--', 'Color', 'k')
            vline(0, ':k');
            hline(0, ':k');
            
            title('DeltaMD1 vs.DeltaMD2');
            % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
            % set(k, 'Interpreter', 'none');
        end
    end
    
    %% Figure 16: plot Tuning stability over 2 bs sessions
    if ~noplot
        figure(3513);
        
        ODn1 = plotarray_base_n_OD(:,baseline_pair(2));
        ODn = plotarray_base_n_OD(:,baseline_pair(1));
        
        crit1 = plotarray_base_n_AmpContra(:,baseline_pair(1))<prctile(plotarray_base_n_AmpContra(:,baseline_pair(1)),75) & plotarray_base_n_AmpContra(:,baseline_pair(1))>prctile(plotarray_base_n_AmpContra(:,baseline_pair(1)),25);
        crit2 = plotarray_base_n_AmpContra(:,baseline_pair(1))>prctile(plotarray_base_n_AmpContra(:,baseline_pair(1)),75); %& plotarray_AmpContra(:,1)<prctile(plotarray_AmpContra(:,1),100);
        crit3 = plotarray_base_n_AmpContra(:,baseline_pair(1))<prctile(plotarray_base_n_AmpContra(:,baseline_pair(1)),25); %& plotarray_AmpContra(:,1)<prctile(plotarray_AmpContra(:,1),100);
        % crit2 = mean(plotarray_OD(:,1:2),2)<0.4;
        % crit2 = plotarray_AmpIpsi(:,1)<prctile(plotarray_AmpIpsi(:,1),10);
        
        
        scatter(ODn1, ODn, 'ok', 'filled'); hold on
        scatter(ODn1(crit2), ODn(crit2), 'or', 'filled')
        scatter(ODn1(crit3), ODn(crit3), 'ob', 'filled')
        lsline
        xlabel('ODI N ');
        ylabel('ODI N + 1 ');
        
        [p r] = corrcoef(ODn1 , ODn );
        [p1 r1] = corrcoef(ODn1 , ODn )
        [p2 r2] = corrcoef(ODn1(crit2) , ODn(crit2) );
        [p3 r3] = corrcoef(ODn1(crit3) , ODn(crit3) );
        
        try
            text(0.5,0.5, ['Full r = ' num2str(p(2)) ', p = ' num2str(r(2))])
            text(0.34,0.34, ['Medium_{50%ile}r = ' num2str(p1(2)) ', p = ' num2str(r1(2))], 'Color' ,'r')
            text(0.25,0.25, ['High_{25%ile}r = ' num2str(p2(2)) ', p = ' num2str(r2(2))], 'Color' ,'r')
            text(0.16,0.16, ['Low_{25%ile}r = ' num2str(p3(2)) ', p = ' num2str(r3(2))], 'Color' ,'b')
        end
        
        line([1 -1],[1 -1],'LineStyle', '--', 'Color', 'k')
        
        
        vline(0, ':k');
        hline(0, ':k');
        
        k = title('plotarray_AmpContra(:,1)>prctile(plotarray_AmpContra(:,1),90)');
        set(k, 'Interpreter', 'none');
        % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
        % set(k, 'Interpreter', 'none');
        hold off
    end
    
    %% Figure 17: plot Tuning stability over 2 bs sessions  14 days apart
    if ~isempty(baseline_pair14)
        if ~noplot
            figure(235113);
            
            ODn1 = plotarray_base_n14_OD(:,baseline_pair14(2));
            ODn = plotarray_base_n14_OD(:,baseline_pair14(1));
            
            crit1 = plotarray_base_n14_AmpContra(:,baseline_pair14(1))<prctile(plotarray_base_n14_AmpContra(:,baseline_pair14(1)),75) & plotarray_base_n14_AmpContra(:,baseline_pair14(1))>prctile(plotarray_base_n14_AmpContra(:,baseline_pair14(1)),25);
            crit2 = plotarray_base_n14_AmpContra(:,baseline_pair14(1))>prctile(plotarray_base_n14_AmpContra(:,baseline_pair14(1)),75); %& plotarray_AmpContra(:,1)<prctile(plotarray_AmpContra(:,1),100);
            crit3 = plotarray_base_n14_AmpContra(:,baseline_pair14(1))<prctile(plotarray_base_n14_AmpContra(:,baseline_pair14(1)),25); %& plotarray_AmpContra(:,1)<prctile(plotarray_AmpContra(:,1),100);
            % crit2 = mean(plotarray_OD(:,1:2),2)<0.4;
            % crit2 = plotarray_AmpIpsi(:,1)<prctile(plotarray_AmpIpsi(:,1),10);
            
            
            scatter(ODn1, ODn, 'ok', 'filled'); hold on
            scatter(ODn1(crit2), ODn(crit2), 'or', 'filled')
            scatter(ODn1(crit3), ODn(crit3), 'ob', 'filled')
            lsline
            xlabel('ODI N(14d) ');
            ylabel('ODI N(14d) + 1 ');
            
            [p r] = corrcoef(ODn1 , ODn );
            [p1 r1] = corrcoef(ODn1 , ODn )
            [p2 r2] = corrcoef(ODn1(crit2) , ODn(crit2) );
            [p3 r3] = corrcoef(ODn1(crit3) , ODn(crit3) );
            
            text(0.5,0.5, ['Full r = ' num2str(p(2)) ', p = ' num2str(r(2))])
            text(0.34,0.34, ['Medium_{50%ile}r = ' num2str(p1(2)) ', p = ' num2str(r1(2))], 'Color' ,'r')
            text(0.25,0.25, ['High_{25%ile}r = ' num2str(p2(2)) ', p = ' num2str(r2(2))], 'Color' ,'r')
            text(0.16,0.16, ['Low_{25%ile}r = ' num2str(p3(2)) ', p = ' num2str(r3(2))], 'Color' ,'b')
            
            
            line([1 -1],[1 -1],'LineStyle', '--', 'Color', 'k')
            
            
            vline(0, ':k');
            hline(0, ':k');
            
            k = title('14d - plotarray_AmpContra(:,1)>prctile(plotarray_AmpContra(:,1),90)');
            set(k, 'Interpreter', 'none');
            % k = title('crit = mean(plotarray_OD(:,1:baseline),2)>0;');
            % set(k, 'Interpreter', 'none');
            hold off
        end
    end
end
% end
%% - - - - DATA SAVING: save relevant variables and backup ALL scripts used (except ana_data and roi_data

% save([adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) '_summary.mat'], '-regexp', '^(?!(roi_data|ana_data|auxdata)$). ');
% save([adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) '_summary.mat'], '-regexp', '^(?!(roi_data||auxdata)$).'); %no roi_data!
if ~fixfit || savedata
    
    %     clear roi_data
    if ~runexclude
        save([adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) '_summary.mat'], '-regexp', '^(?!(auxdata)$).');
    else
        save([adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) '_runexclude_summary.mat'], '-regexp', '^(?!(auxdata)$).');
    end
end

% exportToZip(mfilename,[adata_dir 'Summaries\' mouse{1} '_site' num2str(exp(1)) 'loadtseries_script_backup']);



end
