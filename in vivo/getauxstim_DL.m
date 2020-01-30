function [auxdata ids ids_new frametimes stimarray] = getauxstim_DL(FiID,diri, eyes,level, bscope2)

% bscope2 = 1 ;

if bscope2
     disp('USING NEW AUXDATA CHANNEL CONVENTIONS FOR BSCOPE2!');
end

if nargin<2
    diri=cd;
end

if nargin<3
    eyes=1;
end

cd(diri)

movie_c = 0;
ret_c = 0;
dat = diri(end-9:end);

runextract = 1;

PC = getenv('computername');

if bscope2;
    aux_ch_frame_times = 4;
    aux_ch_running =  15;
else
    aux_ch_frame_times = 3;
end


auxdir=['..\..\Data\' dat '\'];
stimdir=['..\..\_stim\' dat '\'];

auxfile = dir([auxdir '\*'  num2str(FiID) '*.lvd*']);
stimfile = dir([stimdir '\*'  num2str(FiID) '*.txt*']);

auxdata = load_lvd([auxdir auxfile.name]); ids = [];
frametimes = get_frame_times(auxdata(aux_ch_frame_times,:));



if isempty(stimfile)%at some point fix the error in getauxstim.m which cuases it to fail when the stimscript is made the first time
    stimfilerec = dir([stimdir '\*'  num2str(FiID) '*.mat']);
    if isempty(stimfilerec)
        stimfilerec = dir([stimdir '\*'  num2str(FiID) '*']);
        
        if length(stimfilerec) == 2
            for oo = 1:length(stimfilerec)
                [as bs cs ] = fileparts(stimfilerec(oo).name);
                if isempty(cs)
                    founder = oo
                end
            end
        else
            founder = 1;
        end
    else
        founder = 1;
        
    end
    
    try
        stimarray = load([stimdir stimfilerec(founder).name], '-mat');
        
        if ~isequal(stimarray.angles(stimarray.save_order),stimarray.save_angles)
        disp('- - - - - - - - - - - - - - - - - - ')
        disp('Workaraound for Joel''s angle orders. Should only be running when stimarray.save_angles is not in the correct randomized order')
        disp('- - - - - - - - - - - - - - - - - - ')        
        stimarray.save_angles =  stimarray.angles(stimarray.save_order);        
        end
        
        if isfield(stimarray, 'movie_path')
            disp('Movie Stim. NO IDS')
            movie_c = 1;
        elseif isfield(stimarray, 'StimRecorder')
            disp('Retinotopy Stim. NO IDS')
            ret_c = 1;
            
        else
            try
                stimfiler([stimdir stimfilerec.name(1:end-4) '.txt'], stimarray.orientations, ...
                    stimarray.start_orientation, stimarray.sf, stimarray.tf, stimarray.repetitions, ...
                    stimarray.animal_dist_cm, stimarray.stimtype, stimarray.sinwave, stimarray.dateStr, ...
                    stimarray.movieDurationSecs, stimarray.inter_stim_interval, length(frametimes), [], ...
                    reshape(stimarray.save_angles',1,size(stimarray.save_angles,2)*size(stimarray.save_angles,1)), ...
                    reshape(stimarray.save_eyes',1,size(stimarray.save_eyes,2)*size(stimarray.save_eyes,1)), [], ...
                    [], [], [], FiID, [] , [], [], []);
            catch
                stimfiler([stimdir stimfilerec.name(1:end-4) '.txt'], stimarray.orientations, ...
                    stimarray.start_orientation, stimarray.sf, stimarray.tf, stimarray.repetitions, ...
                    stimarray.animal_dist_cm, stimarray.stimtype, stimarray.sinwave, stimarray.dateStr, ...
                    stimarray.movieDurationSecs, stimarray.inter_stim_interval, length(frametimes), [], ...
                    reshape(stimarray.save_angles',1,size(stimarray.save_angles,2)*size(stimarray.save_angles,1)), ...
                    [], [], ...
                    [], [], [], FiID, [] , [], [], []);
            end
            
            
        end
    end
     stimfile = [];
else
    stimfilerec = dir([stimdir '\*'  num2str(FiID) '*.mat*']);
    if isempty(stimfilerec)
        stimfilerec = dir([stimdir '\*'  num2str(FiID) '*']);
        
        if length(stimfilerec) == 2
            for oo = 1:length(stimfilerec)
                [as bs cs ] = fileparts(stimfilerec(oo).name);
                if isempty(cs)
                    founder = oo;
                end
            end
        else
            founder = 1;
        end
             
        % The data aquasition toolbox is removed from the
        % path before loading this mat file to avoid the
        % uncecessary activation of the floating license of
        % said toolbox.
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daq')
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqguis')
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqdemos')
        stimarray = load([stimdir stimfilerec(founder).name], '-mat');
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daq')
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqguis')
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqdemos')
    else
        founder = 1;
        
        % The data aquasition toolbox is removed from the
        % path before loading this mat file to avoid the
        % uncecessary activation of the floating license of
        % said toolbox.
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daq')
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqguis')
        rmpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqdemos')
        stimarray = load([stimdir stimfilerec(founder).name], '-mat');
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daq')
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqguis')
        addpath('C:\MATLAB\R2013b\64bit\toolbox\daq\daqdemos')
        
    end
end

try 
    if ~isequal(stimarray.angles(stimarray.save_order),stimarray.save_angles)
        if isequal(FiID,'60335') || isequal(FiID,'60337') || isequal(FiID,'60338')
            disp('- - - - - - - - - - - - - - - - - - ')
            disp('Workaraound for Joel''s angle orders. File ID is 60335 or 60337 or 60338.')
            disp('- - - - - - - - - - - - - - - - - - ')
            stimarray.save_angles =  stimarray.angles(stimarray.save_order);
        end
    end
end


try
    if movie_c
        ids =  getStimIDs_movie_TR(auxdata, 'eyes', eyes, 'level', level, 'bscope', bscope2);
        ids_new = ids;
    elseif ret_c
        ids =  GET_StimFrames_pieter(auxdata, stimarray.StimSettings, 4);
        ids_new = ids;
    else
        if ~isfield(stimarray, 'stimshutter')
            stimshutter = 0;
        else
            stimshutter = stimarray.stimshutter;
        end
       
        [ids ids_new] = getStimIDs_TR_JB(auxdata, 'stimfile', [stimdir stimfile.name], 'eyes', eyes, 'level', level, 'bscope', bscope2, 'stimshutter', stimshutter);
    end
catch
    disp('Crashed while trying to get StimIDs')
    disp('Maybe the auxfile is too short. Adding a single trial with standard parameters')
    try
        try
            [ids ids_new] =  getStimIDs_TR_JB(auxdata, 'stimfile', [stimdir stimfile.name], 'eyes', eyes, 'level', level, 'special', 1, 'bscope', bscope2);
        catch
            [ids ids_new] =  getStimIDs_DL(auxdata, 'stimfile', [stimdir stimfile.name], 'eyes', eyes, 'level', level, 'special', 1, 'bscope', bscope2);
%         catch
%             [ids ids_new] =  getStimIDs_TR_JJ(auxdata, 'stimfile', [stimdir stimfile.name], 'eyes', eyes, 'level', level, 'special', 1, 'bscope', bscope2);
        end
    catch
        disp('That FAILED');
        disp('getStimIDs was unable to aquire stimuslus ids. Continuing without varialbe ids or ids_new.');
        ids = [];
        ids_new = [];
    end
    
    ids = [];
end

if isempty(ids)
    try
        if movie_c
            ids =  getStimIDs_movie_TR(auxdata, 'eyes', eyes, 'level', level);
        else
            
            beep(); pause(0.3)
            beep();
            pause(0.3);
            beep(); pause(0.3)
            beep();
            
            %             disp('ATTENTION! I HAD TO MODIFY THE AUXFILE TO RUN THROUGH!')
            %             ids =  getStimIDs_TR(auxdata,'eyes', eyes, 'level', level, 'stimfile', [stimdir stimfile.name], 'special', 1);
        end
    catch
        ids = [];
    end
end


% view_tiff_stack_plot_aux(data, [auxdata(7,frametimes)' auxdata(8,frametimes)']);



