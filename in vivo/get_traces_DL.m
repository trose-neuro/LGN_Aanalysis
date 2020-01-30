function [ODI, PO, PD, OS] = get_traces_DL(date, mouse, exp, bscope2,ROIidx, ROIOI, showplots, adata_dir, rdata_dir)
%This function will get the df/f of each stimulus for a single ROI of a
%single experiment.
%If plots are skipped then it is about 7 times faster.
% d
% eg.
% mouse = 'JB_160218';
% exp = 13747;
% date= '2016-04-22';
% bscope2 = 1;
% ROIidx = 1;
% showplot = 1;



%% load ROI and stimulus data
load([adata_dir '\' mouse '\' date '\' mouse '-Adata-' num2str(exp) '.mat']);
cd([rdata_dir '\' mouse '\ImagingData\' date]);

%% extract basal parameters
SamplingFreq1 = regexp(info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
SamplingFreq1 = str2num(SamplingFreq1{1}) / info.level;

%% parameter switchboard
load_group = 1; % dendrite/spine group (e.g. 19, 23, 22)
F0per = 0.9; % duration of prestimulus period taken as F0

%% concatenation of all presstimulus period indices (both eyes)
e1_base = [ids_new(1).stim_boundaries(:,:,2)-1-round(F0per*SamplingFreq1)];%:ids_new(1).stim_boundaries(:,:,2)-1;
e1_base = sort([e1_base(:)]);
e1_baseend = ids_new(1).stim_boundaries(:,:,2)-1;
e1_baseend = sort([e1_baseend(:)]);

e2_base = [ids_new(2).stim_boundaries(:,:,2)-1-round(F0per*SamplingFreq1)];%:ids_new(1).stim_boundaries(:,:,2)-1;
e2_base = sort([e2_base(:)]);
e2_baseend = ids_new(2).stim_boundaries(:,:,2)-1;
e2_baseend = sort([e2_baseend(:)]);

base_idx(:,1) = sort([e1_base ; e2_base]);
base_idx(:,2) = sort([e1_baseend ; e2_baseend]);

%% concatenation of stimulus period indices (both eyes)
e1_stim = [ids_new(1).stim_boundaries(:,:,2)]; %:ids_new(1).stim_boundaries(:,:,2)-1;
e1_stim = sort([e1_stim(:)]);
e1_stimend = ids_new(1).stim_boundaries(:,:,3);
e1_stimend = sort([e1_stimend(:)]);

e2_stim = [ids_new(2).stim_boundaries(:,:,2)]; %:ids_new(1).stim_boundaries(:,:,2)-1;
e2_stim = sort([e2_stim(:)]);
e2_stimend = ids_new(2).stim_boundaries(:,:,3);
e2_stimend = sort([e2_stimend(:)]);

stim_idx(:,1) = sort([e1_stim ; e2_stim]);
stim_idx(:,2) = sort([e1_stimend ; e2_stimend]);

%% extraction and ratioing of grenn/red spine and dendrite signal
if isempty(ROIidx)
    cell_OI = ROIs(ROIOI).activity; %cell_OI = cell of interest
    %r_cell = ROIs(1).activity_r;
else
    cell_OI = ROIs(ROIidx).activity;
    %r_cell = ROIs(ROIidx).activity_r;
end
if mean(cell_OI(:))<100 % If the mean red flourescence is over 100 ratiometric analysis is done (basically if the cell has double construct). 
    noratio = 0;
    disp('not doing ratiometric trace extraction')
    cell_OI_nR = cell_OI;   %cell_OI_nR = cell of interest, non-ratiometric
else
    %gr_cell = g_cell./ r_cell;
    gr_cell = cell_OI;
end



%% extraction and subtraction of baseline ratios


R0_cell_trace = NaN(size(cell_OI_nR));
R0_cell = [];
reclength = length(cell_OI);
timebase = linspace(0, reclength/SamplingFreq1, reclength); %from 0 to reclength in seconds (reclength[frames]/SamplingFreq1[Hz])

for i = 1:length(base_idx)
    R0_cell_trace(base_idx(i,1):base_idx(i,2)) = cell_OI_nR(base_idx(i,1):base_idx(i,2)); %startpoints: base_idx(:,1); endpoints: base_idx(:,2)
    R0_cell(i) = nanmean(cell_OI_nR(base_idx(i,1):base_idx(i,2))); %calculate the baseline value for cell_OI_nR
    x(i) = mean(base_idx(i,1):base_idx(i,2)); %x-axis: pre-stimulus frames
end

R0_cell_interp = interp1(x, smooth(R0_cell,12), 1:reclength ,'linear', 'extrap'); %interpolate all 3 parameters into each other/one graphline
dFF_cell_OI_nR = (cell_OI_nR - R0_cell_interp) ./ R0_cell_interp * 100;
% dFF_gr_cell = (gr_cell - R0_cell) ./ R0_cell * 100;


%% blanck inter trial intervals
% dFF_gr_cell_nan = NaN(1,reclength);
% for i = 1:length(stim_idx)
%     dFF_gr_cell_nan(stim_idx(i,1):stim_idx(i,2)) = dFF_gr_cell(stim_idx(i,1):stim_idx(i,2));
% end

if showplots==1
    %% Raw trace plot
    figure
    plot(timebase, dFF_cell_OI_nR, 'k'); hold all
    axis tight; ylim([min(dFF_cell_OI_nR(500:end-500)) max(dFF_cell_OI_nR(500:end-500))*1.2]);
    title('Uncorrected R/R_0 data');
    xlabel('Time [s]');
    ylabel('\DeltaR/R_0');
    make_nice(stimarray, ids_new, stim_idx, timebase); %%%%% throws error, angles are probably wrong. 
end

%% make a OD/OS plot
tempeyes = stimarray.save_eyes'; % 0 = ipsi; 1 = contra !!!!!!!!!!!!!!
tempeyes = tempeyes(:);
% tempangles = (stimarray.save_angles-90)'/360*(stimarray.orientations / 2)+1; %assume 90° rotation for bscope1
tempangles = (stimarray.save_angles)'/360*(stimarray.orientations / 2)+1; %no rotation assumed
tempangles = tempangles(:);

dFF_cell_OI_nR_indexed(1,:) = dFF_cell_OI_nR;
dFF_cell_OI_nR_indexed(2,:) = timebase;
dFF_cell_OI_nR_indexed(3,:) = NaN(1,length(dFF_cell_OI_nR)); 
dFF_cell_OI_nR_indexed(4,:) = NaN(1,length(dFF_cell_OI_nR));
dFF_cell_OI_nR_indexed(5,:) = NaN(1,length(dFF_cell_OI_nR));
dFF_cell_OI_nR_indexed(6,:) = NaN(1,length(dFF_cell_OI_nR));

for i= 1:length(tempeyes)
    dFF_cell_OI_nR_indexed(3,stim_idx(i,1):length(dFF_cell_OI_nR_indexed)) = tempeyes(i); % 0 = ipsi; 1 = contra;
    
end
dFF_cell_OI_nR_indexed(3,stim_idx(end,1)+round((stimarray.movieDurationSecs+stimarray.inter_stim_interval)*30.1635/4):end) = NaN; % changes all values past the end of the last post tim period to NaN

repcount = repmat([1:stimarray.repetitions]',[1,(stimarray.orientations)])';
repcount = repcount(:);

for i= 1:length(tempangles)
    dFF_cell_OI_nR_indexed(4,stim_idx(i,1):length(dFF_cell_OI_nR_indexed)) = tempangles(i);
    dFF_cell_OI_nR_indexed(5,stim_idx(i,1):length(dFF_cell_OI_nR_indexed)) = repcount(i);
end

% This is a control so that only frames within the stimulation period are used in further analysis.
dFF_cell_OI_nR_indexed(6,:) = ids_new(1).stim_frames;
dFF_cell_OI_nR_indexed(6,dFF_cell_OI_nR_indexed(6,:) > 0.5) = 1;
dFF_cell_OI_nR_indexed(6,dFF_cell_OI_nR_indexed(6,:) <= 0.5) = 0;


for i = 1:length(dFF_cell_OI_nR)
    if isnan(dFF_cell_OI_nR_indexed(3,i))
        dFF_cell_OI_nR_indexed(1,i) = NaN;
    end
end

[~,I]=sort(dFF_cell_OI_nR_indexed(4,:));
dFF_cell_OI_nR_indexed=dFF_cell_OI_nR_indexed(:,I);
[~,I]=sort(dFF_cell_OI_nR_indexed(3,:));
dFF_cell_OI_nR_indexed=dFF_cell_OI_nR_indexed(:,I);



%ylimits(2) = round(max(dFF_gr_cell_indexed(1,dFF_gr_cell_indexed(3,:)==or(0,1))));
% ylimits(2) = max(dFF_gr_cell_indexed(1,500:end-100));
% ylimits(1) = round(min(dFF_gr_cell_indexed(1,dFF_gr_cell_indexed(3,:)==or(0,1))));

ylimits(2)= max(dFF_cell_OI_nR_indexed(1,stim_idx(1,1):stim_idx(end,2)));
ylimits(1)= min(dFF_cell_OI_nR_indexed(1,stim_idx(1,1):stim_idx(end,2)));

if showplots==1
    if stimarray.orientations / 2 == 8
        x_sym_strings = 'VXJLNPRT';
    elseif stimarray.orientations / 2 == 12
        x_sym_strings = 'VWYJKMNOQRSU'; 
    end
    figure('color', 'w');
end

for i = 1:stimarray.orientations/2
    clear peak_dFF
    temptrace{1} = dFF_cell_OI_nR_indexed([1,5,6],(dFF_cell_OI_nR_indexed(4,:)==i & dFF_cell_OI_nR_indexed(3,:)==0)); %  temptrace(1) = ipsi;  temptrace(2) = contra;
    temptrace{2} = dFF_cell_OI_nR_indexed([1,5,6],(dFF_cell_OI_nR_indexed(4,:)==i & dFF_cell_OI_nR_indexed(3,:)==1)); 
    temptrace_mean{1} = nan(stimarray.repetitions, length(temptrace{1}(1,(temptrace{1}(2,:)==1)))+10);
    temptrace_mean{2} = nan(stimarray.repetitions, length(temptrace{2}(1,(temptrace{2}(2,:)==1)))+10);
    
    for eye = 1:2  % 1 is ipsi, 2 is contra
        eye_colour = 'rb'; % ipsi is red , contra is blue
        if showplots==1
            if eye ==1
                subplot(2,stimarray.orientations/2,i); axis off;
            elseif eye == 2
                subplot(2,stimarray.orientations/2,i+stimarray.orientations/2); axis off;
            end
            hold on; 
            patch([0 (stim_idx(1,2)-stim_idx(1,1)) (stim_idx(1,2)-stim_idx(1,1)) 0], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)],eye_colour(eye),'facealpha', 0.5,'edgecolor','none', 'Tag', 'Stim');
        end
        for k = 1:stimarray.repetitions
            clear single_trace
            single_trace(1,:) = temptrace{eye}(1,(temptrace{eye}(2,:)==k));
            single_trace(2,:) = temptrace{eye}(3,(temptrace{eye}(2,:)==k));
            if showplots==1
                plot(single_trace(1,:),'k','linewidth',1.5);
            end
            temptrace_mean{eye}(k,1:length(single_trace(1,:)))=single_trace(1,:);
            single_trace = single_trace(1, single_trace(2,:)==1);
            peak_dFF(eye,k) = max(single_trace(2:round(stimarray.movieDurationSecs*30.1635/4-5))); %not sure this is the best way of determining the end of the stim...
        end

        mean_peak_dFF(eye, i) = mean(peak_dFF(eye,:));
        temptrace_mean{eye} = nanmean(temptrace_mean{eye});

        if showplots==1
            plot(temptrace_mean{eye},'r','linewidth',3); ylim(ylimits);set(gca,'Color','w');
            title(x_sym_strings(i),'FontName', 'OriSymbols', 'FontSize',50,'LineStyle', 'none', 'Color', eye_colour(eye)); axis off
            lhdl = hline(0, '-k'); set(lhdl, 'LineStyle', ':'); 
        end
    end
    
end


%% Calculate and indicators
if stimarray.orientations/2 == 12
    true_angles = [30:30:360];
end

OS(1) = 1-TT_CircularVariance(circshift(mean_peak_dFF(1,:)',-4)')% ipsi orientation selectivity
OS(2) = 1-TT_CircularVariance(circshift(mean_peak_dFF(2,:)',-4)') % contra orientation selectivity
PD(1) = TT_PreferredDirection(circshift(mean_peak_dFF(1,:)',-4)') % ipsi preferred_direction
PD(2) = TT_PreferredDirection(circshift(mean_peak_dFF(2,:)',-4)') % contra preferred_direction
PO(1) = TT_PreferredOrientation(circshift(mean_peak_dFF(1,:)',-4)') % ipsi preferred orientation
PO(2) = TT_PreferredOrientation(circshift(mean_peak_dFF(2,:)',-4)') % contra preferred orientation

quick_pref_ori(1) = find(mean_peak_dFF(1,:) == max(mean_peak_dFF(1,:))); % ipsi
quick_pref_ori(2) = find(mean_peak_dFF(2,:) == max(mean_peak_dFF(2,:))); % contra

% (contra-ipsi)/(contra+ipsi)
ODI = (mean_peak_dFF(2,quick_pref_ori(2))-mean_peak_dFF(1,quick_pref_ori(1)))/(mean_peak_dFF(2,quick_pref_ori(2))+mean_peak_dFF(1,quick_pref_ori(1)))

if showplots==1 
    figure;
    plot(true_angles,circshift(mean_peak_dFF(1,:)',-4)','r'); hold on;
    plot(true_angles,circshift(mean_peak_dFF(2,:)',-4)','b');
    set(gca,'xtick',true_angles)
    xlabel('Direction of moving gratings (°)');
    ylabel('Mean response amplitude ([dR/R])');
    title(['Tuning Curve (ODI=', num2str(ODI), ')']);
    legend(['Ipsi: OS=' num2str(OS(1)) ', PD='  num2str(PD(1)) ', PO='  num2str(PO(1))],['Contra: OS='  num2str(OS(2)) ', PD='  num2str(PD(2)) ', PO='  num2str(PO(2))])
end







%ORI

%1-circular var


%ODI

% plotOri(gr_cell,ids_new(1),'plottype', 'tuningpolar'); %contra
% plotOri(gr_cell,ids_new(2),'plottype', 'tuningpolar'); %ipsi

end

function linkaxesInFigure(varargin)
% linkaxesInFigure - Finds all visible axes in figure and links them for zooming
% 
% Syntax: linkaxesInFigure(varargin)
% varargin - Can be 'x', 'y', 'xy', 'off'  Same functionality as linkaxes
%   Example:
%       figure;
%       subplot(2,1,1);
%       plot(rand(10,1));
%       subplot(2,1,2);
%       plot(1:10);
%       linkaxesInFigure('x')
%
%   See also: linkaxes

% AUTHOR    : Dan Kominsky
% Copyright 2012  Prime Photonics, LC.
%%
  if (nargin == 0) || ~ischar(varargin{1})
    linkAx = 'xy';
  else
    linkAx = lower(varargin{1});
  end
  
  x = findobj(gcf,'Type','axes','Visible','on');
  try
    linkaxes(x,linkAx)
  catch ME
    disp(ME.message)
  end



end 
function make_nice(stimarray, ids, stim_idx, timebase)
%% ORI symols and eye patches
stims = stimarray.orientations / 2;

if stims == 8
    x_sym_strings = 'VXJLNPRT';
    
    x_col = 'ygbrygbr';
    %     x_col = 'kkkkkkkk';
    
elseif stims == 12
    x_sym_strings = 'VWYJKMNOQRSU';
    x_col = 'kkkkkkkkkkkkk';
    x_sym_strings = [x_sym_strings x_sym_strings];
    x_col = [x_col x_col];
elseif stims == 16
    x_sym_strings = 'VWXYJKLMNOPQRSTU';
    x_col = 'kkkkkkkkkkkkkkkkk';
    x_sym_strings = [x_sym_strings x_sym_strings];
    x_col = [x_col x_col];
end

tempeyes = stimarray.save_eyes'; % 0 = ipsi; 1 = contra;  !!!!!!!!!!!
tempeyes = tempeyes(:);

tempangles = (stimarray.save_angles)'/360*stims+1; %assume 90° rotation
tempangles = tempangles(:);

% % Alternative to using stimarray when stimarray if fucked up. 
% eyeseq=stimarray.save_eyes(:)+1;
% k=[1,1]; %counter
% for i = 1:length(eyeseq)
%     tempangles(i) = ids(eyeseq(i)).stimseq_deg(k(eyeseq(i)));
%     k(eyeseq(i))= k(eyeseq(i))+1;
% end
% tempangles = (tempangles)'/360*stims+1;






yl = ylim;
for i = 1:length(stim_idx)
    x = [timebase(stim_idx(i,1)) timebase(stim_idx(i,2))];
    y = [yl(2) yl(2)];
    y2 = [yl(1) yl(1)];
    if tempeyes(i) == 1 ; %contra
        patch([x fliplr(x)],[y fliplr(y2)], 'b','facealpha',0.1,...
            'edgecolor','none', 'Tag', 'Stim');
    elseif tempeyes(i) == 0 ; % ipsi
        patch([x fliplr(x)],[y fliplr(y2)], 'r','facealpha',0.1,...
            'edgecolor','none', 'Tag', 'Stim');
    elseif tempeyes(i) == 2;
        patch([x fliplr(x)],[y fliplr(y2)], 'w','facealpha',0.1,...
            'edgecolor','none', 'Tag', 'Stim');
    end
    %%%% error as tempangles(1) is -1 anc cannot be used as reference for x_sym_strings 
    text(x(1) + diff(x)/3, y(1)-0.1*y(1), x_sym_strings(tempangles(i)), 'FontName', 'OriSymbols', 'FontSize',14,'LineStyle', 'none', 'Color', x_col(tempangles(i)));
end

%% figure format defaults

% Defaults for Cell press. 1 col: 85mm, 1.5 col: 114mm, 2col:174mm
% Defaults for Nature press. 1 col: 89mm, 1.5 col: 136mm, 2col:183mm
% width = 8.9;                  % Width in cm
% height = GR*width;     % Height in cm (golden ratio default  (1 + n.sqrt(5)) / 2 1/GR)
alw = 1;                 % AxesLineWidth
fsz = 18;      % Fontsize
lw = 3;      % LineWidth
opengl software % to get the axes on figures with transparency back!
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fsz)
% set(gcf, 'TickLength', [0.02 0.025]);
set(get(gca, 'Ylabel'), 'FontSize', fsz)
set(get(gca, 'Xlabel'), 'FontSize', fsz)
set(get(gca, 'Title'), 'FontSize', fsz)
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gca, 'visible', 'off');
set(gca, 'Color', 'none');
lhdl = hline(0, '-k'); set(lhdl, 'LineStyle', ':')
end


