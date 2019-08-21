
%% load data
reload = 0;

if reload
    load('/Users/trose/Documents/GitHub/LGN_Analysis/Full_Data.mat')
    load('/Users/trose/Documents/GitHub/LGN_Analysis/in vivo/base_MD_ODI_array.mat')
    load('/Users/trose/Documents/GitHub/LGN_Analysis/in vivo/ana_ODI_Musc.mat')
end

clearvars -except ana_ODI ana_contra ana_ipsi BaseMD_ODI_array data
close all

savedir = '/Users/trose/Documents/GitHub/LGN_Analysis/figures';
savefig = 1;

%% figure settings

% Defaults for Cell press. 1 col: 85mm, 1.5 col: 114mm, 2col:174mm
% Defaults for Nature press. 1 col: 89mm, 1.5 col: 136mm, 2col:183mm

fig.width= 4;              % fig.widthin cm
fig.height = fig.width;         % fig.height in cm (golden ratio default  (1 + n.sqrt(5)) / 2 1/GR)
fig.alw = 1;                % AxesLineWidth
fig.fsz = 18;               % Fontsize
fig.ln = 2;                 % lineWidth for e.g. averages
colorlevels = 3;
fig.coc_c = cbrewer('seq', 'Blues', colorlevels);
fig.coc_i = cbrewer('seq', 'Reds', colorlevels);
coc_a = [0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];
fig.coc_a = flipud(coc_a);
fig.markercol = [0.85 0.85 0.85];
fig.markersz = 15;

fig.basecol = [0 0 0];
fig.MDcol   = [.5 .5 .5];

% opengl software         % to get the axes on figures with transparency back!
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultAxesFontSize', fig.fsz)
set(0,'DefaultAxesTickLength', [0.02 0.025]);

%% histogram settings
bins = 7;
binsOD  = -1:2/bins:1; % bin edges
binsODcentres = conv(binsOD, [0.5 0.5], 'valid');

%% FIND CONTROL AND MD idx
field_number = size(data,2);
md=find([data(:).MD]==1);
cont=find([data(:).MD]==0);

%% Data structure read out for Control and MD
%Control
control_peakramp_AMPA_traces = zeros(length(cont),1000);

for i=1:length(cont)
    %Category read out
    category(i)     = data(cont(i)).ocular_category;
    ODI_AMPA_r(i)   = data(cont(i)).ODI_AMPA_step_peak;
    ODI_NMDA_r(i)   = data(cont(i)).ODI_NMDA_step_peak;
    data_cont(i)    = data(cont(i));
    if strcmp(data(cont(i)).experimentator,'SW')
        srF = 1;
    else
        srF = 2;
    end
    if data(cont(i)).brain_contra_ipsi %normal case. Contra Chrimson/Red
        if length(data(cont(i)).step_red.steps_use_AMPA) == 1
            Contra_AMPA_Amp(i)                  = abs(data(cont(i)).step_red.neg_peak1(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_red.neg_base_mean1(2,data(cont(i)).step_red.steps_use_AMPA) );
            Ipsi_AMPA_Amp(i)                    = abs(data(cont(i)).step_blue.neg_peak2(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_blue.neg_base_mean2(2,data(cont(i)).step_red.steps_use_AMPA) );
            
            control_peakramp_AMPA_traces(i,:)   = data(cont(i)).step_red.ephys_traces_70(1:srF:end,data(cont(i)).step_red.steps_use_AMPA,2)';
        else
            Contra_AMPA_Amp(i)                  = 0;
            Ipsi_AMPA_Amp(i)                    = nanmean(abs(data(cont(i)).step_blue.neg_peak2(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_blue.neg_base_mean2(2,data(cont(i)).step_red.steps_use_AMPA) ) );
            control_peakramp_AMPA_traces(i,:)   = data(cont(i)).step_red.ephys_traces_70(1:srF:end,data(cont(i)).step_red.steps_use_AMPA(1),2)';
        end
    else                              %Contra Chronos / Blue
        if length(data(cont(i)).step_red.steps_use_AMPA) == 1
            Ipsi_AMPA_Amp(i)                  = abs(data(cont(i)).step_red.neg_peak1(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_red.neg_base_mean1(2,data(cont(i)).step_red.steps_use_AMPA) );
            Contra_AMPA_Amp(i)                    = abs(data(cont(i)).step_blue.neg_peak2(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_blue.neg_base_mean2(2,data(cont(i)).step_red.steps_use_AMPA) );
            
            %             control_peakramp_AMPA_traces(i,:)   = data(cont(i)).step_red.ephys_traces_70(1:srF:end,data(cont(i)).step_red.steps_use_AMPA,2)';
        else
            Ipsi_AMPA_Amp(i)                  = 0;
            Contra_AMPA_Amp(i)                    = nanmean(abs(data(cont(i)).step_blue.neg_peak2(2,data(cont(i)).step_red.steps_use_AMPA) ...
                - data(cont(i)).step_blue.neg_base_mean2(2,data(cont(i)).step_red.steps_use_AMPA) ) );
            %             control_peakramp_AMPA_traces(i,:)   = data(cont(i)).step_red.ephys_traces_70(1:srF:end,data(cont(i)).step_red.steps_use_AMPA(1),2)';
            %             control_peakramp_AMPA_traces(i,:)   = data(cont(i)).step_red.ephys_traces_70(:,data(cont(i)).step_red.steps_use_AMPA(1),2)';
        end
    end
    ODI_AMPA_cont_sanity(i)                  = ( Contra_AMPA_Amp(i) - Ipsi_AMPA_Amp(i) ) /  ( Contra_AMPA_Amp(i) + Ipsi_AMPA_Amp(i) );
end

ODI_AMPA_cont_sanity(ODI_AMPA_cont_sanity>1)  = 1;
ODI_AMPA_cont_sanity(ODI_AMPA_cont_sanity<-1) = -1;

[prsort idxc] = sort(ODI_AMPA_cont_sanity);

% view_tiff(control_peakramp_AMPA_traces(idxc,:))

%MD
for i=1:length(md)
    %Category read out
    category_md(i)=data(md(i)).ocular_category;
    ODI_AMPA_r_md(i)=data(md(i)).ODI_AMPA_step_peak;
    ODI_NMDA_r_md(i)=data(md(i)).ODI_NMDA_step_peak;
    data_md(i)=data(md(i));
    
    if strcmp(data(md(i)).experimentator,'SW')
        srF = 1;
    else
        srF = 2;
    end
    if data(md(i)).brain_contra_ipsi %normal case. Contra Chrimson/Red
        if length(data(md(i)).step_red.steps_use_AMPA) == 1 & ~isnan((data(md(i)).step_red.steps_use_AMPA))
            Contra_AMPA_Amp_MD(i)                  = abs(data(md(i)).step_red.neg_peak1(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_red.neg_base_mean1(2,data(md(i)).step_red.steps_use_AMPA) );
            Ipsi_AMPA_Amp_MD(i)                    = abs(data(md(i)).step_blue.neg_peak2(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_blue.neg_base_mean2(2,data(md(i)).step_red.steps_use_AMPA) );
            
            MD_peakramp_AMPA_traces(i,:)   = data(md(i)).step_red.ephys_traces_70(1:srF:end,data(md(i)).step_red.steps_use_AMPA,2)';
        elseif length(data(md(i)).step_red.steps_use_AMPA) > 1 & ~isnan((data(md(i)).step_red.steps_use_AMPA))
            Contra_AMPA_Amp_MD(i)                  = 0;
            Ipsi_AMPA_Amp_MD(i)                    = nanmean(abs(data(md(i)).step_blue.neg_peak2(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_blue.neg_base_mean2(2,data(md(i)).step_red.steps_use_AMPA) ) );
            MD_peakramp_AMPA_traces(i,:)   = data(md(i)).step_red.ephys_traces_70(1:srF:end,data(md(i)).step_red.steps_use_AMPA(1),2)';
        end
    else                              %Contra Chronos / Blue
        if length(data(md(i)).step_red.steps_use_AMPA) == 1 & ~isnan((data(md(i)).step_red.steps_use_AMPA))
            Ipsi_AMPA_Amp_MD(i)                  = abs(data(md(i)).step_red.neg_peak1(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_red.neg_base_mean1(2,data(md(i)).step_red.steps_use_AMPA) );
            Contra_AMPA_Amp_MD(i)                    = abs(data(md(i)).step_blue.neg_peak2(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_blue.neg_base_mean2(2,data(md(i)).step_red.steps_use_AMPA) );
            
            %             control_peakramp_AMPA_traces(i,:)   = data(md(i)).step_red.ephys_traces_70(1:srF:end,data(md(i)).step_red.steps_use_AMPA,2)';
        elseif length(data(md(i)).step_red.steps_use_AMPA) > 1 & ~isnan((data(md(i)).step_red.steps_use_AMPA))
            Ipsi_AMPA_Amp_MD(i)                  = 0;
            Contra_AMPA_Amp_MD(i)                    = nanmean(abs(data(md(i)).step_blue.neg_peak2(2,data(md(i)).step_red.steps_use_AMPA) ...
                - data(md(i)).step_blue.neg_base_mean2(2,data(md(i)).step_red.steps_use_AMPA) ) );
            %             control_peakramp_AMPA_traces(i,:)   = data(md(i)).step_red.ephys_traces_70(1:srF:end,data(md(i)).step_red.steps_use_AMPA(1),2)';
            %             control_peakramp_AMPA_traces(i,:)   = data(md(i)).step_red.ephys_traces_70(:,data(md(i)).step_red.steps_use_AMPA(1),2)';
        end
    end
    if ~isnan((data(md(i)).step_red.steps_use_AMPA))
        ODI_AMPA_MD_sanity(i)                  = ( Contra_AMPA_Amp_MD(i) - Ipsi_AMPA_Amp_MD(i) ) /  ( Contra_AMPA_Amp_MD(i) + Ipsi_AMPA_Amp_MD(i) );
    end
end


ODI_AMPA_MD_sanity(ODI_AMPA_MD_sanity>1)  = 1;
ODI_AMPA_MD_sanity(ODI_AMPA_MD_sanity<-1) = -1;

[prsort idxm] = sort(ODI_AMPA_MD_sanity);


%% CONTROL
%AMPA fraction discretize
temp=ODI_AMPA_r;
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    AMPA_ODI_baseline_fract(i)=(length(find(disc==i)))/length(temp);
end
disc=[];

%NMDA fraction discretize
temp2=ODI_NMDA_r;
temp2(find(isnan(temp2)))=[];
disc=discretize(temp2,[binsOD]);
for i=1:bins
    NMDA_ODI_baseline_fract(i)=(length(find(disc==i)))/length(temp2);
end
disc=[];

%% MD
%AMPA fraction discretize
temp=ODI_AMPA_r_md;
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    AMPA_ODI_MD_fract(i)=(length(find(disc==i)))/length(temp);
end
disc=[];
%NMDA fraction discretize
temp2=ODI_NMDA_r_md;
temp2(find(isnan(temp2)))=[];
disc=discretize(temp2,[binsOD]);
for i=1:bins
    NMDA_ODI_MD_fract(i)=(length(find(disc==i)))/length(temp2);
end
disc=[];

% figure
n = 1;
figure(n);
fig.ImageDescription = 'AMPA ODI baseline';

cl = bar(binsODcentres,AMPA_ODI_MD_fract,'barwidth', 1, 'FaceColor', 'w','EdgeColor', 'k', 'LineWidth', fig.alw);
set(cl, 'FaceColor', 'w', 'EdgeColor', 'k',  'LineWidth', fig.alw);
xlim([-1 1]);
set_fig_properties(n, fig);
xlabel('ODI AMPA'); ylabel('fraction');

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end

% figure

n = n + 1;
figure(n);
fig.ImageDescription = 'NMDA ODI baseline';

cl = bar(binsODcentres,NMDA_ODI_MD_fract,'barwidth', 1, 'FaceColor', 'w','EdgeColor', 'k', 'LineWidth', fig.alw);
set(cl, 'FaceColor', 'w', 'EdgeColor', 'k',  'LineWidth', fig.alw);
xlim([-1 1]);
set_fig_properties(n, fig);
xlabel('ODI NMDA'); ylabel('fraction');

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end


%% compare to in vivo ODIs

% in vivo discretize
% baseline

temp=BaseMD_ODI_array(:,1);
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    Invivo_base_fract(i)=(length(find(disc==i)))/length(temp);
end
disc=[];

% MD
temp=BaseMD_ODI_array(:,2);
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    Invivo_MD_fract(i)=(length(find(disc==i)))/length(temp);
end
disc=[];

% figure
n = n + 1;
figure(n);
fig.ImageDescription = 'In Vivo ODI baseline';

cl = bar(binsODcentres,Invivo_base_fract,'barwidth', 1, 'FaceColor', 'w','EdgeColor', 'k', 'LineWidth', fig.alw);

xlim([-1 1]);
xlabel('ODI'); ylabel('fraction');
set_fig_properties(n, fig);

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end


%% MIs in vivo in vitro
% figure
n = n + 1;
figure(n);
fig.ImageDescription = 'MI comparison';

cl = bar([nanmean(abs(BaseMD_ODI_array(:,1)))  nanmean(abs(ODI_AMPA_r))],'barwidth', .8, 'FaceColor', 'w','EdgeColor', 'k', 'LineWidth', fig.alw); hold on;
clp = plotSpread({abs(BaseMD_ODI_array(:,1)), abs(ODI_AMPA_r)}, 'ShowMM', 4);

set(clp{1}, 'MarkerFaceColor', fig.markercol, 'MarkerEdgeColor', fig.markercol, 'MarkerSize', fig.markersz);
set(clp{2}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'LineWidth', fig.alw);

set(gca,'xticklabel',{'in vivo', 'in vitro'});

xlabel('MI'); ylabel('Monocularity (|ODI|)');

set_fig_properties(n, fig);

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end



%% muscimol ODI
% figure
n = n + 1;
figure(n);
fig.ImageDescription = 'MI comparison AMPA / NMDA';

cl = plot([ones(size(ODI_AMPA_r))' ones(size(ODI_AMPA_r))'*2]', [abs(ODI_AMPA_r)' abs(ODI_NMDA_r)']','-o');
clp = plotSpread({abs(ODI_AMPA_r), abs(ODI_NMDA_r)}, 'ShowMM', 4);

delete(clp{1});

set(cl, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', fig.markersz / 2, 'LineWidth', fig.alw, 'Color', fig.markercol);
set(clp{2}, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r', 'LineWidth', fig.alw * 2);

set(gca,'xticklabel',{'AMPA', 'NMDA'});

xlabel('in vitro'); ylabel('Monocularity (|ODI|)');

set_fig_properties(n, fig);

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end

Figure(n).p = signrank(ODI_AMPA_r, ODI_NMDA_r);
Figure(n).test = 'Wilcoxon signed rank';
Figure(n).ImageDescription = fig.ImageDescription;

%% muscimol ODI
% figure
n = n + 1;
figure(n);
fig.ImageDescription = 'MI comparison Muscimol';

% cl = plot([ones(size(ana_ODI(:,1))) ones(size(ana_ODI(:,1)))*2]', [abs(ana_ODI(:,1)) abs(ana_ODI(:,2))]','-');

cl = bar([nanmean(abs(ana_ODI(:,1))) nanmean(abs(ana_ODI(:,2)))],'barwidth', .8); hold on;
clp = plotSpread({abs(ana_ODI(:,1)), abs(ana_ODI(:,2))}, 'ShowMM', 4);

% delete(clp{1});
% set(cl, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', fig.markersz / 2, 'LineWidth', fig.alw, 'Color', fig.markercol);
set(cl, 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', fig.alw);

set(clp{1}, 'MarkerFaceColor', fig.markercol, 'MarkerEdgeColor', fig.markercol, 'MarkerSize', fig.markersz);
set(clp{2}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'LineWidth', fig.alw);

set(gca,'xticklabel',{'control', 'muscimol'});

xlabel('in vivo'); ylabel('Monocularity (|ODI|)');

set_fig_properties(n, fig);

Figure(n).p = ranksum(abs(ana_ODI(:,1)), abs(ana_ODI(:,2)));
Figure(n).test = 'Mann-Whitney U-test';
Figure(n).ImageDescription = fig.ImageDescription;

%% muscimol AMP
% figure
n = n + 1;
figure(n);
fig.ImageDescription = 'MI comparison Muscimol';

cl = bar([nanmean(ana_contra) nanmean(ana_ipsi)],'barwidth', .8, 'FaceColor', 'w','EdgeColor', 'k', 'LineWidth', fig.alw); hold on;
clp = plotSpread({ana_contra(:,1), ana_contra(:,2), ana_ipsi(:,1), ana_ipsi(:,2)}, 'ShowMM', 4);

set(cl, 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', fig.alw);
set(clp{1}, 'MarkerFaceColor', fig.markercol, 'MarkerEdgeColor', fig.markercol, 'MarkerSize', fig.markersz);
set(clp{2}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Color', 'k', 'LineWidth', fig.alw);

delete(clp{1});

set(gca,'xticklabel',{'contra', 'contra musc', 'ipsi', 'ipsi musc'});

xlabel('in vivo'); ylabel('\DeltaF/F_0');

set_fig_properties(n, fig);

if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end


%% BINO fractions - ALL

control_binofrac = size(find(abs(ODI_AMPA_r)<1),2) / size(ODI_AMPA_r,2);
MD_binofrac = size(find(abs(ODI_AMPA_r_md)<1),2) / size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

n1 = size(find(abs(ODI_AMPA_r)<1),2);
N1 = size(ODI_AMPA_r,2);
n2 = size(find(abs(ODI_AMPA_r_md)<1),2);
N2 = size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

[pval chi2stat tbl] = propstat(n1, N1, n2, N2);

%% BINO fractions - Bino Only

control_binofrac = size(find(abs(ODI_AMPA_r)<1),2) / size(ODI_AMPA_r,2);
MD_binofrac = size(find(abs(ODI_AMPA_r_md)<1),2) / size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

n1 = size(find(abs(ODI_AMPA_r)<1),2);
N1 = size(ODI_AMPA_r,2);
n2 = size(find(abs(ODI_AMPA_r_md)<1),2);
N2 = size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

[pval chi2stat tbl] = propstat(n1, N1, n2, N2);

%% bar graph for categories
mincat = 1;
maxcat = 5;

ctidx = find([data(:).ocular_category] > mincat - 1 & [data(:).ocular_category] < maxcat + 1);
[tbl,chi2,p,labels] = crosstab([data(ctidx).ocular_category]', [data(ctidx).MD]');

%figure
n = n + 1;
figure(n);
fig.ImageDescription = 'Category Comparison MD';

cl = bar(tbl ./ sum(tbl) .* 100,'barwidth', 1.5 ,'EdgeColor', 'k', 'LineWidth', fig.alw); hold on;
% set(gca,'xticklabel',{'contra', 'contra musc', 'ipsi', 'ipsi musc'});
% cl(1).CData = repmat([0 0 0], size(cl(1).CData,1), 1);
% cl(2).CData = repmat([1 1 1], size(cl(2).CData,1), 1);

cl(1).FaceColor = fig.basecol;
cl(2).FaceColor = fig.MDcol;

set_fig_properties(n, fig);
xlabel('Binocularity Categories'); ylabel('Fraction (%)');

Figure(n).p = p;
Figure(n).test = 'Crosstabulated Chi2';
Figure(n).ImageDescription = fig.ImageDescription;


if savefig; fig2svg([savedir filesep 'fig' num2str(n) '.svg']); end

%% Ipsi Contra Amplitudes

%% rearrange figures
tilefigs([],1,[],[],[],[],[],2)


view_tiff(control_peakramp_AMPA_traces(idxc,:))
view_tiff(MD_peakramp_AMPA_traces(idxm,:))
