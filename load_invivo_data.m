
%% load data

load('/Users/trose/Documents/GitHub/LGN_Analysis/Full_Data.mat')
load('/Users/trose/Documents/GitHub/LGN_Analysis/in vivo/base_MD_ODI_array.mat')
load('/Users/trose/Documents/GitHub/LGN_Analysis/in vivo/ana_ODI_Musc.mat')

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

opengl software         % to get the axes on figures with transparency back!
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultAxesFontSize', fig.fsz)
set(0,'DefaultAxesTickLength', [0.02 0.025]);

%% histogram settings
bins = 7;
binsOD  = -1:2/bins:1; % bin edges

%% FIND CONTROL AND MD idx
field_number = size(data,2);
md=find([data(:).MD]==1);
cont=find([data(:).MD]==0);

%% Category read out for Control and MD
%Control
for i=1:length(cont)
    %Category read out
    category(i)=data(cont(i)).ocular_category;
    ODI_AMPA_r(i)=data(cont(i)).ODI_AMPA_step_peak;
    ODI_NMDA_r(i)=data(cont(i)).ODI_NMDA_step_peak;
    data_cont(i)=data(cont(i));
end
%MD
for i=1:length(md)
    %Category read out
    category_md(i)=data(md(i)).ocular_category;
    ODI_AMPA_r_md(i)=data(md(i)).ODI_AMPA_step_peak;
    ODI_NMDA_r_md(i)=data(md(i)).ODI_NMDA_step_peak;
    data_md(i)=data(md(i));
end

%% CONTROL
%AMPA fraction discretize
temp=ODI_AMPA_r;
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    fract(i)=(length(find(disc==i)))/length(temp)
end
disc=[];
%NMDA fraction discretize
temp2=ODI_NMDA_r;
temp2(find(isnan(temp2)))=[];
disc=discretize(temp2,[binsOD]);
for i=1:bins
    fract2(i)=(length(find(disc==i)))/length(temp2)
end
disc=[];

%% MD

%AMPA fraction discretize
temp=ODI_AMPA_r_md;
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    fract(i)=(length(find(disc==i)))/length(temp)
end
disc=[];
%NMDA fraction discretize
temp2=ODI_NMDA_r_md;
temp2(find(isnan(temp2)))=[];
disc=discretize(temp2,[binsOD]);
for i=1:bins
    fract2(i)=(length(find(disc==i)))/length(temp2)
end
disc=[];
%Plot fractions
figure;
test=[fract;fract2]
AMPA_ODI_MD_fract = test(2,1);

%% BINO fractions - TR19

control_binofrac = size(find(abs(ODI_AMPA_r)<1),2) / size(ODI_AMPA_r,2)
MD_binofrac = size(find(abs(ODI_AMPA_r_md)<1),2) / size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2)

n1 = size(find(abs(ODI_AMPA_r)<1),2);
N1 = size(ODI_AMPA_r,2);
n2 = size(find(abs(ODI_AMPA_r_md)<1),2);
N2 = size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

[pval chi2stat tbl] = propstat(n1, N1, n2, N2);

%% compare to in vivop ODIs
load('/Users/trose/Documents/GitHub/LGN_Analysis/in vivo/base_MD_ODI_array.mat')
bins = 7
binsOD  = -1:2/bins:1; % bin edges
%AMPA fraction discretize
% baseline
temp=BaseMD_ODI_array(:,1);
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    fract(i)=(length(find(disc==i)))/length(temp)
end
disc=[];

% MD
temp=BaseMD_ODI_array(:,2);
temp(find(isnan(temp)))=[];
disc=discretize(temp,[binsOD]);
for i=1:bins
    fract2(i)=(length(find(disc==i)))/length(temp)
end
disc=[];


Invivo_base_fract = fract;
Invivo_MD_fract = fract2;

%Plot fractions
figure;
test=[Invivo_base_fract;AMPA_ODI_baseline_fract]

bar(test);
ylabel('Fraction');
set(gca,'XTickLabel',{'AMPA ODI bins','NMDA ODI bins'});
set(gca,'box','off');
set(gcf,'color','w');
axis square;
ylim([0 1])


