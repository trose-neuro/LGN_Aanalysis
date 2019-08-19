% %% LOAD data, you need uipickfiles function and COMBINE SW/MF DATA
% %LOAD MF First
% directory='I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
% filename1=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
% load(char(filename1));%load mat file
% MFdata=LGN;
% %% LOAD SW DATA
% filename2=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
% load(char(filename2));%load mat file
% SWdata=LGN;
% data=[SWdata;MFdata];
%
% %% OR load JB Data
% directory='I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
% load([directory '\Full_Data.mat'])

%% Plot spreadsheet for desired cells
%Define idx of cells
% cellnr=[10 32 46 47 48 55 77 123 205 206 208 210];
% cellnr=[1:12];
% 
% % cellnr=[1];
% 
% for i=1:length(cellnr)
%     cellidx=cellnr(i);
%     if data(cellidx).experimentator=='SW'
%         spreadsheetLGN(data,cellidx,1)
%     else
%         spreadsheetLGN(data,cellidx,2)
%     end
% end

clearvars -except data

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

%%
bins = 7;
binsOD  = -1+(2/bins/2):2/bins:1; % bin centers
binsOD  = -1:2/bins:1; % bin edges

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
%Plot fractions
figure;
test=[fract;fract2]
bar(test);
ylabel('Fraction');
set(gca,'XTickLabel',{'AMPA ODI bins','NMDA ODI bins'});
set(gca,'box','off');
set(gcf,'color','w');
axis square;
ylim([0 1])

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
bar(test);
ylabel('Fraction');
set(gca,'XTickLabel',{'AMPA ODI bins','NMDA ODI bins'});
set(gca,'box','off');
set(gcf,'color','w');
axis square;
ylim([0 1])

%% BINO fractions - TR19

control_binofrac = size(find(abs(ODI_AMPA_r)<1),2) / size(ODI_AMPA_r,2)
MD_binofrac = size(find(abs(ODI_AMPA_r_md)<1),2) / size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2)

n1 = size(find(abs(ODI_AMPA_r)<1),2);
N1 = size(ODI_AMPA_r,2);
n2 = size(find(abs(ODI_AMPA_r_md)<1),2);
N2 = size(ODI_AMPA_r_md(~isnan(ODI_AMPA_r_md)),2);

[pval chi2stat tbl] = propstat(n1, N1, n2, N2);


%% Trace view



%% Control fiber amplitude
[c_A_favg c_N_favg i_A_favg i_N_favg c_A_irravg c_N_irravg i_A_irravg i_N_irravg c_A_f c_N_f i_A_f i_N_f c_A_irr c_N_irr i_A_irr i_N_irr] = failure_extr_ns(data_cont);
diff_irr=abs(c_A_irravg)-abs(i_A_irravg);
diff_irr_N=abs(c_N_irravg)-abs(i_N_irravg);
%% MD fiber amplitude
[c_A_favg_md c_N_favg_md i_A_favg_md i_N_favg_md c_A_irravg_md c_N_irravg_md i_A_irravg_md i_N_irravg_md c_A_f_md c_N_f_md i_A_f_md i_N_f_md c_A_irr_md c_N_irr_md i_A_irr_md i_N_irr_md] = failure_extr_ns(data_md);
diff_irr_md=abs(c_A_irravg_md)-abs(i_A_irravg_md);
diff_irr_N_md=abs(c_N_irravg_md)-abs(i_N_irravg_md);
%% Categories
%Control
Bino_f=find(abs(category==3));
Contra_f=find(abs(category==1));
Ipsi_f=find(abs(category==2));
Ipsis_f=find(abs(category==4));
Contras_f=find(abs(category==5));
%MD
Bino_f_md=find(abs(category_md==3));
Contra_f_md=find(abs(category_md==1));
Ipsi_f_md=find(abs(category_md==2));
Ipsis_f_md=find(abs(category_md==4));
Contras_f_md=find(abs(category_md==5));
%% Define different fiber strengths and read out irradaiance and release probability
%call function fiber_strength
%CONTROL
[allcurra allcurra_i allcurra_NMDA allcurra_NMDA_i low_fibampa high_fibampa low_fibampa_i high_fibampa_i low_fibnmda high_fibnmda low_fibnmda_i ...
    high_fibnmda_i prob_curra prob_curra_i prob_currn prob_curra_ni] = fiber_strength(c_A_f, c_N_f, i_A_f, i_N_f);
%% MD
[allcurra_md allcurra_i_md allcurra_NMDA_md allcurra_NMDA_i_md low_fibampa_md high_fibampa_md low_fibampa_i_md high_fibampa_i_md low_fibnmda_md high_fibnmda_md low_fibnmda_i_md ...
    high_fibnmda_i_md prob_curra_md prob_curra_i_md prob_currn_md prob_curra_ni_md] = fiber_strength(c_A_f_md, c_N_f_md, i_A_f_md, i_N_f_md);
%% Max ramp amplitudes for fiber fraction
%call function maxramp
%Control
[c_A_rmax c_N_rmax i_A_rmax i_N_rmax] = maxramp_ns(data_cont);
%MD
[c_A_rmax_md c_N_rmax_md i_A_rmax_md i_N_rmax_md] = maxramp_ns(data_md);
%% Plot cumulative sum of max currents
figure;cdfplot(c_A_rmax);hold on; cdfplot(i_A_rmax);
hold on; cdfplot(c_A_rmax_md);hold on; cdfplot(i_A_rmax_md);
axis square;
grid off;
title off;
legend('contra', 'ipsi','contra MD', 'ipsi MD');
xlabel('Synaptic input');
set(gca,'box','off');
set(gcf,'color','w');

%% Plot cumulative sum of max failures
%CDF Plot AMPA
figure;cdfplot(abs(vertcat(allcurra{1,:})));hold on; cdfplot(abs(vertcat(allcurra_i{1,:})));
hold on;cdfplot(abs(vertcat(allcurra_md{1,:})));hold on; cdfplot(abs(vertcat(allcurra_i_md{1,:})));
axis square;
grid off;
axis square;
legend('contra', 'ipsi','contra MD','ipsi MD');
xlabel('Synaptic input (pA)');
ylabel('Cumulative Probability');
set(gcf,'color','w');
title('');
xlim([0 1000]);
legend('boxoff');
set(gca,'box','off');
axis square;
%NMDA
figure; cdfplot(vertcat(allcurra_NMDA{1,:}));hold on; cdfplot(vertcat(allcurra_NMDA_i{1,:}));
hold on;cdfplot(vertcat(allcurra_NMDA_md{1,:}));hold on; cdfplot(vertcat(allcurra_NMDA_i_md{1,:}))
axis square;
grid off;
axis square;
legend('contra', 'ipsi','contra MD','ipsi MD');
xlabel('Synaptic input (pA)');
ylabel('Cumulative Probability');
set(gcf,'color','w');
title('');
xlim([0 500]);
legend('boxoff');
set(gca,'box','off');
axis square;
%% % Whisker plot Contra/Ipsi strength
%using the average values
allcur_coc=[nonzeros(c_A_favg);nonzeros(i_A_favg);nonzeros(c_N_favg);nonzeros(i_N_favg)];
allcur_cat=[ones(length(nonzeros(c_A_favg)),1);ones(length(nonzeros(i_A_favg)),1)*2;ones(length(nonzeros(c_N_favg)),1)*3;ones(length(nonzeros(i_N_favg)),1)*4];
allcur_coc_md=[nonzeros(c_A_favg_md);nonzeros(i_A_favg_md);nonzeros(c_N_favg_md);nonzeros(i_N_favg_md)];
allcur_cat_md=[ones(length(nonzeros(c_A_favg_md)),1)*5;ones(length(nonzeros(i_A_favg_md)),1)*6;ones(length(nonzeros(c_N_favg_md)),1)*7;ones(length(nonzeros(i_N_favg_md)),1)*8];
%subdivide into bino ipsi bino contra
%Plot whisker
figure;hold on;p1=plot(ones(length(nonzeros(c_A_favg)),1),abs(nonzeros(c_A_favg)),'ko','MarkerSize',3);p1.Color(4) = 0.01;
hold on;p2=plot(ones(length(nonzeros(i_A_favg)),1)*2,abs(nonzeros(i_A_favg)),'ko','MarkerSize',3);p2.Color(4) = 0.01;
hold on;p3=plot(ones(length(nonzeros(c_N_favg)),1)*3,abs(nonzeros(c_N_favg)),'ko','MarkerSize',3);p3.Color(4) = 0.01;
hold on;p4=plot(ones(length(nonzeros(i_N_favg)),1)*4,abs(nonzeros(i_N_favg)),'ko','MarkerSize',3);p4.Color(4) = 0.01;
hold on;boxplot(allcur_coc,allcur_cat,'Notch','on','Labels',{'Contra AMPAf','Ipsi AMPAf','Contra NMDAf','Ipsi NMDAf' },'Whisker',1, 'OutlierSize', 1)
axis square;
set(gcf,'color','w');
set(gca,'box','off');
ylabel('Synaptic input (pA)');
xtickangle(45);
xlim([0 5]);
ylim([0 280]);
figure;
hold on;p5=plot(ones(length(nonzeros(c_A_favg_md)),1)*1,abs(nonzeros(c_A_favg_md)),'ko','MarkerSize',3);p5.Color(4) = 0.01;
hold on;p6=plot(ones(length(nonzeros(i_A_favg_md)),1)*2,abs(nonzeros(i_A_favg_md)),'ko','MarkerSize',3);p6.Color(4) = 0.01;
hold on;p7=plot(ones(length(nonzeros(c_N_favg_md)),1)*3,abs(nonzeros(c_N_favg_md)),'ko','MarkerSize',3);p7.Color(4) = 0.01;
hold on;p8=plot(ones(length(nonzeros(i_N_favg_md)),1)*4,abs(nonzeros(i_N_favg_md)),'ko','MarkerSize',3);p8.Color(4) = 0.01;
boxplot(allcur_coc_md,allcur_cat_md,'Notch','on','Labels',{'Contra AMPAf MD','Ipsi AMPAf MD','Contra NMDAf MD','Ipsi NMDAf MD' },'Whisker',1, 'OutlierSize', 1)
axis square;
set(gcf,'color','w');
set(gca,'box','off');
ylabel('Synaptic input (pA)');
xtickangle(45);
xlim([0 5]);
ylim([0 280]);

%% Fiber number
%Control
%using average values
Fib_AMPA_f_c=(c_A_rmax./c_A_favg);
Fib_AMPA_f_i=(i_A_rmax./i_A_favg);
Fib_NMDA_f_c=(c_N_rmax./c_N_favg);
Fib_NMDA_f_i=(i_N_rmax./i_N_favg);
%inverted
Fib_AMPA_f_c_in=(c_A_favg./c_A_rmax);
Fib_AMPA_f_i_in=(i_A_favg./i_A_rmax);
Fib_NMDA_f_c_in=(c_N_favg./c_N_rmax);
Fib_NMDA_f_i_in=(i_N_favg./i_N_rmax);
%MD
%using average values
Fib_AMPA_f_c_md=(c_A_rmax_md./c_A_favg_md);
Fib_AMPA_f_i_md=(i_A_rmax_md./i_A_favg_md);
Fib_NMDA_f_c_md=(c_N_rmax_md./c_N_favg_md);
Fib_NMDA_f_i_md=(i_N_rmax_md./i_N_favg_md);
%inverted
Fib_AMPA_f_c_in_md=(c_A_favg_md./c_A_rmax_md);
Fib_AMPA_f_i_in_md=(i_A_favg_md./i_A_rmax_md);
Fib_NMDA_f_c_in_md=(c_N_favg_md./c_N_rmax_md);
Fib_NMDA_f_i_in_md=(i_N_favg_md./i_N_rmax_md);

%% Categorize fiber numer/fraction according to input category
%average fiber numbers
%AMPA
bin=Fib_AMPA_f_i(Bino_f)+Fib_AMPA_f_c(Bino_f);
bin_ipsi=Fib_AMPA_f_i(Bino_f);
bin_contra=Fib_AMPA_f_c(Bino_f);
con=Fib_AMPA_f_c(Contra_f);
ipsi=Fib_AMPA_f_i(Ipsi_f);
ipsi_sil=Fib_AMPA_f_c(Ipsis_f);
con_sil=Fib_AMPA_f_i(Contras_f);
%NMDA
bin_NMDA=Fib_NMDA_f_i(Bino_f)+Fib_NMDA_f_c(Bino_f);
bin_ipsi_NMDA=Fib_NMDA_f_i(Bino_f);
bin_contra_NMDA=Fib_NMDA_f_c(Bino_f);
con_NMDA=Fib_NMDA_f_c(Contra_f);
ipsi_NMDA=Fib_NMDA_f_i(Ipsi_f);
ipsi_sil_NMDA=Fib_NMDA_f_i(Ipsis_f);
con_sil_NMDA=Fib_NMDA_f_c(Contras_f);
%average fiber fraction
%AMPA
bin_in=Fib_AMPA_f_i_in(Bino_f)+Fib_AMPA_f_c_in(Bino_f);
bin_ipsi_in=Fib_AMPA_f_i_in(Bino_f);
bin_contra_in=Fib_AMPA_f_c(Bino_f);
con_in=Fib_AMPA_f_c_in(Contra_f);
ipsi_in=Fib_AMPA_f_i_in(Ipsi_f);
ipsi_sil_in=Fib_AMPA_f_c_in(Ipsis_f);
con_sil_in=Fib_AMPA_f_i_in(Contras_f);
%NMDA
bin_in_NMDA=Fib_NMDA_f_i_in(Bino_f)+Fib_NMDA_f_c_in(Bino_f);
bin_ipsi_in_NMDA=Fib_NMDA_f_i_in(Bino_f);
bin_contra_in_NMDA=Fib_NMDA_f_c_in(Bino_f);
con_in_NMDA=Fib_NMDA_f_c_in(Contra_f);
ipsi_in_NMDA=Fib_NMDA_f_i_in(Ipsi_f);
ipsi_sil_in_NMDA=Fib_NMDA_f_c_in(Ipsis_f);
con_sil_in_NMDA=Fib_NMDA_f_i_in(Contras_f);

%MD
%AMPA
bin_md=Fib_AMPA_f_i_md(Bino_f_md)+Fib_AMPA_f_c_md(Bino_f_md);
bin_ipsi_md=Fib_AMPA_f_i_md(Bino_f_md);
bin_contra_md=Fib_AMPA_f_c_md(Bino_f_md);
con_md=Fib_AMPA_f_c_md(Contra_f_md);
ipsi_md=Fib_AMPA_f_i_md(Ipsi_f_md);
ipsi_sil_md=Fib_AMPA_f_c_md(Ipsis_f_md);
con_sil_md=Fib_AMPA_f_i_md(Contras_f_md);
%NMDA
bin_NMDA_md=Fib_NMDA_f_i_md(Bino_f_md)+Fib_NMDA_f_c_md(Bino_f_md);
bin_ipsi_NMDA_md=Fib_NMDA_f_i_md(Bino_f_md);
bin_contra_NMDA_md=Fib_NMDA_f_c_md(Bino_f_md);
con_NMDA_md=Fib_NMDA_f_c_md(Contra_f_md);
ipsi_NMDA_md=Fib_NMDA_f_i_md(Ipsi_f_md);
ipsi_sil_NMDA_md=Fib_NMDA_f_i_md(Ipsis_f_md);
con_sil_NMDA_md=Fib_NMDA_f_c_md(Contras_f_md);
%average fiber fraction
% %AMPA
% bin_in=Fib_AMPA_f_i_in(Bino_f)+Fib_AMPA_f_c_in(Bino_f);
% bin_ipsi_in=Fib_AMPA_f_i_in(Bino_f);
% bin_contra_in=Fib_AMPA_f_c(Bino_f);
% con_in=Fib_AMPA_f_c_in(Contra_f);
% ipsi_in=Fib_AMPA_f_i_in(Ipsi_f);
% ipsi_sil_in=Fib_AMPA_f_c_in(Ipsis_f);
% con_sil_in=Fib_AMPA_f_i_in(Contras_f);
% %NMDA
% bin_in_NMDA=Fib_NMDA_f_i_in(Bino_f)+Fib_NMDA_f_c_in(Bino_f);
% bin_ipsi_in_NMDA=Fib_NMDA_f_i_in(Bino_f);
% bin_contra_in_NMDA=Fib_NMDA_f_c_in(Bino_f);
% con_in_NMDA=Fib_NMDA_f_c_in(Contra_f);
% ipsi_in_NMDA=Fib_NMDA_f_i_in(Ipsi_f);
% ipsi_sil_in_NMDA=Fib_NMDA_f_c_in(Ipsis_f);
% con_sil_in_NMDA=Fib_NMDA_f_i_in(Contras_f);
%% Average and combine
%Contol
com_AMPAf=[nanmean(con) nanmean(ipsi) nanmean(bin) nanmean(bin_contra) nanmean(bin_ipsi) nanmean(ipsi_sil) nanmean(con_sil)];
sem_AMPAf=[nanstd(con)/sqrt(length(con)) nanstd(ipsi)/sqrt(length(ipsi)) nanstd(bin)/sqrt(length(bin)) nanstd(bin_contra)/sqrt(length(bin_contra)) nanstd(bin_ipsi)/sqrt(length(bin_ipsi)) nanstd(ipsi_sil)/sqrt(length(ipsi_sil)) nanstd(con_sil)/sqrt(length(con_sil))];

com_NMDAf=[nanmean(con_NMDA) nanmean(ipsi_NMDA) nanmean(bin_NMDA) nanmean(bin_contra_NMDA) nanmean(bin_ipsi_NMDA) nanmean(ipsi_sil_NMDA) nanmean(con_sil_NMDA)];
sem_NMDAf=[nanstd(con_NMDA)/sqrt(length(con_NMDA)) nanstd(ipsi_NMDA)/sqrt(length(ipsi_NMDA)) nanstd(bin_NMDA)/sqrt(length(bin_NMDA)) nanstd(bin_contra_NMDA)/sqrt(length(bin_contra_NMDA)) nanstd(bin_ipsi_NMDA)/sqrt(length(bin_ipsi_NMDA)) nanstd(ipsi_sil_NMDA)/sqrt(length(ipsi_sil_NMDA)) nanstd(con_sil_NMDA)/sqrt(length(con_sil_NMDA))];
%MD
com_AMPAf_md=[nanmean(con_md) nanmean(ipsi_md) nanmean(bin_md) nanmean(bin_contra_md) nanmean(bin_ipsi_md) nanmean(ipsi_sil_md) nanmean(con_sil_md)];
sem_AMPAf_md=[nanstd(con_md)/sqrt(length(con_md)) nanstd(ipsi_md)/sqrt(length(ipsi_md)) nanstd(bin_md)/sqrt(length(bin_md)) nanstd(bin_contra_md)/sqrt(length(bin_contra_md)) nanstd(bin_ipsi_md)/sqrt(length(bin_ipsi_md)) nanstd(ipsi_sil_md)/sqrt(length(ipsi_sil_md)) nanstd(con_sil_md)/sqrt(length(con_sil_md))];

com_NMDAf_md=[nanmean(con_NMDA_md) nanmean(ipsi_NMDA_md) nanmean(bin_NMDA_md) nanmean(bin_contra_NMDA_md) nanmean(bin_ipsi_NMDA_md) nanmean(ipsi_sil_NMDA_md) nanmean(con_sil_NMDA_md)];
sem_NMDAf_md=[nanstd(con_NMDA_md)/sqrt(length(con_NMDA_md)) nanstd(ipsi_NMDA_md)/sqrt(length(ipsi_NMDA_md)) nanstd(bin_NMDA_md)/sqrt(length(bin_NMDA_md)) nanstd(bin_contra_NMDA_md)/sqrt(length(bin_contra_NMDA_md)) nanstd(bin_ipsi_NMDA_md)/sqrt(length(bin_ipsi_NMDA_md)) nanstd(ipsi_sil_NMDA_md)/sqrt(length(ipsi_sil_NMDA_md)) nanstd(con_sil_NMDA_md)/sqrt(length(con_sil_NMDA_md))];

%% Plot bar graph for fiber numbers
y_fin=[nanmean(con) nanmean(con_md); nanmean(ipsi) nanmean(ipsi_md); nanmean(bin) nanmean(bin_md);...
    nanmean(bin_contra) nanmean(bin_contra_md); nanmean(bin_ipsi) nanmean(bin_ipsi_md);...
    nanmean(ipsi_sil) nanmean(ipsi_sil_md);nanmean(con_sil) nanmean(con_sil_md)];
y_sem=[sem_AMPAf(1) sem_AMPAf_md(1);sem_AMPAf(2) sem_AMPAf_md(2);sem_AMPAf(3) sem_AMPAf_md(3);...
    sem_AMPAf(4) sem_AMPAf_md(4);sem_AMPAf(5) sem_AMPAf_md(5);sem_AMPAf(6) sem_AMPAf_md(6);...
    sem_AMPAf(7) sem_AMPAf_md(7)];

figure;
ctrs = 1:7;
data = y_fin;
b=y_sem;
figure(1)
hBar = bar(ctrs, data);
for k1 = 1:size(y_fin,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, b', '.k')
ylabel('Fiber number');
set(gca,'XTickLabel',{'contra', 'ipsi','bino','bino contra','bino ipsi','ipsi silent' 'contra silent'});
xtickangle(45);
axis square;
set(gcf,'color','w');
set(gca,'box','off');
ylim([0 100]);
legend('Control', 'MD')


%% Plot bar graph for fiber numbers NMDA
y_fin_NMDA=[com_NMDAf(1) com_NMDAf_md(1);com_NMDAf(2) com_NMDAf_md(2);com_NMDAf(3) com_NMDAf_md(3);...
    com_NMDAf(4) com_NMDAf_md(4);com_NMDAf(5) com_NMDAf_md(5);com_NMDAf(6) com_NMDAf_md(6);...
    com_NMDAf(7) com_NMDAf_md(7)];
y_NMDA_sem=[sem_NMDAf(1) sem_NMDAf_md(1);sem_NMDAf(2) sem_NMDAf_md(2);sem_NMDAf(3) sem_NMDAf_md(3);...
    sem_NMDAf(4) sem_NMDAf_md(4);sem_NMDAf(5) sem_NMDAf_md(5);sem_NMDAf(6) sem_NMDAf_md(6);...
    sem_NMDAf(7) sem_NMDAf_md(7)];

figure
ctrs = 1:7;
data = y_fin_NMDA;
b=y_NMDA_sem;
figure(1)
hBar = bar(ctrs, data);
for k1 = 1:size(y_fin_NMDA,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, b', '.k')
ylabel('Fiber number');
set(gca,'XTickLabel',{'contra', 'ipsi','bino','bino contra','bino ipsi','ipsi silent' 'contra silent'});
xtickangle(45);
axis square;
set(gcf,'color','w');
set(gca,'box','off');
ylim([0 100]);
legend('Control', 'MD')

%%


