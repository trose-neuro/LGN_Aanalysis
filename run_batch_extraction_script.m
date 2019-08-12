% initial traces & morphology processing
clear all
close all
user = 0; % SW
Analysis_mini_ramp

clear all 
close all
user = 1; % MF
Analysis_mini_ramp

% combine final structure and add localization info
clear all 
close all
Add_loc_info_to_DATAstructure

%
Axon_ratio_calc_save('I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2', 'JB_181126_1','Martin',1,'190124LH')


% reload data
clear all 
close all
directory='D:\LGN project'; % use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
load([directory '\Full_Data.mat'])

% overview step protocol ODI distribution
AMPA_NMDA_stepprot_ODI_plot(data,'color_cat')

% brows cells
cellnr=[1:length(data)];
cellnr=find([data(:).ocular_category]==1 & [data(:).ODI_AMPA_step]<1);

for i=1:length(cellnr)
    cellidx=cellnr(i);
    spreadsheetLGN(data,cellidx)
    pause
    close
end
disp('End of cell series')

find([data(:).ocular_category]==1 & [data(:).ODI_AMPA_step]<1)

spreadsheetLGN(data,16)


