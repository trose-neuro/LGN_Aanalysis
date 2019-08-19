JData =  '/Volumes/archive_bonhoeffer_group$-1/Juliane Jaepel-Schael/Shared/master_Juliane.m'

load(JData, '-mat');

saver = analysis_MD_dLGN(1,master_selector, fig, batchopt, saver)