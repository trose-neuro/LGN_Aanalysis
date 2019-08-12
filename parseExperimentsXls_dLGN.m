function [batchopt] = parseExperimentsXls_dLGN(path,user)

[xls_num,xls_txt]=xlsread(path);

loadcol        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'BatchAnalyze')));
mousecol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExperimentalDay')));
expcol         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsSW')));
expcols2       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsSW_slice_nr')));
category       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Category_SW')));
injection_order= find(~cellfun(@isempty, strfind(xls_txt(1,:),'contra_ipsi_SW')));
tracing        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Tracing_SW')));
morphox        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'xf_SW')));
morphoy        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'yf_SW')));
morphoz        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'zf_SW')));
MD             = find (~cellfun(@isempty, strfind(xls_txt(1,:),'MD_contra_ipsi_SW')));
Clear_sl_nr    = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Cleared_brain_slice_number_SW')));
Hemisphere     = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Hemisphere_SW')));
eye_inj_order  = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Injection_order'))); 
photodiode_flag = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Photodiode signal flag SW'))); 



if user==1
expcol          = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsMF')));
expcols2        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingsMF_slice_nr')));
category        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Category_MF')));
injection_order = find(~cellfun(@isempty, strfind(xls_txt(1,:),'contra_ipsi_MF')));
tracing         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Tracing_MF')));
morphox         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'xf_MF')));
morphoy         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'yf_MF')));
morphoz         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'zf_MF')));
MD              = find (~cellfun(@isempty, strfind(xls_txt(1,:),'MD_contra_ipsi_MF')));
Clear_sl_nr     = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Cleared_brain_slice_number_MF')));
Hemisphere      = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Hemisphere_MF')));
eye_inj_order   = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Injection_order')));
photodiode_flag = find (~cellfun(@isempty, strfind(xls_txt(1,:),'Photodiode signal flag MF'))); 
end

loaddrivecol  = find(~cellfun(@isempty, strfind(xls_txt(1,:),'loaddrive')));
animalname    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Animal_ID')));


k = 1;

batchopt.XLS.txt = xls_txt;
batchopt.XLS.num = xls_txt;

for i = 2:size(xls_txt,1)
    ana{k}= xls_num(i-1,loadcol-1);
    
    
    
    if ~ana{k}
        disp(['skipping experiments ' xls_txt{i,mousecol} ' (no batchload flag)']);
        continue
    end
    
    batchopt.mouse{k} = xls_txt(i,mousecol);
    batchopt.mouseID{k} = xls_txt(i,animalname);
    batchopt.eye_inj_order{k} = xls_num(i-1,eye_inj_order-1);
    
    expcellids{k}                = xls_txt(i,expcol);
    expcellids2{k}               = xls_txt(i,expcols2);
    expcellids3{k}               = xls_txt(i,category);
    expcellids4{k}               = xls_txt(i,injection_order);
    expcellids5{k}               = xls_txt(i,tracing);
    expcellids6{k}               = xls_txt(i,morphox);
    expcellids7{k}               = xls_txt(i,morphoy);
    expcellids8{k}               = xls_txt(i,morphoz);
    expcellids9{k}               = xls_txt(i,MD);
    expcellids10{k}              = xls_txt(i,Clear_sl_nr);
    expcellids11{k}              = xls_txt(i,Hemisphere);
    expcellids13{k}              = xls_txt(i,photodiode_flag);
   % expcellids12{k}              = xls_txt(i,eye_inj_order);
    
    %     spontcellids{k}                = xls_txt(i,spontcol);
    %     sftfcellids{k}                = xls_txt(i,sftfcol);
    %     puffcellids{k}                = xls_txt(i,airpuffcol);
    
    %batchopt.exp_ids{k}          = (expcellids{k}{1});
    
    
    batchopt.exp_ids{k}         = str2num((expcellids{k}{1}));
    batchopt.slice{k}           = str2num((expcellids2{k}{1}));
    batchopt.category{k}        = str2num((expcellids3{k}{1}));
    batchopt.injection_order{k} = str2num((expcellids4{k}{1}));
    batchopt.tracing{k}         = str2num((expcellids5{k}{1}));
    
    batchopt.morphox{k}         = str2num((expcellids6{k}{1}));
    batchopt.morphoy{k}         = str2num((expcellids7{k}{1}));
    batchopt.morphoz{k}         = str2num((expcellids8{k}{1}));
    batchopt.MD{k}              = str2num((expcellids9{k}{1}));
    batchopt.clear_sl_nr{k}     = str2num((expcellids10{k}{1}));
    batchopt.hemisphere{k}      = (expcellids11{k}{1});
    batchopt.photodiode_flag{k}      = str2num((expcellids13{k}{1}));
    
    %batchopt.eye_inj_order{k}  = str2num((expcellids12{k}{1}));
    %     batchopt.spont_ids{k}          = str2num((spontcellids{k}{1}));
    %     batchopt.sftf_ids{k}          = str2num((sftfcellids{k}{1}));
    %     batchopt.puff_ids{k}          = str2num((puffcellids{k}{1}));
    %
    %     batchopt.samesite{k}         = xls_num(i-1,samesitecol-1);
    %     batchopt.baseline{k}         = xls_num(i-1,baselinecol(1)-1);
    %     try
    %         batchopt.baselinepair{k}     = eval(cell2mat(xls_txt(i,basepaircol(1))));
    %         batchopt.baselinepair14{k}   = eval(cell2mat(xls_txt(i,basepair14col)));
    %         batchopt.recovery{k}         = xls_num(i-1,recoverycol-1);
    %     end
    batchopt.loaddrive{k}        = xls_txt(i,loaddrivecol);
    k = k+1;
end

