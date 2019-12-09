function  [resp_maps, traces, illumination, resp_windows]=LGN_scracm_new(list, idx, voltage, pathName, fc, show, ramp_rtrace, user, filterephys,adata_dir,map_trace);



%define temporal windows
base_start          =   1000+1;%1;
base_end            =   1000+99;%99;

redpeak_start       =   1000+100;%100;
redpeak_end         =   1000+199;%350; %looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   1000+200;%351;
bluepeak_end        =   1000+299;%400;

% base_start2          =   1000;+%bluepeak_start-100;
% base_end2            =   1000;+%bluepeak_start-2;

% window_ampa = bluepeak_end-bluepeak_start;
% window_nmda = bluepeak_end-bluepeak_start;
window_ampa = 98;
window_nmda = 98;

BL_AMPA_P1 = [redpeak_start-1-window_ampa,   redpeak_start-1];
BL_AMPA_P2 = [bluepeak_start-1-window_ampa,  bluepeak_start-1];
BL_NMDA_P1 = [redpeak_start-1-window_nmda,   redpeak_start-1];
BL_NMDA_P2 = [bluepeak_start-1-window_nmda,  bluepeak_start-1];

RW_AMPA_P1 = [redpeak_start+3,               redpeak_start+3+window_ampa];
RW_AMPA_P2 = [bluepeak_start+3,              bluepeak_start+3+window_ampa];
RW_NMDA_P1 = [redpeak_start+3,               redpeak_start+3+window_nmda];
RW_NMDA_P2 = [bluepeak_start+3,              bluepeak_start+3+window_nmda];

% BL_AMPA_P1 = [redpeak_start-28 redpeak_start-1];
% BL_AMPA_P2 = [bluepeak_start-28 bluepeak_start-1];
% BL_NMDA_P1 = [redpeak_start-99 redpeak_start-1];
% BL_NMDA_P2 = [bluepeak_start-99 bluepeak_start-1];
% 
% RW_AMPA_P1 = [redpeak_start+3 redpeak_start+30];
% RW_AMPA_P2 = [bluepeak_start+3 bluepeak_start+30];
% RW_NMDA_P1 = [redpeak_start+3 redpeak_start+107];
% RW_NMDA_P2 = [bluepeak_start+3 bluepeak_start+107];

%IMPORTANT
%neg_peak1: atm shortened the window!
%% TR2019: filtering
% filterephys = 1;        % filtering yes/no?
cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel. Use Bessel at > 4 order to prevent ripples)

if filterephys;
    disp('- - - - - - - -')
    disp(['Filtering: ' num2str(order) ' pole ' type '-Filter w/ ' num2str(cutoff) ' Hz cutoff']);
    disp('- - - - - - - -')
end

if user==0%SW
    vclamp=voltage(idx);
    AMPArep = 0;
    NMDArep = 0;
    for j=1:length(idx)% how many maps in total; loop across ramps per cell
        counter=1;
        %         if show==1
        %             subplot(2,(length(idx))-2,j);
        %         end
        load([char(pathName) filesep list(idx(j)).name],'-mat');
        sr = header.ephys.ephys.sampleRate;%check sample rate
        srF = 1/(1000/sr);
        
        samples_per_sweep = header.ephys.ephys.traceLength*sr;
        timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
        
        traces=data.ephys.trace_1;%raw ephys trace
        if filterephys % TR2019: filtering
            traces = lowpassfilt(traces, order, cutoff, sr, type);
        end
        
        photodiode=data.acquirer.trace_1;
%         ind_traces=reshape(traces,[length(traces)/64 64]);
        ind_traces=reshape(traces,[length(traces)/49 49]);
                
        baseline=ind_traces(base_start*srF:base_end*srF,:);
        
        %%
        
        mapPattern=header.mapper.mapper.mapPatternArray;
        if 0%map_trace
            plot_stim_resp_map(ind_traces,baseline,mapPattern,redpeak_start,bluepeak_start,bluepeak_end,srF)
        end
        
        %%
        
        %         photodiode=reshape(photodiode,[length(traces)/64 64]);
        photodiode=reshape(photodiode,[length(traces)/49 49]);
        bs_photodiode=photodiode-mean(photodiode(base_start*srF:base_end*srF,:));
        double_pulse(:,j)=mean(bs_photodiode((redpeak_start+20)*srF:(redpeak_end-20)*srF))>0.025;
        %         double_pulse(:,j)=mean(mean(photodiode(20100:20400,1:10)))>1;% find which traces are single vs double laser stimulation
        %         dp=mean(mean(photodiode(20100:20400,1:10)))>1;% find which traces are single vs double laser stimulation per map
        
        soma_position=header.mapper.mapper.soma1Coordinates;
        %MAP PATTERN FOR EACH ACQUISTION
        
        map_pattern=header.mapper.mapper.mapPatternArray;
        map_pattern_all(:,:,j)=header.mapper.mapper.mapPatternArray;
        %         map_negpeak1=map_pattern;
        %         map_negpeak2=map_pattern;
        %         map_pospeak1=map_pattern;
        %         map_pospeak2=map_pattern;
        %         map_negfail1=map_pattern;
        %         map_negfail2=map_pattern;
        %         map_posfail1=map_pattern;
        %         map_posfail2=map_pattern;
        
        
        for i=1:size(ind_traces,2);
            fitting_used = 0;
            
            traces_clip=ind_traces(:,i);
            bs=traces_clip(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
            bs_std=std(bs);%std of baseline trace
            bs_traces=traces_clip-mean(traces_clip(base_start*srF:base_end*srF,:));%subtract baseline
            
            % shape and variation of photodiode signal needs to be saved for quality control!
            temp = [mean(bs_photodiode'); std(bs_photodiode')];
            bs_photodiode_mean_std{i,j} = temp(:,1:10:end); % downsample by 10
            PD_pulse1(:,i)= mean(bs_photodiode(RW_AMPA_P1(1) *srF:RW_AMPA_P1(2)*srF,i));
            PD_pulse2(:,i)= mean(bs_photodiode(RW_AMPA_P2(1) *srF:RW_AMPA_P2(2)*srF,i));
            
            [r, c]=find(map_pattern_all(:,:,j)==i);
            if user==0%SW
                PD1_r(r,c,j)=PD_pulse1(:,i);
                PD2_b(r,c,j)=PD_pulse2(:,i);
                IR1_r(r,c,j)=(12.19*PD_pulse1(:,i)-0.4319)/100;
                IR2_b(r,c,j)=(7.232*PD_pulse2(:,i)-0.9951)/100;
            end
            clear PD_pulse1 PD_pulse2
            %% new trace extraction by JB
            
            if vclamp(j) %if AMPA
                pulse1_BL = bs_traces(BL_AMPA_P1(1)*srF:BL_AMPA_P1(2)*srF,:);
                pulse1_RW = bs_traces(RW_AMPA_P1(1)*srF:RW_AMPA_P1(2)*srF,:);
                pulse2_BL = bs_traces(BL_AMPA_P2(1)*srF:BL_AMPA_P2(2)*srF,:);
                pulse2_RW = bs_traces(RW_AMPA_P2(1)*srF:RW_AMPA_P2(2)*srF,:);
                
            else
                pulse1_BL = bs_traces(BL_NMDA_P1(1)*srF:BL_NMDA_P1(2)*srF,:);
                pulse1_RW = bs_traces(RW_NMDA_P1(1)*srF:RW_NMDA_P1(2)*srF,:);
                if  max(pulse1_RW) < fc*bs_std % if no first pulse is detected dont do fitting
                    pulse2_BL = bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:);
                    pulse2_RW = bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:);
                else
                    xt=1:50000; % assuming  SW
                    A=max(pulse1_RW);
                    t1=find(bs_traces==A)+(10*srF); % start 10 after peak
                    t=t1:(redpeak_end-1)*srF;
                    t=t';  % check dif SW MF
                    curr_t=bs_traces(t);
                    
                    try
                        % fit with constraints
                        % single exponential
                        fo = fitoptions('Method','NonlinearLeastSquares', ...
                            'Lower',[0 -1], ...
                            'Upper', [max(pulse1_RW)*10 0],...
                            'StartPoint',[max(pulse1_RW)*3 -0.0005]);
                        ft = fittype('a*exp(b*(x))', 'options', fo);
                        [f1 gof1]=fit(t,curr_t,ft);
                        yf1=f1.a*exp(f1.b*xt);
                        %                     figure; plot(bs_traces); hold on; plot(yf)
                        
                        try % double exponential decay (if it fails or is worse exp1 is used)
                            fo = fitoptions('Method','NonlinearLeastSquares', ...
                                'MaxIter', 1000,...
                                'Lower',[0 -1 0 -1], ...
                                'Upper', [max(pulse1_RW)*10 0 max(pulse1_RW)*10 0],...
                                'StartPoint',[max(pulse1_RW)*1 -0.001 max(pulse1_RW)*1 -0.000005]);
                            ft = fittype('a*exp(b*(x)) + c*exp(d*(x))', 'options', fo);
                            [f2 gof2]=fit(t,curr_t,ft);
                            yf2=f2.a*exp(f2.b*xt) + f2.c*exp(f2.d*xt);
                            %                     figure; plot(bs_traces); hold on; plot(yf)
                        end
                        try
                            if gof2.adjrsquare>gof1.adjrsquare
                                gof = gof2;
                                f = f2;
                                yf = yf2;
                            else
                                gof = gof1;
                                f = f1;
                                yf = yf1;
                            end
                        catch
                            gof = gof1;
                            f = f1;
                            yf = yf1;
                        end
                        
                        diff_bs_traces = bs_traces;
                        for m=1:10000
                            diff_bs_traces(m,:)=bs_traces(m)-yf(m);
                        end
                        
                        if show==1
                            fitfig = figure;
                            plot(bs_traces);hold on;plot(yf);plot(diff_bs_traces);
                            
                            set(fitfig, 'Name', ['FIT:' char(pathName) ]);
                            %%red vertical lines
                            hold on;
                            y1=get(gca,'ylim');
                            x1= redpeak_start*srF;
                            hold on;
                            p1=plot([x1 x1],y1,'--','Color','r');
                            p1.Color(4) = 0.3;
                            hold on;
                            y1=get(gca,'ylim');
                            x1=redpeak_end*srF;
                            hold on;
                            p2=plot([x1 x1],y1,'--','Color','r');
                            p2.Color(4) = 0.3;
                            hold on;
                            %%blue vertical lines
                            y1=get(gca,'ylim');
                            x1=bluepeak_start*srF;
                            hold on;
                            p3=plot([x1 x1],y1,'--','Color','b');
                            p3.Color(4) = 0.3;
                            hold on;
                            y1=get(gca,'ylim');
                            x1=bluepeak_end  *srF;
                            hold on;
                            p4=plot([x1 x1],y1,'--','Color','b');
                            p4.Color(4) = 0.3;
                            xlim([1,10000])
                        end
                        
                        if gof.adjrsquare>0.9
                            pulse2_BL = diff_bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:);
                            pulse2_RW = diff_bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:);
                            fitting_used = 1;
                        else
                            pulse2_BL = bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:);
                            pulse2_RW = bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:);
                        end
                    catch
                        pulse2_BL = bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:);
                        pulse2_RW = bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:);
                    end
                end
            end
            
            [r, c]=find(map_pattern_all(:,:,j)==i);
            if vclamp(j)==1 && double_pulse(j)==1 % AMPA double pulse
                if i ==1
                    AMPArep = AMPArep + 1;
                end
                pulseBL_and_RW{1,1,1,AMPArep}(r,c,:) = pulse1_BL; % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW{2,1,1,AMPArep}(r,c,:) = pulse1_RW;
                pulseBL_and_RW{1,2,1,AMPArep}(r,c,:) = pulse2_BL; % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW{2,2,1,AMPArep}(r,c,:) = pulse2_RW;
                
                pulseBL_and_RW_m{1,1,1,AMPArep}(r,c,:) = mean(pulse1_BL); % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW_m{2,1,1,AMPArep}(r,c,:) = mean(pulse1_RW);
                pulseBL_and_RW_m{1,2,1,AMPArep}(r,c,:) = mean(pulse2_BL); % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW_m{2,2,1,AMPArep}(r,c,:) = mean(pulse2_RW);
                
            elseif vclamp(j)==0 && double_pulse(j)==1 % NMDA double pulse
                if i ==1
                    NMDArep = NMDArep + 1;
                end
                pulseBL_and_RW{1,1,2,NMDArep}(r,c,:) = pulse1_BL; % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW{2,1,2,NMDArep}(r,c,:) = pulse1_RW;
                pulseBL_and_RW{1,2,2,NMDArep}(r,c,:) = pulse2_BL; % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW{2,2,2,NMDArep}(r,c,:) = pulse2_RW;
                
                pulseBL_and_RW_m{1,1,2,NMDArep}(r,c,:) = mean(pulse1_BL); % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW_m{2,1,2,NMDArep}(r,c,:) = mean(pulse1_RW);
                pulseBL_and_RW_m{1,2,2,NMDArep}(r,c,:) = mean(pulse2_BL); % pulseBL_and_RW dim: BL/RW, pulse1/2,AMPA/NMDA, rep,
                pulseBL_and_RW_m{2,2,2,NMDArep}(r,c,:) = mean(pulse2_RW);
                
            end
            
            %ephys_traces
            ephys_traces(:,counter,j)=bs_traces(1:srF*1000);
            [r, c]=find(map_pattern_all(:,:,j)==i);
            ephys_traces_grid(:,r,c,j)=bs_traces(1:srF*1000); 
            
            if fitting_used % exist('gof','var')==1
                gof_fit(counter,j)=gof.adjrsquare; clear gof.adjrsquare
                fit_param{counter,j} = f;
                fit_traces(:,counter,j)=yf(1:srF*1000);
            else
                gof_fit(counter,j) = NaN;
                fit_param{counter,j} = NaN;
                fit_traces(:,counter,j)=nan(1,srF*1000);
            end
            
            counter=counter+1;
            
            
        end
        if 0% map_trace==1
            figure;
            plot(ephys_traces(:,:,j),'k'); hold on; xlim([1 500*srF])
            plot(mean(ephys_traces(:,:,j),2),'r'); hold on;
            if vclamp(j)==1
                plot([BL_AMPA_P1(1)*srF,BL_AMPA_P1(1)*srF],[ylim],'--k');
                plot([BL_AMPA_P1(2)*srF,BL_AMPA_P1(2)*srF],[ylim],'--k');
                
                plot([RW_AMPA_P1(1)*srF,RW_AMPA_P1(1)*srF],[ylim],'--r');
                plot([RW_AMPA_P1(2)*srF,RW_AMPA_P1(2)*srF],[ylim],'--r');
                
                plot([BL_AMPA_P2(1)*srF,BL_AMPA_P2(1)*srF],[ylim],'--k');
                plot([BL_AMPA_P2(2)*srF,BL_AMPA_P2(2)*srF],[ylim],'--k');
                
                plot([RW_AMPA_P2(1)*srF,RW_AMPA_P2(1)*srF],[ylim],'--b');
                plot([RW_AMPA_P2(2)*srF,RW_AMPA_P2(2)*srF],[ylim],'--b');
            else
                plot([BL_NMDA_P1(1)*srF,BL_NMDA_P1(1)*srF],[ylim],'--k');
                plot([BL_NMDA_P1(2)*srF,BL_NMDA_P1(2)*srF],[ylim],'--k');
                
                plot([RW_NMDA_P1(1)*srF,RW_NMDA_P1(1)*srF],[ylim],'--r');
                plot([RW_NMDA_P1(2)*srF,RW_NMDA_P1(2)*srF],[ylim],'--r');
                
                plot([BL_NMDA_P2(1)*srF,BL_NMDA_P2(1)*srF],[ylim],'--k');
                plot([BL_NMDA_P2(2)*srF,BL_NMDA_P2(2)*srF],[ylim],'--k');
                
                plot([RW_NMDA_P2(1)*srF,RW_NMDA_P2(1)*srF],[ylim],'--b');
                plot([RW_NMDA_P2(2)*srF,RW_NMDA_P2(2)*srF],[ylim],'--b');
            end
            
        end
    end
    
    traces=[];
end

tempr = squeeze(pulseBL_and_RW(2,1,1,:));
tempb = squeeze(pulseBL_and_RW(2,2,1,:));
for i = 1:length(tempr)
    temp2(:,:,:,i)=squeeze(cat(3,tempr{i},tempb{i}));
end
Map_trace_m = mean(temp2,4);
    
if 0%map_trace
%     tempb = squeeze(pulseBL_and_RW(1,1,1,:));

    
    figure;
    Rect    = [0.1, 0.1, 0.8, 0.8];
    AxisPos = myPlotPos(size(Map_trace_m,1), size(Map_trace_m,2), Rect);
    for i = 1:size(Map_trace_m,1)*size(Map_trace_m,1)
        axes('Position', AxisPos(i, :))
        row = floor((i-1)/size(Map_trace_m,1))+1;
        col = rem((i-1),size(Map_trace_m,1))+1;
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
    set(gcf,'color','w');
    
end
try; autoArrangeFigures; end

pulseBL_and_RW_mean = cellfun(@(x) mean(x,3),pulseBL_and_RW,'UniformOutput', 0);
pulseBL_and_RW_peak(:,:,1,:) = cellfun(@(x) min(x,[],3),pulseBL_and_RW(:,:,1,:),'UniformOutput', 0);
if size(pulseBL_and_RW,3) == 2
    pulseBL_and_RW_peak(:,:,2,:) = cellfun(@(x) max(x,[],3),pulseBL_and_RW(:,:,2,:),'UniformOutput', 0);
end
pulseBL_and_RW_std = cellfun(@(x) std(x,[],3),pulseBL_and_RW,'UniformOutput', 0);

% pulseBL_and_RW_test_std = squeeze(cellfun(@(a,b) abs(a)>b.*fc,permute(pulseBL_and_RW_peak(2,:,:,:),[2,3,4,1]), permute(pulseBL_and_RW_std(1,:,:,:),[2,3,4,1]), 'UniformOutput', 0));
pulseBL_and_RW_test_std = cellfun(@(a,b) abs(a)>b.*fc,permute(pulseBL_and_RW_peak(2,:,:,:),[2,3,4,1]), permute(pulseBL_and_RW_std(1,:,:,:),[2,3,4,1]), 'UniformOutput', 0);

% pulseBL_and_RW_test_rs = squeeze(cellfun(@(a,b)  myfun(a,b)<(0.05/sqrt(size(a,1))),squeeze(pulseBL_and_RW(1,:,:,:)), squeeze(pulseBL_and_RW(2,:,:,:)), 'UniformOutput', 0));
%pulse1/2,AMPA/NMDA, rep

output_map = cellfun(@(a,b) a.*b,permute(pulseBL_and_RW_peak(2,:,:,:),[2,3,4,1]), pulseBL_and_RW_test_std, 'UniformOutput', 0);
% 
% red_scracm_AMPA = reshape([output_map{1,1,:}],8,8,[]);
% red_scracm_NMDA = reshape([output_map{1,2,:}],8,8,[]);
% blue_scracm_AMPA = reshape([output_map{2,1,:}],8,8,[]);
% blue_scracm_NMDA = reshape([output_map{2,2,:}],8,8,[]);

red_scracm_AMPA = reshape([output_map{1,1,:}],7,7,[]);
blue_scracm_AMPA = reshape([output_map{2,1,:}],7,7,[]);

if size(pulseBL_and_RW,3) == 2
    red_scracm_NMDA = reshape([output_map{1,2,:}],7,7,[]);
    blue_scracm_NMDA = reshape([output_map{2,2,:}],7,7,[]);
end

red_scracm_AMPA(red_scracm_AMPA==0) = nan;
blue_scracm_AMPA(blue_scracm_AMPA==0) = nan;
if size(pulseBL_and_RW,3) == 2
    red_scracm_NMDA(red_scracm_NMDA==0) = nan;
    blue_scracm_NMDA(blue_scracm_NMDA==0) = nan;
end

red_scracm_AMPA_mean = nanmean(red_scracm_AMPA,3);
blue_scracm_AMPA_mean = nanmean(blue_scracm_AMPA,3);
if size(pulseBL_and_RW,3) == 2
    red_scracm_NMDA_mean = nanmean(red_scracm_NMDA,3);
    blue_scracm_NMDA_mean = nanmean(blue_scracm_NMDA,3);
end
% calc estimated sample size required to detect events of 10pA
% mean(mean([pulseBL_and_RW_mean{1,:,1,:}]))
% mean(mean([pulseBL_and_RW_mean{2,:,1,:}]))
% nout = sampsizepwr('t2',[0 mean(mean([pulseBL_and_RW_std{1,:,1,:}]))],-5,0.95);

%
all_traces_AMPA = reshape(ephys_traces(:,:,find([vclamp==1 & double_pulse==1])),size(ephys_traces,1),[]); % AMPA
if size(pulseBL_and_RW,3) == 2
    all_traces_NMDA = reshape(ephys_traces(:,:,find([vclamp==0 & double_pulse==1])),size(ephys_traces,1),[]); % NMDA
end
% figure; plot(all_traces_AMPA,'k');

PD1_r_AMPA = PD1_r(:,:,find([vclamp==1 & double_pulse==1]));
PD2_b_AMPA = PD2_b(:,:,find([vclamp==1 & double_pulse==1]));
if size(pulseBL_and_RW,3) == 2
    PD1_r_NMDA = PD1_r(:,:,find([vclamp==0 & double_pulse==1]));
    PD2_b_NMDA = PD2_b(:,:,find([vclamp==0 & double_pulse==1]));
end
% same for irradiance (calc is wrong so disabled)
% IR1_r_AMPA = IR1_r(:,:,find([vclamp==1 & double_pulse==1]));
% IR1_r_NMDA = IR1_r(:,:,find([vclamp==0 & double_pulse==1]));
% IR2_b_AMPA = IR2_b(:,:,find([vclamp==1 & double_pulse==1]));
% IR2_b_NMDA = IR2_b(:,:,find([vclamp==0 & double_pulse==1]));

% organize output vars
resp_maps.red_scracm_AMPA_mean = red_scracm_AMPA_mean;
resp_maps.blue_scracm_AMPA_mean = blue_scracm_AMPA_mean;
if size(pulseBL_and_RW,3) == 2
    resp_maps.red_scracm_NMDA_mean = red_scracm_NMDA_mean;
    resp_maps.blue_scracm_NMDA_mean = blue_scracm_NMDA_mean;
else
    resp_maps.red_scracm_NMDA_mean = [];
    resp_maps.blue_scracm_NMDA_mean = [];
end

traces.all_traces_AMPA = all_traces_AMPA;
traces.ampa_map_order = map_pattern_all(:,:,find([vclamp==1 & double_pulse==1]));
traces.all_traces_AMPA_grid = ephys_traces_grid(:,:,:,find([vclamp==1 & double_pulse==1]));
traces.Map_trace_m = Map_trace_m;
if size(pulseBL_and_RW,3) == 2
    traces.all_traces_NMDA = all_traces_NMDA;
    traces.nmda_map_order = map_pattern_all(:,:,find([vclamp==0 & double_pulse==1]));
    traces.all_traces_NMDA_grid = ephys_traces_grid(:,:,:,find([vclamp==0 & double_pulse==1]));
else
    traces.all_traces_NMDA = [];
    traces.nmda_map_order = [];
    traces.all_traces_NMDA_grid = [];
end

illumination.PD1_r_AMPA = PD1_r_AMPA;
illumination.PD2_b_AMPA = PD2_b_AMPA;
if size(pulseBL_and_RW,3) == 2
    illumination.PD1_r_NMDA = PD1_r_NMDA;
    illumination.PD2_b_NMDA = PD2_b_NMDA;
else
    illumination.PD1_r_NMDA = [];
    illumination.PD2_b_NMDA = [];
end

resp_windows.BL_AMPA_P1 = BL_AMPA_P1;
resp_windows.BL_AMPA_P2 = BL_AMPA_P2;
if size(pulseBL_and_RW,3) == 2
    resp_windows.BL_NMDA_P1 = BL_NMDA_P1;
    resp_windows.BL_NMDA_P2 = BL_NMDA_P2;
else
    resp_windows.BL_NMDA_P1 = [];
    resp_windows.BL_NMDA_P2 = [];
end

resp_windows.RW_AMPA_P1 = RW_AMPA_P1;
resp_windows.RW_AMPA_P2 = RW_AMPA_P2;
if size(pulseBL_and_RW,3) == 2
    resp_windows.RW_NMDA_P1 = RW_NMDA_P1;
    resp_windows.RW_NMDA_P2 = RW_NMDA_P2;
else
    resp_windows.RW_NMDA_P1 = [];
    resp_windows.RW_NMDA_P2 = [];
end

% outputs
% resp_maps
% traces
% irradiance



%%
% [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = ...
%     extract_events(-1*ephys_traces(:,:,idx70), [BL_AMPA_P2(1)*srF:BL_AMPA_P2(2)*srF], [RW_AMPA_P2(1)*srF:RW_AMPA_P2(2)*srF], [], 1);


end

function out=myfun(x,y)
for g=1:size(x,1)
    for gg = 1:size(x,2)
        out(g,gg) = ranksum(squeeze(x(g,gg,:)),squeeze(y(g,gg,:)));
    end
end
end

function plot_stim_resp_map(ind_traces,baseline,mapPattern,redpeak_start,bluepeak_start,bluepeak_end,srF)

flipFlag = 0;
bs_traces=ind_traces-mean(baseline);%subtract baseline
array=bs_traces;
traceLength = size(array,1 );
newArray = [];
for n=1:numel(mapPattern)
    newArray(:,find(mapPattern==n)) = array(:,n);
end
startTime = 1;
stopTime = traceLength;
showStart = (redpeak_start-100)*srF;
showStop = (bluepeak_end+200)*srF;
array = newArray(startTime:stopTime, :);
[rows,cols] = size(array);
totalTime = (rows-1)/(srF*1000);
xTimeAxis = linspace(0, totalTime, rows)';
traceAxis = ( 1 : cols );
[sizeX, sizeY] = size(mapPattern);
yFactor = 200; % offset, in pA
scaleBarTxt = 'pA';
if flipFlag == 1
    yFactor = -yFactor;
end

offsetVector = yFactor * ( 0 : cols-1 );
offsetArray = repmat(offsetVector, rows, 1);
array = array-offsetArray;

% set up the figure -------------------------------------------------------------
x = .14;
y = .14;
w = .8;
h = .8;
F=figure;
subplotRows = 1; subplotCols = sizeY; plotnum = 0;
for N = 1:sizeY % plot traces
    startInd = (N-1)*sizeX + 1;
    endInd = N*sizeX;
    plotnum = plotnum+1;
    
    %     hsub(plotnum) = subplot(subplotRows,subplotCols,plotnum);
    
    %     pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    %      pos2 = 0.02;
    %     pos3 = 0.05;
    %     pos4 = 0.96;
    
    pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    pos2 = 0.02;
    pos3 = 0.038;
    pos4 = 0.96;
    
    hsub(N) = axes('Position', [pos1 pos2 pos3 pos4]);
    
    %     set(gca, 'Position', );
    
    plot(xTimeAxis(showStart:showStop), array(showStart:showStop,startInd:endInd),'Color','k');
    hold on;
    y1=get(gca,'ylim');
    x1= redpeak_start/(srF*100);
    hold on;
    p1=plot([x1 x1],y1,'--','Color','r');
    p1.Color(4) = 0.3;
    hold on;
    %                 y1=get(gca,'ylim');
    %                 x1=redpeak_end/(srF*100);
    %                 hold on;
    %                 p2=plot([x1 x1],y1,'--','Color','r');
    %                 p2.Color(4) = 0.3;
    %                 hold on;
    %blue vertical lines
    y1=get(gca,'ylim');
    x1=bluepeak_start/(srF*100);
    hold on;
    p3=plot([x1 x1],y1,'--','Color','b');
    p3.Color(4) = 0.3;
    hold on;
    y1=get(gca,'ylim');
    x1=bluepeak_end/(srF*100);
    hold on;
    p4=plot([x1 x1],y1,'--','Color','b');
    p4.Color(4) = 0.3;
    
    minval = min(mean(array(1:100,startInd:endInd)));
    maxval = max(mean(array(1:100,startInd:endInd)));
    tweakFactor = abs(maxval - minval)*0.05;
    yrange = [minval-tweakFactor maxval+tweakFactor];
    set(gca, 'YLim', yrange);
    set(gca, 'XLim', [(showStart-200)/(srF*1000) (showStop+200)/(srF*1000)]);
    xlabel('Seconds');
    set(gcf,'color','w');
    
    
end
set(findobj(gcf, 'Type', 'axes'), 'Visible', 'off');

% scalebar lines
Y = mean(array(:,end))+yFactor/4;
% Y = min(get(gca, 'YLim'));
hscalebar = line([.1 .2], [Y Y]);
set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');
hscalebar = line([.1 .1], [Y Y+yFactor/2]);
set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');

% scalebar text
ht(1) = text(.12, Y+yFactor/6, '100 ms');
ht(2) = text(.12, Y+yFactor/3, [num2str(yFactor/2) ' ' scaleBarTxt]);
set(ht, 'Color', 'k', 'FontSize', 8, 'Tag', 'scaleBarText');

end

function [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = extract_events(bs_traces, baseline_window, response_window, fc, show_plot)
if nargin<4
    show_plot = 0;
end

mean_resp = nanmean(bs_traces(response_window,:));
mean_base = nanmean(bs_traces(baseline_window,:));
bs_std = std(bs_traces(baseline_window,:),1);
resp_thresh = bs_std*4;

try
    test = signrank(mean_base,mean_resp)<0.05;
catch
    test = 0;
end

peaks = max(bs_traces(response_window,:));

if test
    if show_plot
        % To see how the distribution of peak hights differse between
        % baseline and stim period
        peak_base = min(bs_traces(baseline_window,:));
        figure; hist([peaks;peak_base]',100)
    end
    
    % error occurs here, check!
    resp_idx = peaks>mean(resp_thresh); % which trials show a response
    
    % Fit exp decay in order to only measure steady state responses
    try
        % single exponential
        sucs_peaks = peaks(find(resp_idx))';
        x = [0:sum(resp_idx)-1]';
        fo = fitoptions('Method','NonlinearLeastSquares', ...
            'Lower',[0 -1 0], ...
            'Upper', [max(sucs_peaks)*10 0 max(sucs_peaks)],...
            'StartPoint',[max(sucs_peaks)*3 -0.0005 min(sucs_peaks)]); % adjust start settings
        ft = fittype('a*exp(b*(x))+c', 'options', fo);
        [f gof]=fit(x,sucs_peaks,ft);
        yf=f.a.*exp(f.b.*x)+f.c;
        if show_plot
            figure; plot(x,yf); hold on; plot(x,sucs_peaks,'ro')
        end
        
        R2_exp = sum((sucs_peaks-yf).^2);
        R2_linear = sum((sucs_peaks-mean(sucs_peaks)).^2);
        
        if R2_exp<R2_linear
            % set respones before steady state to nans
            steady_state=find(yf<((yf(1)-f.c)/4+f.c),1); % define steady state at point where 3/4 of the dacay has taken place
            temp = find(resp_idx,steady_state);
            steady_state = temp(end); % needs to be adjusted to mean steady state based on all trials instead of all sucesses
            if isempty(steady_state)
                steady_state = 1;
            end
        else
            steady_state=1;
        end
        
        resp_prop = sum(resp_idx(steady_state:end))/length(resp_idx(steady_state:end));
        temp = peaks(steady_state:end);
        avg_amp = mean(temp(find(resp_idx(steady_state:end))));
        
    catch
        steady_state = nan;
        resp_prop = sum(resp_idx)/length(resp_idx);
        avg_amp = mean(peaks(find(resp_idx)));
    end
    
else
    resp_idx = zeros(size(bs_traces,2),1);
    steady_state = nan;
    resp_prop = 0;
    avg_amp = 0;
end

end

