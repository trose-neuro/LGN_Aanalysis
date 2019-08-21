function AMPA_NMDA_stepprot_ODI_plot(data,scatter_color)

%% automatic vs human classification missmatch all cells

figure;

% region lines
patch([-1 1 1 -1],[-1 -1 1 1]                       ,'b','FaceColor','none', 'EdgeColor', 'b', 'LineStyle', '--'); hold on %  true bino
patch([1.4 1.6 1.6 1.4],[1.4 1.4 1.6 1.6]           ,'r','FaceColor','none', 'EdgeColor', 'r', 'LineStyle', '--') % contra only
patch([-1.4 -1.6 -1.6 -1.4],[-1.4 -1.4 -1.6 -1.6]   ,'g','FaceColor','none', 'EdgeColor', 'g', 'LineStyle', '--') % ipsi only
patch([1.4 1.6 1.6 1.4],[-1 -1 1 1]                 ,'m','FaceColor','none', 'EdgeColor', 'm', 'LineStyle', '--') % contra & silent ipsi
patch([-1.4 -1.6 -1.6 -1.4],[-1 -1 1 1]             ,'c','FaceColor','none', 'EdgeColor', 'c', 'LineStyle', '--') % ipsi & contra ipsi
patch([-1 1 1 -1],[-1.4 -1.4 -1.6 -1.6]             ,'k','FaceColor','none', 'LineStyle', '--') % ipsi & contra ipsi
patch([-1 1 1 -1],[1.4 1.4 1.6 1.6]                 ,'k','FaceColor','none', 'LineStyle', '--') % ipsi & contra ipsi
patch([-1.4 -1.6 -1.6 -1.4],[1.4 1.4 1.6 1.6]       ,'k','FaceColor','none', 'LineStyle', '--') % contra only
patch([1.4 1.6 1.6 1.4],[-1.4 -1.4 -1.6 -1.6]       ,'k','FaceColor','none', 'LineStyle', '--') % contra only

number_of_nans = 0;
number_non_resp = 0;
for cellidx = 1:length(data)
    % make temporary var
    if data(cellidx).ocular_category == 0 %%|| data(cellidx).ocular_category == 6
        data_temp(cellidx).ocular_category = nan;
        number_non_resp = number_non_resp +1;
    else
        data_temp(cellidx).ocular_category = data(cellidx).ocular_category;
    end
    
    % shift the monocular cells by 0.5 +/- a bit to more to clearly display them in the scatter plot
    if abs(data(cellidx).ODI_AMPA_step)<1
%         data_temp(cellidx).ODI_AMPA_step = data(cellidx).ODI_AMPA_step;
        data_temp(cellidx).ODI_AMPA_step = data(cellidx).ODI_AMPA_step_peak;
    elseif data(cellidx).ODI_AMPA_step == 1
        data_temp(cellidx).ODI_AMPA_step = 1.5 + (rand*0.2-0.1); %+(rand*0.2-0.1);
    elseif data(cellidx).ODI_AMPA_step == -1
        data_temp(cellidx).ODI_AMPA_step = -1.5 + (rand*0.2-0.1); %+(rand*0.2-0.1);
    elseif isnan(data(cellidx).ODI_AMPA_step)
        data_temp(cellidx).ODI_AMPA_step = nan;
    end
    
    if abs(data(cellidx).ODI_NMDA_step)<1
        data_temp(cellidx).ODI_NMDA_step = data(cellidx).ODI_NMDA_step;
%         data_temp(cellidx).ODI_NMDA_step = data(cellidx).ODI_NMDA_step_peak;
    elseif data(cellidx).ODI_NMDA_step == 1
        data_temp(cellidx).ODI_NMDA_step = 1.5 + (rand*0.2-0.1); %+(rand*0.2-0.1);
    elseif data(cellidx).ODI_NMDA_step == -1
        data_temp(cellidx).ODI_NMDA_step = -1.5 + (rand*0.2-0.1); %+(rand*0.2-0.1);
    elseif isnan(data(cellidx).ODI_NMDA_step)
        data_temp(cellidx).ODI_NMDA_step = nan;
    end
    
    if strcmp(scatter_color,'color_cat')
        % set color based on human annotation
        if data_temp(cellidx).ocular_category == 1
            data_temp(cellidx).cell_cat_col = 'r'; % RED: contra_AMPA_contra_NMDA
        elseif data_temp(cellidx).ocular_category == 2
            data_temp(cellidx).cell_cat_col = 'g'; % GREEN: ipsi_AMPA_ipsi_NMDA
        elseif data_temp(cellidx).ocular_category == 3
            data_temp(cellidx).cell_cat_col = 'b'; % Blue: bino_AMPA_bino_NMDA
        elseif data_temp(cellidx).ocular_category == 4
            data_temp(cellidx).cell_cat_col = 'm'; % Magenta: contra_AMPA_bino_NMDA_ipsi silent
        elseif data_temp(cellidx).ocular_category == 5
            data_temp(cellidx).cell_cat_col = 'c'; % Cyan: ipsi_AMPA_bino_NMDA_contra silent
        elseif data_temp(cellidx).ocular_category == 6
            data_temp(cellidx).cell_cat_col = 'k'; % Black: responsive but strange
        end
    elseif strcmp(scatter_color,'Rs')
        Rs_step_AMPA = nanmean(data(cellidx).step_blue.neg_Rs(2,:));
        if ~isempty(data(cellidx).step_blue.pos_Rs)
            Rs_step_NMDA = nanmean(data(cellidx).step_blue.pos_Rs(2,:));
        else
            Rs_step_NMDA = nan;
        end
        Rs_change(cellidx) = (Rs_step_NMDA-Rs_step_AMPA)/(Rs_step_AMPA)*100;
        Rs_AMPA(cellidx) = Rs_step_AMPA;
        Rs_NMDA(cellidx) = Rs_step_NMDA;
%         data_temp(cellidx).cell_cat_col = log(Rs_change(cellidx));
        data_temp(cellidx).cell_cat_col = log(Rs_NMDA(cellidx));
    elseif strcmp(scatter_color,'MD')
        if data(cellidx).MD == 0
            data_temp(cellidx).cell_cat_col = 'b'; 
        elseif data(cellidx).MD == 1
            data_temp(cellidx).cell_cat_col = 'r'; 
        end
    end
    
end


if  strcmp(scatter_color,'color_cat') | strcmp(scatter_color,'MD')
    for cellidx = 1:length(data)
        % plot
        if data(cellidx).photodiode_flag | isempty(data_temp(cellidx).cell_cat_col)
            number_of_nans = number_of_nans+1;
        elseif ~isnan(data_temp(cellidx).ocular_category) && ~isnan(data_temp(cellidx).ODI_AMPA_step) ...
                && ~isnan(data_temp(cellidx).ODI_NMDA_step)
            scatter([data_temp(cellidx).ODI_AMPA_step],[data_temp(cellidx).ODI_NMDA_step],...
                20,data_temp(cellidx).cell_cat_col,'filled');
        elseif isnan(data_temp(cellidx).ocular_category) && ~isnan(data_temp(cellidx).ODI_AMPA_step) ...
                && ~isnan(data_temp(cellidx).ODI_NMDA_step)
            scatter([data_temp(cellidx).ODI_AMPA_step],[data_temp(cellidx).ODI_NMDA_step],20,'k');
        else
            number_of_nans = number_of_nans+1;
        end
    end
    
    if strcmp(scatter_color,'MD')
        legend({'no MD', 'MD'},'location','northeastoutside')
    end
    number_of_nans
    number_non_resp
    
elseif strcmp(scatter_color,'Rs')
    scatter([data_temp(:).ODI_AMPA_step],[data_temp(:).ODI_NMDA_step],...
        20,[data_temp(:).cell_cat_col],'filled')
    h = colorbar;
%     ylabel(h, {'log(Rs %change from AMPA',' to NMDA step prot)'})
    ylabel(h, {'log(Rs) during NMDA step prot.'})    
end

% make nice
ylabel('Step prot. NMDA ODI'); xlabel('Step prot. AMPA ODI');
xlim([-1.7 1.7]); ylim([-1.7 1.7]);
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
yticklabels({'-1' '-0.999' '-0.5' '0' '0.5' '0.999' '1'})
xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'-1' '-0.999' '-0.5' '0' '0.5' '0.999' '1'})
set(gca,'TickDir','out'); box on
set(gcf,'color','w');
title('Human vs Automatic ODI class');



clear data_temp
axis square;
