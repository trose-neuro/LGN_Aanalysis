% Add localization to Data structure

%% load data strucutre
% data_loc = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';
% cd(data_loc);
% 
% load('Data_MF_10-19-Jul-2019','LGN'); data_MF = LGN; clear LGN
% load('Data_SW_11-19-Jul-2019','LGN'); data_SW = LGN; clear LGN
% data = cat(1,data_SW,data_MF); clear data_MF data_SW

%% Load most recent structures
% data_loc = 'I:\Martin Fernholz\LGN _Project_Common\AnalyzedData';
data_loc = 'D:\LGN project';

cd(data_loc);
d = dir('Data_SW*.mat');
load(d(find([d(:).datenum]==max([d(:).datenum]))).name,'LGN'); data_SW = LGN; clear LGN
d = dir('Data_MF*.mat');
load(d(find([d(:).datenum]==max([d(:).datenum]))).name,'LGN'); data_MF = LGN; clear LGN

data = cat(1,data_SW,data_MF); clear data_MF data_SW

% data = data_SW;
skip_morph = 1;

%% get schematic stack

schematic_file_loc = 'I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\';
schematic_file_name = 'LGN_schematic_stack.tif';
cd(schematic_file_loc)

LGN_schematic = load_tiff('LGN_schematic_stack.tif');
xlimit = size(LGN_schematic,2);
ylimit = size(LGN_schematic,1);
zlimit = size(LGN_schematic,3)/2;
LGN_schematic = logical(LGN_schematic);
LGN_schematic_4dstack(:,:,:,1) = LGN_schematic(:,:,1:2:end-1); % ipsi
LGN_schematic_4dstack(:,:,:,2) = LGN_schematic(:,:,2:2:end); % contra

LGN_schamtic_zproject(:,:,1) = max(LGN_schematic_4dstack(:,:,:,1),[],3); % ipsi
LGN_schamtic_zproject(:,:,2) = max(LGN_schematic_4dstack(:,:,:,2),[],3); % contra

[ipsi_perim(:,1), ipsi_perim(:,2)] = find(bwperim(LGN_schamtic_zproject(:,:,1)));
[contra_perim(:,1), contra_perim(:,2)]= find(bwperim(LGN_schamtic_zproject(:,:,2)));

% loc in schematic
remove_cell_idx = []; n = 0;

for i = 1:length(data)
    if strcmp(data(i).experimentator,'SW')
        experimenter_name_full = 'Simon';
    elseif strcmp(data(i).experimentator,'MF')
        experimenter_name_full = 'Martin';
    end
    
    spacial_data_loc = ['I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\' data(i).animal_name ...
        '\ephys ' experimenter_name_full '\slice ' num2str(data(i).slice_nr) '\Schematic\'];
    data(i).schematic_file_loc  = [schematic_file_loc schematic_file_name];
    
    % load schematic location position
    try
        cd(spacial_data_loc)
        try 
            temp = ReadImageJROI('cell_loc_schematic.zip'); 
        catch
            temp = ReadImageJROI('cell_loc_schematic.roi');
        end
        temp = [temp{:}];
        find_cell = find(strcmp({temp.strName},[data(i).experimentator data(i).cellname]));
        data(i).schematic_loc_xyz = [temp(find_cell).vnRectBounds(2:-1:1)  temp(find_cell).vnPosition(2)];
        if isempty(data(i).schematic_loc_xyz)
            abcd(1000); %cause error
        end
        data(i).schematic_imagej_zipfile_loc = [spacial_data_loc 'cell_loc_schematic.zip'];
    catch
        disp(['Mouse ' data(i).animal_name ' Cell ' [data(i).patching_date(1:6) data(i).experimentator data(i).cellname] ' location in schematic not found'])
        n = n+1;
        remove_cell_idx(n) = i;
        data(i).schematic_loc_xyz = [];
        data(i).schematic_imagej_zipfile_loc = [];
    end
    
    
    
    % schematic zone
    if ~isempty(data(i).schematic_loc_xyz)
        xyz = data(i).schematic_loc_xyz;
        cell_isin_contraZ = LGN_schematic_4dstack(xyz(2),xyz(1),xyz(3),2); % contra
        cell_isin_ipsiZ = LGN_schematic_4dstack(xyz(2),xyz(1),xyz(3),1); % ipsi
        if ~cell_isin_contraZ && ~cell_isin_ipsiZ
            data(i).schematic_loc_zone = 0;
        elseif cell_isin_contraZ && ~cell_isin_ipsiZ
            data(i).schematic_loc_zone = 1;
        elseif ~cell_isin_contraZ && cell_isin_ipsiZ
            data(i).schematic_loc_zone = 2;
        elseif cell_isin_contraZ && cell_isin_ipsiZ
            data(i).schematic_loc_zone = 3;
        end
    else
        data(i).schematic_loc_zone = [];
    end
    
    % get 2p overview stack
    overview_data_loc = ['I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\' data(i).animal_name ...
        '\ephys ' experimenter_name_full '\slice ' num2str(data(i).slice_nr) '\2P_Overview\'];
    overview_filename = '*Overview*';
    folder_contents = dir(overview_data_loc);
    
    n = 0;
    for ii = 3: length(folder_contents)
       if ~folder_contents(ii).isdir && contains(folder_contents(ii).name,'Overview') && contains(folder_contents(ii).name,'tif')
           n = n+1;
           overview_filenames{n} = folder_contents(ii).name;
       end
    end
    
    % for now we just use the first Overview stack found
    try
        overview_filenames = flip(sort(overview_filenames));
        clear overview_stack
        if length(imfinfo([overview_data_loc overview_filenames{1}]))>1
            
            overview_stack = load_tiff([overview_data_loc overview_filenames{1}]);
        else
            
            overview_stack = imread([overview_data_loc overview_filenames{1}]);
        end
        overview_stack = double(overview_stack);
        overview_stack(:,:,1) = (overview_stack(:,:,1)-min(min(overview_stack(:,:,1)))+1); %./(max(max(overview_stack(:,:,1)))-min(min(overview_stack(:,:,1))));
        overview_stack(:,:,2) = (overview_stack(:,:,2)-min(min(overview_stack(:,:,2)))+1); %./(max(max(overview_stack(:,:,2)))-min(min(overview_stack(:,:,2))));
        overview_stack(:,:,3) = zeros(size(overview_stack,1),size(overview_stack,2));

        clear overview_stack_adj
        temp = double(overview_stack(:,:,1));
        temp = temp./max(temp(:));
        overview_stack_adj(:,:,1) = imadjust(temp);
        temp = double(overview_stack(:,:,2));
        temp = temp./max(temp(:));
        overview_stack_adj(:,:,2) = imadjust(temp);
        
        overview_stack_adj(:,:,3)=zeros(size(overview_stack,1),size(overview_stack,2));
        %     imshow(overview_stack_adj(:,:,:))
        
        data(i).Overview2p = imresize(overview_stack_adj,0.5);
    catch
        data(i).Overview2p = [];
    end
    
    % confocal stack localization and FlogR
    try              
        FlourLogRatio_data_loc = ['I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\' data(i).animal_name ...
            '\ephys ' experimenter_name_full '\slice ' num2str(data(i).slice_nr) '\Cleared\'];
        FlourLogRatio_filename = [data(i).patching_date(1:6) data(i).hemisphere];
        
        temp = load([FlourLogRatio_data_loc FlourLogRatio_filename '_Cellinfo.mat']);
        current_info = [];
        for ii = 1:size(temp.Cell_info,2)
            if strcmp(temp.Cell_info{ii}.cell_name,[data(i).experimentator data(i).cellname])
                current_info = temp.Cell_info{ii};
            end
        end
        
        if ~isempty(current_info) && (size(current_info,1)==1 && size(current_info,2)==1)
            data(i).FlourLogRatio_mean = current_info.LogRatioF;
            data(i).FlourLogRatio_data = current_info;
            data(i).Confocal_data_loc = {[FlourLogRatio_data_loc FlourLogRatio_filename '.tif'] [FlourLogRatio_data_loc FlourLogRatio_filename '.zip'] [FlourLogRatio_data_loc FlourLogRatio_filename '_Cellinfo.mat']};
        else
            data(i).FlourLogRatio_mean = nan;
            data(i).FlourLogRatio_data = [];
            data(i).Confocal_data_loc = [];
        end
        
    catch
        disp(['Cell axon flourescence file for cell ' [data(i).patching_date(1:6) data(i).experimentator data(i).cellname] ' not found'])
        data(i).FlourLogRatio_mean = nan;
        data(i).FlourLogRatio_data = [];
        data(i).Confocal_data_loc = [];
    end
end

idx_xzymissing = cellfun(@isempty , {data.schematic_loc_xyz});
{data(idx_xzymissing).animal_name; data(idx_xzymissing).patching_date; data(idx_xzymissing).experimentator; data(idx_xzymissing).slice_nr; data(idx_xzymissing).cellname}'

%% add alternative measure to DOI
if ~skip_morph
    
    for i = 1:size(data,1)
        
        if ~isempty(data(i).morphology)
            dendrite_tree = data(i).morphology.traces;
            cell_center = [mean(data(i).morphology.somaX), mean(data(i).morphology.somaY), mean(data(i).morphology.somaZ)];
            
            % first perform 2D trace interpolation
            
            % define the target interpolations density [pt/px]
            target_density = 5;
            
            % interpolate in 2d
            [node_1, node_2] = find(dendrite_tree.dA);
            node_number = length(node_1);
            
            % start a counter for  each of the nodes
            node_counter = zeros(node_number+1,1);
            distance = zeros(node_number,1);
            
            % allocate memory for the interpolated segments
            clear cell
            segment_cell = cell(node_number, 1);
            
            %     figure
            for nodes = 1:node_number
                % count the start
                node_counter(node_1(nodes)) = node_counter(node_1(nodes)) + 1;
                
                % get the node coordinates for start and end
                node_start = [dendrite_tree.X(node_1(nodes)), ...
                    dendrite_tree.Y(node_1(nodes)),...
                    dendrite_tree.Z(node_1(nodes))];
                node_end = [dendrite_tree.X(node_2(nodes)), ...
                    dendrite_tree.Y(node_2(nodes)),...
                    dendrite_tree.Z(node_2(nodes))];
                % get the distance between the start and end
                node_distance = norm(node_end - node_start);
                distance(nodes) = node_distance;
                % if the counter value is > 2 or the distance is too long, skip
                if node_counter(node_1(nodes)) > 1 || node_distance > 2
                    continue
                end
                
                % calculate the number of points for the segment
                interp_points = round(node_distance*target_density);
                % assemble the segment
                segment_cell{nodes} = [linspace(node_start(1), node_end(1), interp_points);...
                    linspace(node_start(2), node_end(2), interp_points);...
                    linspace(node_start(3), node_end(3), interp_points)];
                
                %         plot(segment_cell{nodes}(1,:),segment_cell{nodes}(2,:),'o')
                %         hold('on')
            end
            
            dendrite_tree_interp{i} = cat(2,segment_cell{:})';
            
            
            %
            
            %coordinates in x and y
            XYd=[XP YP]';
            
            for degr=1:360;
                %Rotation matrix
                R=[cos(degr) -sin(degr); sin(degr) cos(degr)];
                
                %New coordinates
                new_XYd=R*XYd;
                
                %counting in the 4 quadrants
                a1_counts= sum((new_XYd(1,:) >= a1x(end)) & (new_XYd(1,:) <= a1x(1)) & (new_XYd(2,:) >= a1y(1)) & (new_XYd(2,:) <= a1y(end)));
                a2_counts= sum((new_XYd(1,:) >= a2x(1)) & (new_XYd(1,:) <= a2x(end)) & (new_XYd(2,:) >= a2y(end)) & (new_XYd(2,:) <= a2y(1)));
                b1_counts= sum((new_XYd(1,:) >= b1x(1)) & (new_XYd(1,:) <= b1x(end)) & (new_XYd(2,:) >= b1y(1)) & (new_XYd(2,:) <= b1y(end)));
                b2_counts= sum((new_XYd(1,:) >= b2x(end)) & (new_XYd(1,:) <= b2x(1)) & (new_XYd(2,:) >= b2y(end)) & (new_XYd(2,:) <= b2y(1)));
                
                %sum the two planes
                com_a=a1_counts+ a2_counts;
                com_b=b1_counts+ b2_counts;
                
                %keep all itirratively
                plane_a(:,degr)=com_a;
                plane_b(:,degr)=com_b;
                
                %Counting intersections
                delta_counts(:,degr)=abs(com_a-com_b);
            end
            
            %find maximum counts
            max_delta=find(delta_counts==max(delta_counts));
            
            %use maximum plane a and b
            com_planes=[plane_a(:,max_delta(1)) plane_b(:,max_delta(1))];
            
            for i = 1:length(max_delta)
                com_planes(i,:)=[plane_a(:,max_delta(i)) plane_b(:,max_delta(i))];
                DOi(i)=min(com_planes(i,:))/max(com_planes(i,:));
            end
            %
            
            
            
            try
                %%% rotated based in PCA (in 2d)
                [coeff,~,~,~,explained] = pca(dendrite_tree_interp{i}(:,1:2));
                dendrite_tree_interp_rotPC{i} = [(coeff*dendrite_tree_interp{i}(:,1:2)')' dendrite_tree_interp{i}(:,3) ]; % rotate only along x and y
                ratio(i) = (explained(2))/(explained(1)); % explained variance ratio between PCs
                
                
                %%% rotated based on center of mass (asimetry)
                centermass = mean(dendrite_tree_interp{i}(:,1:2));
                centermass_angle = -atan(centermass(1)/centermass(2))-pi/2;
                rot_matrix = [cos(centermass_angle) sin(centermass_angle); -sin(centermass_angle) cos(centermass_angle)];
                xyrot = [rot_matrix*dendrite_tree_interp{i}(:,1:2)']';
                centermass_new = rot_matrix*centermass';
                dendrite_tree_interp_rotasim{i} = [xyrot dendrite_tree_interp{i}(:,3) ]; % rotate only along x and y
                
                %             figure; hold on; plot([centermass(1)],[centermass(2)],'bx'); plot([0 centermass(1)],[0 centermass(2)],'b');
                %             plot([0 centermass_new(1)],[0 centermass_new(2)],'r'); plot([centermass_new(1)],[centermass_new(2)],'rx');
                %
                %             figure; scatter(dendrite_tree_interp{i}(:,1),dendrite_tree_interp{i}(:,2),5); hold on
                %             scatter(xyrot(:,1),xyrot(:,2),5,'r')
                
                for k = 1:length(centermass)
                    
                    if ~isempty(dendrite_tree_interp_rotasim{i}(dendrite_tree_interp_rotasim{i}(:,k)>0))
                        side1 = abs(mean(dendrite_tree_interp_rotasim{i}(dendrite_tree_interp_rotasim{i}(:,k)>0,k))); % number fo points
                        %                     side1 = abs(mean(dendrite_tree_interp_rot{i}(dendrite_tree_interp_rot{i}(:,k)>0))); % average point positions
                    else
                        side1=0;
                    end
                    if ~isempty(dendrite_tree_interp_rotasim{i}(dendrite_tree_interp_rotasim{i}(:,k)<0))
                        side2 = abs(mean(dendrite_tree_interp_rotasim{i}(dendrite_tree_interp_rotasim{i}(:,k)<0,k)));
                        %                     side2 = abs(mean(dendrite_tree_interp_rot{i}(dendrite_tree_interp_rot{i}(:,k)<0)));
                    else
                        side2=0;
                    end
                    asim(k) = abs((side1 - side2)/(side1 + side2));
                    
                end
                
                output_param(i,:) = [ratio(i) asim(1)];
            catch
                ratio(i) = nan;
                output_param(i,:) = [nan nan];
            end
        else
            ratio(i) = nan;
            output_param(i,:) = [nan nan];
        end
        
    end
    
    figure;
    scatter(output_param(:,1),output_param(:,2),'b','filled'); hold on
    xlabel(['Elongation (varexp PC2/PC1)']); ylabel(['Asimetry'])
    xlim([0 1]);ylim([0 1]); hold off
    
    for cell = 1:size(output_param,1)
        if ~isempty(data(cell).morphology)
            data(cell).morphology.pca_based_elongation_measure =  output_param(cell,1);
            data(cell).morphology.pca_based_asymetry_measure =  output_param(cell,2);
        end
    end
end

%%

disp('Saveing data')
% cd('I:\Martin Fernholz\LGN _Project_Common\AnalyzedData')
cd('D:\LGN project')
delete('Full_Data.mat')
save('Full_Data','data','-v7.3');

