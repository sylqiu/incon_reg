% run inconisistent shape registration for 50x50 meshes

prefix = 'HDM 2016_preprocessed/';

data = dir([prefix,'*.mat']);
data_list = ({data.name})';
% data_list = natsortfiles({data.name})';

%%
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 1; % turn higher if you want stronger intensity matching
param.demons_stepsize = 5;% turn higher if you want stronger intensity matching
% param.landmark_iter = 1;
% param.overall_iter = 20;
landmark_iter = 1;
overall_iter = 20;

Hdiff_all = zeros(length(data_list),length(data_list));
Kdiff_all = zeros(length(data_list),length(data_list));
Hndiff_all = zeros(length(data_list),length(data_list)); % difference of the normalized H
Kndiff_all = zeros(length(data_list),length(data_list)); % difference of the normalized K
mu_all = zeros(length(data_list),length(data_list));

tic;

for i = 45%1:length(data_list)
    %% load source data and preprocess
    disp(['i = ',num2str(i)]);
    
    source_mat =  load([prefix,data_list{i}]);
    source_face = source_mat.f;
    source_vertex = source_mat.v;
    source_gauss_curvature = source_mat.K;
    source_mean_curvature = source_mat.H;
    landmark_source_index = source_mat.lm;
    
    source_mean_curvature_normalized = intensity_normalization(source_mean_curvature);
    source_gauss_curvature_normalized = intensity_normalization(source_gauss_curvature);
    
    num_landmark = length(landmark_source_index);
    
    % build laplacian of the source domain for coefficient smoothing
    L = cotmatrix(source_vertex, source_face);
    F2Vm = F2V(source_vertex', source_face');
    V2Fm = V2F(source_vertex', source_face');
        
    % assume one boundary component                           
    [~, source_boundary_index] = meshboundaries(source_face);
    
    %% Target
    for j = 49%1:length(data_list)
        %% load target data and preprocess
        target_mat = load([prefix,data_list{j}]);

        target_face = target_mat.f;
        target_vertex = target_mat.v;
        target_gauss_curvature = target_mat.K;
        target_mean_curvature = target_mat.H;
        landmark_target_index = target_mat.lm;

        % choose matching intensity to be gauss curvature
        [target_intensity, min_, max_] = intensity_normalization(target_gauss_curvature);
        source_intensity = intensity_normalization(source_gauss_curvature, min_, max_);
    
        target_mean_curvature_normalized = intensity_normalization(target_mean_curvature);
        target_gauss_curvature_normalized = target_intensity;
        
%         disp('furtherst two points in the target will be used for conformal flattening \n');
        distance_table = repmat(reshape(target_vertex(landmark_source_index, :), ...
            [1, num_landmark, 3]), [num_landmark, 1, 1]) - ...
            repmat(reshape(target_vertex(landmark_source_index, :), ...
            [num_landmark, 1, 3]), [1, num_landmark, 1]);
        distance_table = sum(distance_table.^2, 3);
        [~, sorted_] = sort(distance_table(:));
        lscm_target_ind1 = mod(sorted_(end), num_landmark);
        lscm_target_ind2 = round(sorted_(end)/num_landmark + 0.5);
        lscm_target_ind = [landmark_target_index(lscm_target_ind1), landmark_target_index(lscm_target_ind2)];
        lscm_source_ind = [landmark_source_index(lscm_target_ind1), landmark_source_index(lscm_target_ind2)];

        % conformal flattening
        flat_source_vertex = lscm(source_vertex, source_face, lscm_source_ind, [0, 0; 1, 0]);
        flat_target_vertex = lscm(target_vertex, target_face, lscm_target_ind, [0, 0; 1, 0]);
%         disp('Conformal flattening done\n');
        landmark_target_pos = flat_target_vertex(landmark_target_index,1:2);
        % assume one boundary component                           
        [~, target_boundary_index] = meshboundaries(target_face);

        %% Inconsistent registration
        initial_landmark_error = compute_landamrk_err(flat_source_vertex, landmark_source_index,...
            landmark_target_pos);
        [~,~,source_intensity_grid, target_intensity_grid, ~,~] = combine_to_same_grid(...
                                                flat_source_vertex, source_intensity,...
                                                flat_target_vertex, target_intensity,...
                                                source_boundary_index, target_boundary_index);
        initial_intensity_error = sum(abs(source_intensity_grid(:) - target_intensity_grid(:)));

        % algo begins
        source_vertex_reg_pre = flat_source_vertex;
        landmark_err = [initial_landmark_error];
        intensity_err = [initial_intensity_error];
        fprintf('L1 norm intensity difference % f ', initial_intensity_error);
        fprintf('Eulidean landmark difference % f \n', initial_landmark_error);
        for iter = 1:overall_iter
            fprintf('Iter %d \n', iter);
            for l_iter = 1:landmark_iter
                source_vertex_reg_pre = reg_landmark(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                            landmark_source_index, landmark_target_pos,...
                                            source_boundary_index, F2Vm, V2Fm, L, param);
%                 figure(21); gpp_plot_mesh(source_face, source_vertex_reg_pre); title('Landmark matching step');
%                 drawnow;
            end
            [source_vertex_reg, ie] = reg_intensity(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                        source_intensity, flat_target_vertex, target_intensity,...
                                        source_boundary_index, target_boundary_index,...
                                        F2Vm, V2Fm, L, param);
%             figure(31); gpp_plot_mesh(source_face, source_vertex_reg); title('Intensity matching step');
%             drawnow;
            % calculate the landmark and intensity difference
            le = compute_landamrk_err(source_vertex_reg, ...
                                landmark_source_index,...
                                landmark_target_pos);
            fprintf('L1 norm intensity difference % f ', ie);
            fprintf('Eulidean landmark difference % f \n', le);
            landmark_err = [landmark_err, le];
            intensity_err = [intensity_err, ie];
            %
            source_vertex_reg_pre = source_vertex_reg;

        end
        % algo ends
        %% Get the results
        [source_vertex_reg_intersect_index, correspondence_mask, ...
            source_vertex_reg_3D, source_face_reg, displace, dist,...
            target_intensity_reg, intensity_diff] = ...
            prepare_result(source_face, source_vertex, source_vertex_reg,...
                            target_face, target_vertex, flat_target_vertex,...
                            source_intensity, target_intensity);
        %                
        % interpolate target curvature to registered
        target_mean_curvature_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_mean_curvature, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_mean_curvature_reg(isnan(target_mean_curvature_reg)) = 0;
        target_gauss_curvature_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_gauss_curvature, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_gauss_curvature_reg(isnan(target_gauss_curvature_reg)) = 0;
        
        % curvature difference over the common domain
        Hdiff = (abs(source_mean_curvature(source_vertex_reg_intersect_index) - target_mean_curvature_reg(source_vertex_reg_intersect_index)));
        Kdiff = (abs(source_gauss_curvature(source_vertex_reg_intersect_index) - target_gauss_curvature_reg(source_vertex_reg_intersect_index)));
        
        % interpolate target curvature with normalization to registered
        target_mean_curvature_normalized_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_mean_curvature_normalized, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_mean_curvature_normalized_reg(isnan(target_mean_curvature_normalized_reg)) = 0;
        target_gauss_curvature_normalized_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_gauss_curvature_normalized, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_gauss_curvature_normalized_reg(isnan(target_gauss_curvature_normalized_reg)) = 0;
        
        % normalized curvature difference over the common domain (use this or the version rescaled with min_, max_?)
        Hndiff = (abs(source_mean_curvature_normalized(source_vertex_reg_intersect_index) - target_mean_curvature_normalized_reg(source_vertex_reg_intersect_index)));
        Kndiff = (abs(source_gauss_curvature_normalized(source_vertex_reg_intersect_index) - target_gauss_curvature_normalized_reg(source_vertex_reg_intersect_index)));
        
        % BC over the common domain
%         mu = F2Vm*abs(compute_bc(source_face, flat_source_vertex, source_vertex_reg, 2));   
%         mu = (mu(source_vertex_reg_intersect_index,:));
        % better to be face-based?
        mu = abs(compute_bc(source_face, flat_source_vertex, source_vertex_reg, 2));   
        mu = mu(source_face_reg,:);
        
        % scalar quantities capturing the shape difference
        % use mean instead of sum because the size of the overlapping region is different for different pairs of meshes?
        Hdiff_all(i,j) = mean(Hdiff); 
        Kdiff_all(i,j) = mean(Kdiff);
        Hndiff_all(i,j) = mean(Hndiff);
        Kndiff_all(i,j) = mean(Kndiff);
        mu_all(i,j) = mean(mu); 
                        
        %%
% %         result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
% %                             target_face, target_vertex, flat_target_vertex,...
% %                             source_intensity, target_intensity,...
% %                             landmark_source_index, landmark_target_index,...
% %                             source_vertex_reg_intersect_index, correspondence_mask, ...
% %                             source_vertex_reg_3D, source_face_reg, displace, dist,...
% %                             target_intensity_reg, intensity_diff,...
% %                             landmark_err, intensity_err);
                        
    end
end

toc;