% Thesis Sanna Bakels, 4480279
% Step 3: Perform the mapping from the latent state to the source state

clear; close all; clc;

% Using the Matlab Toolbox EEGLAB from:
% Delorme, Arnaud, and Scott Makeig. "EEGLAB: an open source toolbox for 
% analysis of single-trial EEG dynamics including independent component 
% analysis." Journal of neuroscience methods 134.1 (2004): 9-21. 

% And the FieldTrip toolbox from:
% Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J. M. (2011). 
% FieldTrip: open source software for advanced analysis of MEG, EEG, 
% and invasive electrophysiological data. Computational intelligence 
% and neuroscience, 2011, 1-9.

% This code is partly based on the tutorial for source localization
% https://eeglab.org/tutorials/09_source/EEG_sources.html

%% Import data 
% Define the participant number
subject_number = 3;

% Initialize EEGLAB + load the file 
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename', ['EEG_subject' num2str(subject_number) '_dipfitReloc_125s.set'], 'filepath', ['/Users/sannabakels/Documents/Master/Afstuderen/Thesis/Matlab/ExperimentData/Subject ' num2str(subject_number) '/']);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% Addpath to ensure FieldTrip is running
addpath /~/Master/Afstuderen/Thesis/Matlab/Toolbox/fieldtrip-20230907
ft_defaults
savepath

%% Leadfield Matrix calculation
% Converting data from the EEGLAB format to fieldtrip
eegData = eeglab2fieldtrip(EEG,'preprocessing','none');
eegData.elec  = ft_convert_units(eegData.elec, 'mm');

% Loading the head model determined with EEGlab
vol = load('-mat', EEG.dipfit.hdmfile);
headmodel = vol.vol;

if strcmp(headmodel.type, 'bemcp')
    scalp_index = 3;
elseif strcmp(headmodel.type, 'dipoli')
    scalp_index = 1;
end

% Align the electrodes with the headshape
cfg = [];
cfg.method = 'project'; % onto scalp surface
cfg.headshape = headmodel.bnd(scalp_index); % scalp surface
eegData.elec = ft_electroderealign(cfg, eegData.elec);

%% Step 2: Calculate the foward solution
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
tlckEeg = ft_timelockanalysis(cfg, eegData);

cfg             = [];
cfg.headmodel   = headmodel;
cfg.resolution  = 20;
cfg.inwardshift = 5;
sourcemodel = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
cfg.elec        = tlckEeg.elec;

leadfield = ft_prepare_leadfield(cfg)

%% Calculate the leadfield matrices for further calculations in Granger Causality 
% Extract lead field matrices that are not empty
nonEmptyIndices = find(~cellfun('isempty', leadfield.leadfield));
nonEmptyPosData = cell(length(nonEmptyIndices), 1);

% Extract 'sourcemodel.pos' from non-empty cells and store in the cell array
for i = 1:length(nonEmptyIndices)
    cellIndex = nonEmptyIndices(i);
    nonEmptyPosData{i} = leadfield.leadfield{cellIndex};
end

numPosData = length(nonEmptyPosData);
posMatrix = zeros(numPosData, 3);
for j = 1:numPosData 
    posMatrix(j, :) = leadfield.pos(nonEmptyIndices(j), :);
end

% Calculate the lead field matrices for the magnitude of the vector
n_sources = length(nonEmptyPosData);
n_electrodes = length(nonEmptyPosData{1});
leadfield_mag = zeros(n_sources, n_electrodes);

for i = 1:n_sources
    matrix = nonEmptyPosData{i};

    for j = 1:n_electrodes
        cubeRootMag = matrix(j, 1)^2 + matrix(j,2)^2 + matrix(j,3)^2;
        % mag = cubeRootMag^(1/3);^2
        leadfield_mag(i, j) = sqrt(cubeRootMag);
    end
end


%% Interpolate the sourcemodel and the atlas
atlas = ft_read_atlas('/Users/sannabakels/Documents/Master/Afstuderen/Thesis/Matlab/Toolbox/fieldtrip-20230907/template/atlas/aal/ROI_MNI_V4.nii')

% and call ft_sourceinterpolate:
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
cfg.inwardshift = -0.5;
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, sourcemodel);

% Initialize a counter to keep track of the last assigned position
last_assigned = 0;

num_tissues = length(sourcemodel2.tissuelabel);
combined_sourcemodel.pos = zeros(length(sourcemodel.pos),3);
num_cells = 0;
leadfield2.pos = posMatrix;
leadfield2.leadfield = leadfield_mag;

figure(1);
ft_plot_mesh(headmodel.bnd(1), 'facecolor', 'white', 'edgecolor', 'none');
alpha(0.4);
camlight;
ft_plot_mesh(leadfield2.pos, 'vertexcolor', 'b');

% Initialize leadfield2.labels
leadfield2.labels = cell(size(leadfield2.pos, 1), 1);

% Loop through each position in leadfield2.pos to add labels
for j = 1:length(leadfield2.pos)
    position = leadfield2.pos(j, :);
    
    % Initialize an empty label for this position
    label = [];
    % Loop through each tissue (e.g., from 1 to the number of tissues)
    for i = 1:num_tissues
        matching_indices = find(sourcemodel2.tissue == i);
        for k = 1:length(matching_indices)
            sourcemodel_pos = sourcemodel.pos(matching_indices(k), :);
            if isequal(position, sourcemodel_pos)
                label = sourcemodel2.tissuelabel(i);
                break; 
            end
        end
        
        if ~isempty(label)
            break;
        end
    end    
    leadfield2.labels{j} = label;
end

% Initialize a counter to keep track of the last assigned position
last_assigned = 0;
num_tissues = length(sourcemodel2.tissuelabel);
combined_sourcemodel.pos = zeros(length(sourcemodel.pos),3);
num_cells = 0;

for tissue_label = 1:num_tissues
    tissue_indices = find(sourcemodel2.tissue == tissue_label);
    num_positions = length(tissue_indices);
    num_cells = num_cells + num_positions;
    if num_positions > 0
        combined_sourcemodel.pos(last_assigned+1 :last_assigned + num_positions, :) = sourcemodel.pos(tissue_indices, :);
        combined_sourcemodel.tissuelabel(last_assigned+1:last_assigned + num_positions) = sourcemodel2.tissuelabel(tissue_label);
        last_assigned = last_assigned + num_positions;
    end
end

nonEmptyIndices2 = find(~cellfun('isempty', leadfield2.labels));
leadfield_final.pos = leadfield2.pos(nonEmptyIndices2, :);
leadfield_final.leadfield = leadfield2.leadfield(nonEmptyIndices2, :);
leadfield_final.labels = leadfield2.labels(nonEmptyIndices2, :);

% Create a text label for each position with corresponding labels
for j = 1:length(leadfield_final.pos)
    position = leadfield_final.pos(j, :);
    label = leadfield_final.labels{j};
    text(position(1), position(2), position(3), label, 'FontSize', 10);
end

%% Sort the data according to the sequence set in atlas
sort_indices = [];

for i = 1:length(sourcemodel2.tissuelabel)
    for j = 1:length(leadfield_final.labels)
        if strcmp(sourcemodel2.tissuelabel{i}, leadfield_final.labels{j})
            sourcemodel2.tissuelabel{i};
            leadfield_final.labels{j};

            % If there's a match, store the index
            sort_indices = [sort_indices, j];
        end
    end
end

leadfield_final.pos = leadfield_final.pos(sort_indices, :);
leadfield_final.leadfield = leadfield_final.leadfield(sort_indices, :);
leadfield_final.labels = leadfield_final.labels(sort_indices, :);

%% Initialize a cell array to store unique labels
unique_labels = {};

for i = 1:length(leadfield_final.labels)
    current_label = leadfield_final.labels{i};
    if ~any(cellfun(@(x) isequal(x, current_label), unique_labels))
        unique_labels{end + 1} = current_label;
    end
end

% Count the number of unique labels
num_unique_labels = length(unique_labels);

%% Store the averaged positions and leadfield values
averaged_positions = [];
averaged_leadfield = [];
unique_labels = {};

% Initialize variables to keep track of the current label and accumulated data
current_label = leadfield_final.labels{1};
accumulated_positions = leadfield_final.pos(1, :);
accumulated_leadfield = leadfield_final.leadfield(1, :);
count = 1;

% Loop through labels starting from the second element
for i = 2:length(leadfield_final.labels)
    if strcmp(leadfield_final.labels{i}, current_label)
        accumulated_positions = accumulated_positions + leadfield_final.pos(i, :);
        accumulated_leadfield = accumulated_leadfield + leadfield_final.leadfield(i, :);
        count = count + 1;
    else
        averaged_positions(end + 1, :) = accumulated_positions / count;
        averaged_leadfield(end + 1, :) = accumulated_leadfield / count;
        unique_labels{end + 1} = current_label;
        current_label = leadfield_final.labels{i};
        accumulated_positions = leadfield_final.pos(i, :);
        accumulated_leadfield = leadfield_final.leadfield(i, :);
        count = 1;
    end
end

averaged_positions(end + 1, :) = accumulated_positions / count;
averaged_leadfield(end + 1, :) = accumulated_leadfield / count;
unique_labels{end + 1} = current_label;

% Update leadfield_final with averaged data
leadfield_final_avg.pos = averaged_positions;
leadfield_final_avg.leadfield = averaged_leadfield;
leadfield_final_avg.labels = unique_labels';

%% Compute the H matrix

% Load matrices
systemMatrix_path = ['BatchValidation_POMOESP_sub' num2str(subject_number) '.mat'];

% Load data for the specified subject
systemMatrix = load(systemMatrix_path);
L = transpose(leadfield_final_avg.leadfield);
C = systemMatrix.Ck;
n_electrodes = length(C(:,1));
n_latents = length(C(1,:));
n_sources = length(L(1,:));

%% l1 regularization including the BIC criterion
% Initialize H
H = zeros(n_sources, n_latents);

% Define the regularization parameter lambda
lambda_max = 0.000006;  % Max value sqp
lambda_values = linspace(lambda_max * 10^-4, lambda_max, 10);

% Define optimization options
% algorithm = 'sqp';
algorithm = 'interior-point';
options = optimoptions('fmincon', 'Algorithm', algorithm, 'Display', 'off');

% Define linear equality constraints + bounds for variables 
Aeq = [];
beq = [];
lb = [];
ub = [];

% Loop through each latent variable (source)
for source_idx = 1:n_latents

    H_source = zeros(n_sources, length(lambda_values));
    BIC_values = zeros(length(lambda_values), 1);

    % Define the objective function for Lasso for this source
    for j = 1:length(lambda_values)
        f = @(x) norm(C(:, source_idx) - L * x, 'fro')^2 + lambda_values(j) * norm(x, 1);

        % Initialize the solution (e.g., with zeros)
        x0 = zeros(size(L, 2), 1);

        % Solve the optimization problem using fmincon
        x = fmincon(f, x0, [], [], Aeq, beq, lb, ub, [], options);
        H_source(:, j) = x;

        % Calculate residuals for this solution
        residuals = C(:, source_idx) - L * x;

        % Calculate the covariance matrix of residuals
        Sigma_hat = cov(residuals);

        % Calculate the effective number of parameters (non-zero coefficients)
        k_lambda = sum(x ~= 0);

        % Calculate the BIC for this source
        n = size(C, 1); % Number of data points (samples)
        m = sum(x ~= 0); % Number of effective parameters
        BIC = -n * m - n * log(det(Sigma_hat)) + k_lambda * log(n);
        BIC_values(j) = BIC;
    end

    % Find the lambda that minimizes the BIC for this source
    [min_BIC, min_lambda_idx] = min(BIC_values);
    H(:, source_idx) = H_source(:, min_lambda_idx);

end

%% Compute source signals z(k)

z_k = H * xek;
