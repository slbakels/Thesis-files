% Thesis Sanna Bakels, 4480279
% Step 4: Perform the connectivity analysis

clear; close all; clc;

% Using the Matlab Toolbox EEGLAB from:


%% Load matrices
% Define the subject number o
subject_number = 3;
XX = 0.309;

% Construct file paths based on the subject number
leadfield_path = fullfile('~/Matlab/Lead field', ['leadfield_subject' num2str(subject_number) '_sorted.mat']);
input_path = fullfile('~/Matlab/ExperimentData/ERPs', ['input_subject' num2str(subject_number) '.mat']);
systemMatrix_path = ['BatchValidation_POMOESP_sub' num2str(subject_number) '.mat'];

% Load data for the specified subject
leadfield_subject = load(leadfield_path);
input_subject = load(input_path);
systemMatrix = load(systemMatrix_path);

%% Extract relevant data
L = transpose(leadfield_subject.leadfield_final_avg.leadfield);
C = systemMatrix.Ck;
H = leadfield_subject.H;
n_electrodes = length(C(:,1));
n_latents = length(C(1,:));
n_sources = length(L(1,:));

%% Granger causality analysys
% Estimation of noise covariance in the source dynamic
% Compute Sigma_w and Sigma_e from the state matrices determined with SID

% Load system matrices
A = systemMatrix.Ae;
B = systemMatrix.Be;
C = systemMatrix.Ce;
D = systemMatrix.De;
x = systemMatrix.xek';
y = systemMatrix.y7';

time = systemMatrix.time;
ts = 1/2048;
u = input_subject.(['input_subject' num2str(subject_number)])';

% Form state-space model
sys = ss(A, B, C, D, ts);
sys_K = ss(systemMatrix.Ak, systemMatrix.Bk, systemMatrix.Ck, systemMatrix.Dk, ts);

[y_sys, t, x_sys] = lsim(sys, u, time);

% Estimate of the state noise covariance
residual_w = systemMatrix.xek - x_sys ;
Sigma_w = cov(residual_w);

% Estimate of the measurement noise covariance
residual_e = systemMatrix.y7 - systemMatrix.yek;
Sigma_e = cov(residual_e);

%% Perform optimization problem for determining Sigma_eta and Sigma_v
% Inititiate matrices
Sigma_eta = zeros(n_sources, n_sources);
Sigma_v = ones(n_electrodes, n_electrodes);

% Define the objective function
objectiveFunction = @(x) norm(Sigma_e - L * (x(1)*eye(size(Sigma_eta)))*L' - diag(x(2:end)), 'fro');

% Initial guess for the optimization variables (including α and diagonal elements of Σv)
alpha_initial_guess = 0.01;
diag_initial_guess =ones(n_electrodes, 1);
initialGuess = [alpha_initial_guess; diag_initial_guess]; 

% Define lower bounds for variables
% Ση ⪰ 0 (element-wise) ;  % Σv ⪰ 0 (element-wise) and Σv is diagonal
% α must be >= epsilon, and all diagonal elements of Σv must be >= 0
epsilon = 1e-3;
lb = [epsilon; zeros(n_electrodes, 1)];

% Solve the optimization problem
options = optimoptions('fmincon', 'StepTolerance', 1e-6, 'Display', 'iter');
[x_opt, fval, exitflag, output] = fmincon(objectiveFunction, initialGuess, [], [], [], [], lb, [], [], options);

% Extract the optimized values
alpha_opt = x_opt(1);
Sigma_eta = alpha_opt * eye(size(Sigma_eta));
Sigma_v = x_opt(2:end);
Sigma_v = diag(Sigma_v);

% Calculate the optimized Sigma_e
Sigma_e_opt = L * Sigma_eta * L' + Sigma_v;

%% Granger causality analysis for the full model
P = dare(A, H', Sigma_w, Sigma_eta);
cov_output = H*P*H' + Sigma_eta;
F = zeros(n_sources, n_sources);

for j = 1:n_sources
    % Remove the j-th row and corresponding column from H and Sigma_eta
    H_r = [H(1:j-1, :); H(j+1:end, :)];
    Sigma_eta_r = Sigma_eta([1:j-1, j+1:end], [1:j-1, j+1:end]);
    
    % Calculate the DARE solution with the modified matrices
    P_r = dare(A, H_r', Sigma_w, Sigma_eta_r);
   
    % Calculate the covariance matrix based on the modified matrices
    cov_output_r = H_r * P_r * H_r' + Sigma_eta_r;

    for i = 1:n_sources
        if i == j
            F(i, j) = NaN;
        elseif i < j
            F(i, j) = log(cov_output_r(i, i) / cov_output(i, i));
        else
            F(i, j) = log(cov_output_r(i-1, i-1) / cov_output(i, i));
        end
    end
end

%% Create a Granger causality plot
labelStrArray = string(leadfield_subject.leadfield_final_avg.labels);
labelStrArray = strrep(labelStrArray, '_', ' ');
source_indices = 1:length(labelStrArray);
threshold = XX * max(abs(F(:)));

% Plot grid positions
grid_positions = leadfield_subject.leadfield_final_avg.pos;
x = grid_positions(:, 1);
y = grid_positions(:, 2);
z = grid_positions(:, 3);

figure(1);
imagesc(abs(F));
colormap('jet');
colorbar;
caxis([0, threshold]);

% Set the current axes explicitly
set(gca, 'XTick', source_indices);
set(gca, 'XTickLabel', labelStrArray);
set(gca, 'YTick', source_indices);
set(gca, 'YTickLabel', labelStrArray);

%% Create a connectivity diagram
indices = 1:length(leadfield_subject.leadfield_final_avg.pos);
num_points = length(labelStrArray)+1;

% Define the angles for evenly spaced points around the circle
angles = linspace(0, 2*pi, num_points);
radius = 1.2; % Increase the radius to place labels outside the circle

% Calculate the coordinates of the points
x = radius * cos(angles);
y = radius * sin(angles);

% Create a figure
figure(2);

% Plot the circle, each region in a distinct colour
plot(x(1:2), y(1:2), 'r.', 'MarkerSize', 10) % red
hold on
plot(x(3:17), y(3:17), 'm.', 'MarkerSize', 10); % pink
plot(x(18:19), y(18:19), '.', 'MarkerSize', 10, 'MarkerEdgeColor', "#006400", 'MarkerFaceColor',  "#006400"); % dark green
plot(x(20:24), y(20:24), '.', 'MarkerSize', 10, 'MarkerEdgeColor',	"#EDB120", 'MarkerFaceColor', 	"#EDB120"); % yellow
plot(x(25:35), y(25:35), '.', 'MarkerSize', 10, 'MarkerEdgeColor',	"#D95319", 'MarkerFaceColor', 	"#D95319"); % orange
plot(x(36:46), y(36:46), '.', 'MarkerSize', 10, 'MarkerEdgeColor',	"#7E2F8E", 'MarkerFaceColor', 	"#7E2F8E"); % purple
plot(x(47:52), y(47:52), 'b.', 'MarkerSize', 10, 'MarkerEdgeColor',	"#4DBEEE", 'MarkerFaceColor', 	"#4DBEEE"); % light blue
plot(x(53:59), y(53:59), 'b.', 'MarkerSize', 10); % blue
plot(x(60:end-1), y(60:end-1), '.', 'MarkerSize', 10, 'MarkerEdgeColor',"#77AC30", 'MarkerFaceColor',	"#77AC30"); % green

% Add labels outside the circle
text_offset = 0.05; % offset for the labels 

for i = 1:num_points-1
    label_x = x(i) + text_offset * cos(angles(i));
    label_y = y(i) + text_offset * sin(angles(i));
    
    text(label_x, label_y, labelStrArray{i}, 'Rotation', 360 * i / num_points, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
arrowhead_size = 0.05; % Adjust this value to your preference
indices_threshold = [];

% Draw arrows to indicate a connection
for i = indices(1):indices(end)
    for j = indices(1):indices(end)
        if i ~= j && abs(F(i, j)) > threshold
            indices_threshold = [indices_threshold; i, j];
            p1 = [x(j) y(j)];
            p2 = [x(i) y(i)];
            dp = p2 - p1;
            magnitude = norm(dp);
            norm_dp = dp / magnitude;
            arrow_length = magnitude;
            quiver(p1(1), p1(2), norm_dp(1), norm_dp(2), arrow_length, 'AutoScaleFactor', arrowhead_size, 'HandleVisibility', 'off');
        end
    end
end

% Add the legend without labels for arrows
legend('Central sulcus', 'Frontal lobe', 'Insular lobe', 'Limbic lobe', 'Occipital lobe', 'Parietal lobe', 'Sub cortical grey nuclei', 'Temporal lobe', 'Cerebellum', 'Location', 'Best');
