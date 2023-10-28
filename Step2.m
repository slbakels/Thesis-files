% Thesis Sanna Bakels, 4480279
% Step 2: Perform the subspace identification methocs PO-MOESP and N4SID

clear; close all; clc;

% Using the Matlab System Identification toolbox and the LTI toolbox from:
% Verhaegen, M., Verdult, V., & Bergboer, N. (2007). Filtering and system identification: 
% an introduction to using matlab software. Delft University of Technology, 68, 163.

%% Import data
% Define the path to the folder containing the .mat file
inputFolderPath = '~/Matlab/ExperimentData/ERPs'; 

subjectNumbers = [3, 4, 5, 6, 8, 9, 10];

% Loop through the subject numbers
for i = 1:numel(subjectNumbers)
    % Generate the file name using the current subject number
    fileNameOutput = sprintf('ERP_subject%d.mat', subjectNumbers(i));
    fileNameInput = sprintf('input_subject%d.mat', subjectNumbers(i));
    fullFilePathOutput = fullfile(inputFolderPath, fileNameOutput);
    fullFilePathInput = fullfile(inputFolderPath, fileNameInput);
    
    % Load the .mat output files and store it in dynamically named variables
    variableNameOutput = ['output_subject' num2str(subjectNumbers(i))];
    assignin('base', variableNameOutput, load(fullFilePathOutput));

    % Load the .mat input files and store it in dynamically named variables
    variableNameInput = ['input_subject' num2str(subjectNumbers(i))];
    assignin('base', variableNameInput, load(fullFilePathInput));
end

%% Load the matrices
time = output_subject3.time;
Ts = 1/2048;

u1 = input_subject4.input_subject4';    y1 = output_subject4.averageData';
u2 = input_subject3.input_subject3';    y2 = output_subject3.averageData';  
u3 = input_subject5.input_subject5';    y3 = output_subject5.averageData';
u4 = input_subject6.input_subject6';    y4 = output_subject6.averageData';
u5 = input_subject8.input_subject8';    y5 = output_subject8.averageData';
u6 = input_subject9.input_subject9';    y6 = output_subject9.averageData';
u7 = input_subject10.input_subject10';  y7 = output_subject10.averageData';

%% PO-MOESP method - multiple batch method
% Step 1 - Data compression and order estimation
s = 20; % Number of block rows, a higher value does not work in the matlab functions

% Calculate the singular values
[S1,R1]=dordpo(u1,y1,s);
[S2,R2]=dordpo(u2,y2,s,R1);
[S3,R3]=dordpo(u3,y3,s,R2);
[S4,R4]=dordpo(u4,y4,s,R3);
[S5,R5]=dordpo(u5,y5,s,R4);
[S6,R6]=dordpo(u6,y6,s,R5);

% Plot of the singular values
figure(1)
subplot(1, 6, 1);   semilogy(1:s,S1,'x');
subplot(1, 6, 2);   semilogy(1:s,S2,'x');
subplot(1, 6, 3);   semilogy(1:s,S3,'x');
subplot(1, 6, 4);   semilogy(1:s,S4,'x');
subplot(1, 6, 5);   semilogy(1:s,S5,'x');
subplot(1, 6, 6);   semilogy(1:s,S6,'x');

%% Step 2 - Estimation of A and C
n = 15; % Order of the model determined from the singular values
[Ae,Ce,Ke] = dmodpo(R4,n); % Estimate of A and C

%% Step 3 - Estimation of B, D and the initial state
[Be,De] = dac2bd(Ae,Ce, u1,y1, u2,y2, u3,y3, u4,y4);

%% Step 4 - Model simulation and validation on a new batch 
% Model simulation
Ak = Ae-Ke*Ce;
Bk = [Be-Ke*De, Ke];
Ck = Ce;
Dk = [De zeros(size(De,1))];
x0k = dinit(Ak, Bk, Ck, Dk, [u7 y7], y7);
[yek, xek] = dltisim(Ak, Bk, Ck, Dk, [u7 y7], x0k);

VAF_PO = vaf(y7, yek);
vaf_mean = mean(VAF_PO, 1);
RMSE_mean = mean(rmse(y7, yek));

% Plot the results
figure(2)
plot(time, y5(:, 1))
hold on
plot(time, yek(:, 1))

%% Check assumptions PO-MOESP
disp(['The order of the model is: ', num2str(n)])
% Observability check
obs_mat = obsv(Ae, Ce);
obs = rank(obs_mat);
disp(['The rank of the observability matrix is: ', num2str(obs)])
% Controllability check
ctrb_mat = ctrb(Ak, Bk);
contr = rank(ctrb_mat);
disp(['The rank of the observability matrix is: ', num2str(contr)])

% Check if eigenvalues are inside the unit circle
eigenvalues = eig(Ak);
if all(abs(eigenvalues) < 1)
    disp('All eigenvalues of A-KC are inside the unit circle.');
else
    disp('At least one eigenvalue of A-KC is outside the unit circle.');
end

%% N4SID method - multiple batches
Ybatch = [y1; y2; y3; y4; y5; y6];
Ubatch = [u1; u2; u3; u4; u5; u6];

% Apply N4SID toolbox
z = iddata(Ybatch, Ubatch, Ts);
nx = 15;
[sys, x0] = n4sid(z, nx, 'Ts', Ts);

% Model simulation and validation
Ae = sys.A; Be = sys.B; Ce = sys.C; De = sys.D; Ke = sys.K;
Ak = Ae-Ke*Ce;
Bk = [Be-Ke*De, Ke];
Ck = Ce;
Dk = [De zeros(size(De,1))];
x0k = dinit(Ak, Bk, Ck, Dk, [u7 y7], y7); % To check if x0 is x0k --> is the case
[yk, xk] = dltisim(Ak, Bk, Ck, Dk, [u7 y7], x0k);

VAF_PO = vaf(y7, yk);
vaf_mean = mean(VAF_PO, 1)

% Plot the results
figure(3)
plot(time, y7(:, 1))
hold on
plot(time, yk(:, 1))

%% Check assumptions N4SID
contr = rank(ctrb_mat);

% Observability check
obs_mat = obsv(Ae, Ce);
obs = rank(obs_mat);
disp(['The rank of the observability matrix is: ', num2str(obs)])

% Controllability check
ctrb_mat = ctrb(Ak, Bk);
disp(['The rank of the controllability matrix is: ', num2str(contr)])

% Check if eigenvalues are inside the unit circle
eigenvalues = eig(Ak);
if all(abs(eigenvalues) < 1)
    disp('All eigenvalues of A-KC are inside the unit circle.');
else
    disp('At least one eigenvalue of A-KC is outside the unit circle.');
end
