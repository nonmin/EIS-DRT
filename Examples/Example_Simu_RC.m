clear all
close all

%% Specify frequency domain f and time domain t profiles
% NF: number of frequency sampling points
% NT: number of time domain sampling points
NF = 70;
%NF_h = 35; % 5 ppd

f = logspace(-1, 6, NF)'; % normal range, 10 ppd
%f = logspace(-1,6, NF_h)'; % half range, 5ppd
%f = logspace(0, 3.5, NF_h)'; % half range, 10 ppd

NT = 240;
t = logspace(-6, 6, NT)';

%% Calculate basis matrix A_real, A_imag
% Given f and t, get A_real, A_imag via cal_Basis(f,t) function
% cal_Basis(f,t) can take some time, A can be pre-computed and stored
% 
% tic
% [A_real, A_imag] = cal_Basis(f,t);
% toc

load simu_A_real_FR_NF.mat; load simu_A_imag_FR_NF.mat; % pre-computed basis matrix for normal range, 10 ppd;
%load simu_A_real_FR_NFh.mat; load simu_A_imag_FR_NFh.mat; % pre-computed basis matrix for normal range, 5 ppd;
%load simu_A_real_HR_NF.mat; load simu_A_imag_HR_NF.mat; % pre-computed basis matrix for half range, 10 ppd;

%% Set-up impedance model
R_infy_simu = 1; %ohm, high-frequency cut-off resistance
R_p_simu = 1; %ohm, overall polarization resisitance
C = 1; % F

% analytical impedace data
Z_simu = R_infy_simu + R_p_simu./(1+2i*pi*R_p_simu*C.*f);
Z_real_simu = real(Z_simu);
Z_imag_simu = imag(Z_simu);

%% Deconvolution
lambda = logspace(-10,1,100); %set-up grid of shrinkage tuning parameter
model = sms_DRT(Z_real_simu,Z_imag_simu,A_real,A_imag,lambda,0);

%% Estmation results
R_infy_est = model.R_infy;
R_p_est = model.R_p;
DRT_est = model.beta;
Z_real_est = model.Z_real;
Z_imag_est = model.Z_imag;

figure(1)
plot(Z_real_simu, -Z_imag_simu, '-o')
hold on
plot(Z_real_est, -Z_imag_est, '-*')
legend('Truth','Estimation')
xlim([0.8,1.8])
axis equal
xlabel('Re(Z)/\Omega')
ylabel('-Im(Z)/\Omega')

tc = (t(2:end) + t(1:end-1))./2; % using center of inteval [t_m, t_{m+1})
figure(2)
semilogx(tc, DRT_est, '-x')
xlabel('\tau /S')
ylabel('G(\tau)')

