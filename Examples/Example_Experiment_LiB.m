clear all
close all

%% Load experimental data and specify frequency domain f and time domain t profiles
% Experimental Lithium ion battery impedance data
% from paper: https://doi.org/10.1016/j.electacta.2015.09.097
load LiB_1.mat

% Full range case

Z_real_exp = Z_prime;
Z_imag_exp = Z_double_prime;
f = freq;

%Full range, half sampling density

% Z_real_exp = Z_prime(1:2:end);
% Z_imag_exp = Z_double_prime(1:2:end);
% f = freq(1:2:end);

% Reduced range case
% 
% Z_real_exp = Z_prime(11:97);
% Z_imag_exp = Z_double_prime(11:97);
% f = freq(11:97);

% Reduced range, half sampling density

% Z_real_exp = Z_prime(11:2:97);
% Z_imag_exp = Z_double_prime(11:2:97);
% f = freq(11:2:97);

% NF: number of frequency sampling points
% NT: number of time domain sampling points
NF = length(f);
NT = 120;
t = logspace(-6, 6, NT)';

%% Calculate basis matrix A_real, A_imag
% Given f and t, get A_real, A_imag via cal_Basis(f,t) function
% cal_Basis(f,t) can take some time, A can be pre-computed and stored

% tic
% [A_real, A_imag] = cal_Basis(f,t);
% toc

load exp_A_real_FR_FS.mat; load exp_A_imag_FR_FS.mat; % pre-computed basis matrix for full range, 10ppd;
% load exp_A_real_FR_HS.mat; load exp_A_imag_FR_HS.mat; % pre-computed basis matrix for full range, 5ppd;
% load exp_A_real_RR_FS.mat; load exp_A_imag_RR_FS.mat; % pre-computed basis matrix for reduced range, 10ppd;
% load exp_A_real_RR_HS.mat; load exp_A_imag_RR_HS.mat; % pre-computed basis matrix for reduced range, 5ppd;

%% Deconvolution
lambda = logspace(-10,1,100); %set-up grid of shrinkage tuning parameter
%model = sms_DRT(Z_real_exp,Z_imag_exp,A_real,A_imag,lambda,0);
model = sms_DRT(Z_real_exp,Z_imag_exp,A_real,A_imag,lambda,1,f);% considering high-frequency inductance 

%% Estmation results
R_infy_est = model.R_infy;
R_p_est = model.R_p;
DRT_est = model.beta;
L_est = model.inductance
Z_real_est = model.Z_real;
Z_imag_est = model.Z_imag;

figure(1)
plot(Z_real_exp, -Z_imag_exp, '-o')
hold on
plot(Z_real_est, -Z_imag_est, '-*')
axis equal
legend('Truth','Estimation')
xlabel('Re(Z)/\Omega')
ylabel('-Im(Z)/\Omega')

tc = (t(2:end) + t(1:end-1))./2; % using center of inteval [t_m, t_{m+1})
figure(2)
semilogx(tc, DRT_est, '-x')
xlabel('\tau /S')
ylabel('G(\tau)')

