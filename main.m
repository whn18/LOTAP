clear
clc
%% add path
addpath(genpath(pwd))

%% generate synthetic (SYN) data
fprintf('*** Generate Synthetic Data ****\n')
T=200;
scale_data=[100,100,10];
scale_core=4;
AR_init=[-0.674363671145833;0.0420496915437676;0.723935677619934];
AR_evolution=zeros(T,length(AR_init));
for i=1:length(AR_init)
    AR_evolution(:,i)=AR_init(i)+0.1*sin((1:T)'*2*pi/200);
end
d=1;
gamma = 0.1;

[X,fft_S,fft_U,fft_V]=TSVD_generate_ARIMA(T,scale_data,scale_core,AR_evolution,d);
% plot(1:T,cellfun(@frob,X))

% add noisy
X_without_noise = X;
for i = 1: T
    new_TDS_data = X{i} + gamma*frob(X{i}/sqrt(prod(scale_data)))*rand(scale_data);
    X{i} = new_TDS_data;
end

%% initial the parameter
trans_length = 80;
predict_length = 20;
predictor_type = 'ARIMA';
core_size = 4;
para = [3 1 0];
phi = 100;
itermax = 3;

%% LOTAP
fprintf('**** Begin LOTAP ****\n')
opts = struct();    
opts.predictor_type = 'ARIMA';  opts.para = para;   opts.phi = phi; 
opts.isclearG = 0;  opts.core_size = 6;    opts.itermax = itermax;
predict_X = cell(1, predict_length);
for t=1: predict_length
    [predict_data, out] = LOTAP(X(1: trans_length+t-1), opts);
    predict_X{t} = predict_data;
end
% LOTAP_err = geterr(predict_X,X_without_noise(trans_length+1:trans_length+predict_length));
% mean(LOTAP_err)
LOTAP_err2 = geterr(predict_X,X(trans_length+1:trans_length+predict_length));
fprintf('MSPE = %5.3f, time cost per time step = %5.3fs.\n',mean(LOTAP_err2), mean(out.time))

%% 
rmpath(genpath(pwd))