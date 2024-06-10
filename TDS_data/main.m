%% load data
clear
clc
load tds_real.mat;
%% initial the parameter
trans_length = 20;
predict_length = 5;
predictor_type = 'ARIMA';
core_size = 4;
para = [3 1 0];
phi = 100;
itermax = 3;

%% 有噪声数据与无噪声数据对比
% err = geterr(X(trans_length+1:trans_length+predict_length),X_without_noise(trans_length+1:trans_length+predict_length));
err = geterr(X_without_noise(trans_length+1:trans_length+predict_length),X(trans_length+1:trans_length+predict_length));

%% 直接用init_out.predict_para进行预测
% predict_X = ARIMA(X(1:20),[3 1 0],init_out.predict_para);
% error1 = frob(predict_X-X_without_noise{20+1})/frob(X_without_noise{20+1});

%% 用ARIMA_update对X计算出predict_para后再进行预测
predict_X2 = cell(1, predict_length);
for t=1: predict_length
    tic;
    predict_para2 = ARIMA_update(X(1:trans_length+t-1),[3 1 0]);
    predict_X2{t} = ARIMA(X(1:trans_length+t-1),[3 1 0],predict_para2);
    ARIMA_RUNTIME = toc;
end
error1 = geterr(predict_X2,X_without_noise(trans_length+1:trans_length+predict_length));
mean(error1)
error2 = geterr(predict_X2,X(trans_length+1:trans_length+predict_length));
mean(error2)

%% 用ARIMA_update对fft_G计算出predict_para后直接对X预测
% predict_para3 = ARIMA_update(fft_G(1:20),[3 1 0]);
% predict_X3 = ARIMA(X(1:20),[3 1 0],predict_para3);
% error3 = frob(predict_X3-X{20+1})/frob(X{20+1});

%% 用ARIMA_update对fft_G计算出predict_para后预测
% predict_para3 = ARIMA_update(init_out.fft_G(1:20),para);
% predict_fft_G = ARIMA(init_out.fft_G(1: 20), [3 1 0], predict_para3);
% predict_fft_data = slice_thr_dot(init_out.fft_U,predict_fft_G,permute(init_out.fft_V, [2,1,3]));
% predict_X4 = ifft(predict_fft_data,[],3);
% error4 = frob(predict_X4-X_without_noise{20+1})/frob(X_without_noise{20+1});

%% thot(phi=100)
phi = 10;
opts = struct();    
opts.predictor_type = 'ARIMA';  opts.para = para;   opts.phi = phi; 
opts.isclearG = 0;  opts.core_size = 6;    opts.itermax = itermax;
predict_X = cell(1, predict_length);
for t=1: predict_length
    [predict_data, out] = thot(X(1: trans_length+t-1), opts);
    predict_X{t} = predict_data;
end
thot_err = geterr(predict_X,X_without_noise(trans_length+1:trans_length+predict_length));
mean(thot_err)
thot_err2 = geterr(predict_X,X(trans_length+1:trans_length+predict_length));
mean(thot_err2)

%% MCAR(small core)
MCAR_predict = cell(1, predict_length);     MCAR_Runtime = zeros(1, predict_length);
itermax = 3;
for i=1: predict_length
    tic;
    [G_initial, predict_para, U_initial] = ostp_initial(X(1:trans_length+i-1), predictor_type, [4,4,10], para, phi, itermax);
    G_predict = ARIMA(G_initial, para, predict_para);   
    MCAR_predict{i} = tmprod(G_predict, U_initial, 1:3);
    MCAR_Runtime(i)=toc;
end
MCAR_ERR1 = geterr(MCAR_predict, X(trans_length+1:trans_length+predict_length));
mean(MCAR_ERR1)

%% MCAR(large core)
MCAR_predict2 = cell(1, predict_length);     MCAR_Runtime2 = zeros(1, predict_length);
itermax = 3;
for i=1: predict_length
    tic;
    [G_initial, predict_para, U_initial] = ostp_initial(X(1:trans_length+i-1), predictor_type, [40,40,10], para, phi, itermax);
    G_predict = ARIMA(G_initial, para, predict_para);   
    MCAR_predict2{i} = tmprod(G_predict, U_initial, 1:3);
    MCAR_Runtime2(i)=toc;
end
MCAR_ERR2 = geterr(MCAR_predict2, X(trans_length+1:trans_length+predict_length));
mean(MCAR_ERR2)
% E2 = geterr(MCAR_predict2, X_without_noise(trans_length+1:trans_length+predict_length));
% mean(E2)
