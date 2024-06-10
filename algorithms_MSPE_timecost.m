clear
clc
% load TDS_data\tds.mat
% load USHCN_data\USHCN.mat
load NASDAQ100_data\NASDAQ.mat
% load CCDS_data\CCDS.mat

%% TSVD-AR
predict_TSVD_AR = cell(1, predict_length);
Timecost_TSVD_AR = zeros(1, predict_length);
std_opts.itermax = 3;
for t=1: predict_length
    tstart=tic;
    [predict_data, out] = thot(X(1: trans_length+t-1), std_opts);
    Timecost_TSVD_AR(t) = toc(tstart);
    predict_TSVD_AR{t} = predict_data;
end
MSPE_TSVD_AR = geterr(predict_TSVD_AR, X(trans_length+1:trans_length+predict_length));
% MSPE_TSVD_AR = mean(MSPE_TSVD_AR);
mean(MSPE_TSVD_AR)
Timecost_TSVD_AR = mean(Timecost_TSVD_AR);

%% classical ARIMA
% predict_ARIMA = cell(1, predict_length);
% Timecost_ARIMA = zeros(1, predict_length);
% para = std_opts.para;
% for t=1: predict_length
%     tic;
%     predict_para = ARIMA_update(X(1:trans_length+t-1),para);
%     predict_ARIMA{t} = ARIMA(X(1:trans_length+t-1),para,predict_para);
%     Timecost_ARIMA(t) = toc;
% end
% MSPE_ARIMA = geterr(predict_ARIMA, X(trans_length+1:trans_length+predict_length));
% MSPE_ARIMA = mean(MSPE_ARIMA);
% Timecost_ARIMA = mean(Timecost_ARIMA);

%% classical point ARIMA
predict_Point_ARIMA = cell(1, predict_length);
Timecost_Point_ARIMA = zeros(1, predict_length);
para = std_opts.para;
for t=1: predict_length
    tic;
    predict_para = POINT_ARIMA_update(X(1:trans_length+t-1),para);
    predict_Point_ARIMA{t} = POINT_ARIMA(X(1:trans_length+t-1),para,predict_para);
    Timecost_Point_ARIMA(t) = toc;
end
MSPE_Point_ARIMA = geterr(predict_Point_ARIMA, X(trans_length+1:trans_length+predict_length));
% MSPE_Point_ARIMA = mean(MSPE_Point_ARIMA);
Timecost_Point_ARIMA = mean(Timecost_Point_ARIMA);

%% MOAR
Timecost_MOAR = zeros(1, predict_length);
predict_MOAR = cell(1, predict_length);  
itermax = 3;
core_size = [4,4,fix(size(X{1},3)/2)];
% core_size = [4,4,size(X{1},3)];
for t=1: predict_length
    tic;
    [G_initial, predict_para, U_initial] = ostp_initial(X(1:trans_length+t-1), 'ARIMA', core_size, std_opts.para, 1e-2, itermax);
    G_predict = ARIMA(G_initial, para, predict_para);   
    predict_MOAR{t} = tmprod(G_predict, U_initial, 1:3);
    Timecost_MOAR(t)=toc;
end
MSPE_MOAR = geterr(predict_MOAR, X(trans_length+1:trans_length+predict_length));
% MSPE_MOAR = mean(MSPE_MOAR);
Timecost_MOAR = mean(Timecost_MOAR);

%% MCAR
Timecost_MCAR = zeros(1, predict_length);
predict_MCAR = cell(1, predict_length);
itermax = 3;
core_size = [4,4,fix(size(X{1},3)/2)];
% core_size = [4,4,size(X{1},3)];
for t=1: predict_length
    tic;
    [G_initial, predict_para, U_initial] = ostp_initial(X(1:trans_length+t-1), 'ARIMA', core_size, std_opts.para, 1e2, itermax);
    G_predict = ARIMA(G_initial, para, predict_para);   
    predict_MCAR{t} = tmprod(G_predict, U_initial, 1:3);
    Timecost_MCAR(t)=toc;
end
MSPE_MCAR = geterr(predict_MCAR, X(trans_length+1:trans_length+predict_length));
% MSPE_MCAR = mean(MSPE_MCAR);
Timecost_MCAR = mean(Timecost_MCAR);

%% BHT-ARIMA
tt = 3;     
bht_core_size = [4,4,fix(size(X{1},3)/2),tt];    
% bht_core_size = [4,4,size(X{1},3),tt];   
itermax = 3; 
predict_BHT_ARIMA = cell(1, predict_length);    Timecost_BHT_ARIMA = zeros(1, predict_length);
BHT_X = cell(1, predict_length + trans_length);
for i = tt: (predict_length+trans_length)
    BHT_X{i} = zeros([size(X{i}) tt]);
    for j = 1:tt
        BHT_X{i}(:,:,:,j) = X{i-j+1};
    end
end
for i=1: predict_length
    predict_BHT_ARIMA{i} = X{1};
    tic;
    [G_initial, predict_para, U_initial] = ostp_initial(BHT_X(tt:trans_length+i-1), 'ARIMA', bht_core_size, std_opts.para, 1e3, itermax);
    G_predict = ARIMA(G_initial, para, predict_para);   
    BHT_predict = tmprod(G_predict, U_initial, 1:4);
    predict = BHT_predict(:,:,:,1);
    predict_BHT_ARIMA{i}(:,:,:) = BHT_predict(:,:,:,1); 
    Timecost_BHT_ARIMA(i)=toc;
end
MSPE_BHT_ARIMA = geterr(predict_BHT_ARIMA, X(trans_length+1:trans_length+predict_length));
% MSPE_BHT_ARIMA = mean(MSPE_BHT_ARIMA);
Timecost_BHT_ARIMA = mean(Timecost_BHT_ARIMA);
