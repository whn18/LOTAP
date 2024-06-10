clear
clc
T=500;
scale_data=[100,100,10];
scale_core=4;
AR_init=[-0.674363671145833;0.0420496915437676;0.723935677619934];
% AR_init=[-0.5707 -0.2281 0.8581];
% AR_init=[-1.13,-0.55,0.313];
AR_evolution=zeros(T,length(AR_init));
for i=1:length(AR_init)
    AR_evolution(:,i)=AR_init(i)+0.1*sin((1:T)'*2*pi/200);
end
d=1;
gamma = 0.1;

[X,fft_S,fft_U,fft_V]=TSVD_generate_ARIMA(T,scale_data,scale_core,AR_evolution,d);
plot(1:100,cellfun(@frob,X))

%% add noisy
X_without_noise = X;
for i = 1: T
    new_TDS_data = X{i} + gamma*frob(X{i}/sqrt(prod(scale_data)))*rand(scale_data);
    X{i} = new_TDS_data;
end

init_out = struct();
init_out.fft_U = fft_U;     init_out.fft_V = fft_V;     init_out.fft_S = fft_S;
init_out.phi = 100;
predict_para = cell(1,3);   predict_para{1} = AR_init;
para = [3 1 0];
init_out.predict_para = predict_para;
f1 = 0; f2 = 0;
for t = 1: T
    f2 = f2 + init_out.phi * frob(X{t}-...
        slice_thr_dot(fft_U,fft_S{t},permute(fft_V,[2,1,3])))^2;
    fft_S_hat = ARIMA(fft_S(1: t-1), para, predict_para);
    if fft_S_hat ~= -999
        f1 = f1 + frob(fft_S{t}-fft_S_hat)^2;
    end
end
init_out.f1 = f1;    init_out.f2 = f2;
init_out.f = (f1+f2)/scale_data(3);

save('tds.mat','init_out','X_without_noise','X');
