function [predict_data, out] = thot(data, opts)
%% Prepare for initialization
% Read parameters from the input structure |opts|.
%
% * |opts.predictor_type|: predictor operator, such as 'ARIMA', 'AR'
% * |opts.para|: Hyperparameters corresponding to the prediction operator, e.g. para=[p d q] corresponding to the ARIMA operator
% * |opts.phi|: regularisation parameter
% * |opts.isclearG|: Whether to zero out the non-diagonal elements of fft_S after each iteration. 1 is zeroed out, default is 0.
% * |opts.itermax|: ：最大迭代次数
% * |opts.core_size|: ：最大迭代次数
% * |opts.isReal|: ：输入数据是否为实数。如果为实数则为1，否则则为0，默认为0.
if ~isfield(opts, 'predictor_type'); opts.predictor_type = 'ARIMA'; end
if ~isfield(opts, 'para'); opts.para = [3 1 0]; end
if ~isfield(opts, 'phi'); opts.phi = 1e2; end
if ~isfield(opts, 'isclearG'); opts.isclearG = 0; end
if ~isfield(opts, 'itermax'); opts.itermax = 5; end
if ~isfield(opts, 'core_size'); opts.core_size = 4; end
if ~isfield(opts, 'isReal'); opts.isReal = 0; end

[n1, n2, n3] = size(data{1});
out = struct();
T = length(data);
fft_data = cell(1, T);
time = zeros(1,opts.itermax+1);
tic;
for t = 1: T
    fft_data{t} = fft(data{t}, [], 3);
end
fft_S = cell(1, T);
fft_U = zeros(n1, opts.core_size, n3);
fft_V = zeros(n2, opts.core_size, n3);
for i = 1: n3
    [U, S, V] = svd(fft_data{T}(:,:,i),'econ');
    fft_S{T}(:,:,i) = S(1:opts.core_size,1:opts.core_size);
    fft_U(:,:,i) = U(:,1:opts.core_size);
    fft_V(:,:,i) = V(:,1:opts.core_size);
end
for t=1:T-1
%     for i = 1: n3
%         [~, S, ~] = svd(fft_data{t}(:,:,i),'econ');
%         fft_S{t}(:,:,i) = S(1:opts.core_size,1:opts.core_size);
%     end
    fft_S{t} = slice_thr_dot(permute(fft_U,[2,1,3]),fft_data{t},fft_V);
end
time(1) = toc;

%% initialize the input; predictor is "ar.m", "var.m", "arima.m";
% update_para is "ar_update.m", "var_update.m", "arima_update.m";
update_para_type = [opts.predictor_type '_update'];
predictor = str2func(opts.predictor_type);           update_para = str2func(update_para_type);
predict_para = update_para(fft_S, opts.para);
% predict_para = update_para(cell_tensor_fdiag(fft_S), opts.para);


%% begin iteration
ftol = 0;
out.f_hist = zeros(1, opts.itermax);
out.f1_hist = zeros(1, opts.itermax);
out.f2_hist = zeros(1, opts.itermax);
out.fft_U_hist = cell(1, opts.itermax+1);
out.fft_V_hist = cell(1, opts.itermax+1);
out.fft_U_hist{1} = fft_U;  out.fft_V_hist{1} = fft_V;
for iter = 1: opts.itermax
    tic;
    % 记录当前函数值
    f_now2 = 0;     % 记录S与S_hat的差距
    f_now1 = 0;     % 记录X_t与U*S_t*V^H的差距
    for t = 1: T
        f_now2 = f_now2 + opts.phi * frob(fft_data{t}-...
            slice_thr_dot(fft_U,fft_S{t},permute(fft_V,[2,1,3])))^2;
        fft_S_hat = predictor(fft_S(1: t-1), opts.para, predict_para);
        if fft_S_hat ~= -999
            f_now1 = f_now1 + frob(fft_S{t}-fft_S_hat)^2;
        end
    end
    out.f1_hist(iter) = f_now1/n3;
    out.f2_hist(iter) = f_now2/n3;
    out.f_hist(iter) = (f_now1+f_now2)/n3;    
    
    % 判断迭代是否停止
    if iter > 1 && abs(out.f_hist(iter) - out.f_hist(iter-1)) / abs(out.f_hist(1)) < ftol
        break;
    end
    
    % update fft_S
    fft_U_H = permute(fft_U, [2,1,3]);
    for t = 1: T
        fft_S_hat = predictor(fft_S(1: t-1), opts.para, predict_para);
        G = slice_thr_dot(fft_U_H, fft_data{t}, fft_V);
        if fft_S_hat == -999
            fft_S{t} = G;
        else
            fft_S{t} = (fft_S_hat+opts.phi*G)/(1+opts.phi);
        end
    end
    
    % 若opts.isclearG=1则将fft_S的非对角线元素清零
    if opts.isclearG
        for t=1: T
            for i = 1: n3
                fft_S{t}(:,:,i) = diag(diag(fft_S{t}(:,:,i)));
            end
        end
    end
    
    % update fft_U
    for i = 1: n3
        SVD_matric = zeros(n1,opts.core_size);
        for t = 1: T
            SVD_matric = SVD_matric + fft_data{t}(:,:,i)*fft_V(:,:,i)*fft_S{t}(:,:,i)';
        end
        [L, ~, R] = svd(SVD_matric,'econ');
        fft_U(:,:,i) = L*R';
    end
    out.fft_U_hist{iter+1} = fft_U;
    
    % update fft_V
    for i = 1: n3
        SVD_matric = zeros(n2,opts.core_size);
        for t = 1: T
            SVD_matric = SVD_matric + fft_data{t}(:,:,i)'*fft_U(:,:,i)*fft_S{t}(:,:,i);
        end
        [L, ~, R] = svd(SVD_matric,'econ');
        fft_V(:,:,i) = L*R';
    end

    out.fft_V_hist{iter+1} = fft_V;

    
    % update predict_para
    predict_para = update_para(fft_S, opts.para);
    % predict_para = update_para(cell_tensor_fdiag(fft_S), opts.para);
    time(iter+1)=toc;
end

%% predict X_{T+1} and output the result
predict_fft_S = predictor(fft_S(1: T), opts.para, predict_para);
% predict_fft_data = slice_thr_dot(fft_U,predict_fft_S,permute(fft_V, [2,1,3]));
% predict_data = ifft(predict_fft_data,[],3);
predict_data = predictor(data(1:T), opts.para, predict_para);
out.fft_U = fft_U;  out.fft_V = fft_V;
out.U = ifft(fft_U,[],3);   out.V = ifft(fft_V,[],3);
G_bar = fft_S;  G_bar{T+1} = predict_fft_S;
out.fft_S = G_bar;
for t = 1: T+1
    G_bar{t} = ifft(G_bar{t},[],3);
end
out.G = G_bar;
out.predict_para = predict_para;
% out.iter = iter;
out.time = time;

end

function diag_G = cell_tensor_fdiag(G)  
    T = length(G);
    n = min(size(G{1},1),size(G{1},2));
    n3 = size(G{1},3);
    diag_G = cell(T,1);

    for t = 1:T
        G_t = reshape(diag(reshape(G{t}, n, [])), n, n3);
        diag_G{t} = G_t;
    end
%     for t = 1: T
%         G_t = zeros(n,n3);
%         for i = 1: n
%             G_t(i,:) = G{t}(i,i,:);
%         end
%         diag_G{t} = G_t;
%     end
end

