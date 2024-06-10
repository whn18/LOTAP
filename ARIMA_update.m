function predict_para = ARIMA_update(data, para)
% ARIMA model : x_t = a_1*x_{t-1}+...+a_p*x_{t-p}-b_1*e_{t-1}-...-b_q*e_{t-q}+e_t
% input:    data is a cell vector of tensor.
%           para = [p d q] is the order of the autoregressive model.
% output:   predict_para = {[a_1, a_2, ..., a_p], [b_1, ..., b_q], {e_1,...,e_T}}, 
%                        where a_i, b_j are the coefficients of AR and MA,
%                        T = (length(data_diff) - p), e_i is the residual of
%                        model.

%initialize
p = para(1);    d = para(2);    q = para(3);
data_length = length(data);
predict_para = cell(1, 3);

% tensor difference
data_diff = data;   % data_diff is the d order difference of data
for i = 1: d
    for j = 1: data_length - i
        data_diff{j} = data_diff{j+1} - data_diff{j};
    end
end
data_diff = data_diff(1: data_length - d);

% calculate the coefficients of AR, i.e. a_i
data_diff_length = length(data_diff);
autocorr = zeros(1, p+1);
for l = 0: p    % calculate the auto related coefficient 
    product = 0;
    for t = 1: data_diff_length-l
        product = product + sum(data_diff{t}.*conj(data_diff{t+l}), 'all');
    end
%     autocorr(l+1) = product;
    autocorr(l+1) = product/(data_diff_length-l);
end
R = double(toeplitz(autocorr(1: p)));   % R is the sample autocorrelation matrix
a = linsolve(R, autocorr(2:p+1)');  % R*a = autocorr, where a is the coefficients of AR
predict_para{1} = a;

% calculate the coefficients of MA, i.e. b_i
    % calculate the residual
res = cell(1, data_diff_length - p);
for i = p+1: data_diff_length
    product_1 = data_diff{i};
    for j = 1: p
        product_1 = product_1 - a(j)*data_diff{i-j};
    end
    res{i-p} = product_1;
end
b_autocorr = zeros(1, q+1);
for l = 0: q    % calculate the auto related coefficient 
    product = 0;
    for t = 1: length(res)-l
        product = product + sum(res{t}.*conj(res{t+l}), 'all');
    end
%     b_autocorr(l+1) = product;
    b_autocorr(l+1) = product/(length(res)-l);
end
if q == 0   
    return
elseif q == 1
    b1 = (-b_autocorr(1)+sqrt(b_autocorr(1)^2-4*b_autocorr(2)^2))/(2*b_autocorr(2));
    b = [b1];
    eps = res;
    for l = 2: data_diff_length - p
        eps{l} = res{l} + b1 * eps{l-1};
    end
    predict_para{2} = b;    predict_para{3} = eps;
else
    return
end

end

