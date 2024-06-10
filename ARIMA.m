function predict_data = ARIMA(data, para, predict_para)
% ARIMA model : x_t = a_1*x_{t-1}+...+a_p*x_{t-p}-b_1*e_{t-1}-...-b_q*e_{t-q}+e_t
% input:    data = {G_1, ..., G_n} is a cell vector of tensor.
%           para = [p d q] is the order of the ARIMA model.
%           predict_para = {[a_1, a_2, ..., a_p], [b_1, ..., b_q], {e_1,...,e_T}} 
%                           where a_i, b_j are the coefficients of AR and MA,
%                           T = (length(data_diff) - p), e_i is the residual of
%                           model.
% output:   if data.length < para(1)+para(2), predict_data = -999
%           else, predict_data is the prediction calculated by ar model

p = para(1);    d = para(2);    q = para(3);
a = predict_para{1};    b = predict_para{2};    eps = predict_para{3};
data_length = length(data);
if data_length < d + p
    predict_data = -999;
    return
end
data_diff = data;   % data_diff is the d order difference of data
diff_begin = cell(1, d); % diff_begin is used for recovering the origin data
                         % diff_begin = {x1, x2-x1, x3-2*x2+x1, ...}
% difference-process
for i = 1: d
    diff_begin{i} = data_diff{1};
    for j = 1: data_length - i
        data_diff{j} = data_diff{j+1} - data_diff{j};
    end
end
data_diff = data_diff(1: data_length - d);
data_diff_length = length(data_diff);

% calculate the predict data of data_diff, i.e.predict_data_diff
if data_diff_length < p
    predict_data = -999;
    return
else
    predict_data_diff = zeros(size(data{1}));   % the predict_data without inverse-difference-process
    for i = 1:p
        predict_data_diff = predict_data_diff + a(i) * data_diff{end+1-i};
    end
    if q==1
        predict_data_diff = predict_data_diff - b(1)*eps{end};
    end
end

% inverse-difference-process
data_diff_length = data_diff_length + 1;
new_data_diff = cell(1, data_diff_length);
new_data_diff(1:data_diff_length-1) = data_diff;    new_data_diff{end} = predict_data_diff;
for i = 1: d
    new_data_diff{1} = new_data_diff{1} + diff_begin{d+1-i};
    for j = 2: length(new_data_diff)
        new_data_diff{j} = new_data_diff{j} + new_data_diff{j-1};
    end
    temp = cell(1, length(new_data_diff)+1);    temp{1} = diff_begin{d+1-i};    temp(2:end) = new_data_diff;
    new_data_diff = temp;
end
predict_data = new_data_diff{end};
    

end

