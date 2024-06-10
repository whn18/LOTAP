function data_err = geterr(data_predict, data_true)
% geterr input：   两个cell数组，长度一样，cell中的元素为张量
% geterr 计算预测值与真实值的误差

data_length = length(data_predict);
data_err = zeros(1, data_length);
for i = 1: data_length
    data_err(i) = frob(data_predict{i}-data_true{i})/frob(data_true{i});
end

end

