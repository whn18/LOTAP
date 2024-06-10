function data_err = geterr(data_predict, data_true)
% geterr input��   ����cell���飬����һ����cell�е�Ԫ��Ϊ����
% geterr ����Ԥ��ֵ����ʵֵ�����

data_length = length(data_predict);
data_err = zeros(1, data_length);
for i = 1: data_length
    data_err(i) = frob(data_predict{i}-data_true{i})/frob(data_true{i});
end

end

