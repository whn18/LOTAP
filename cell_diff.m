function [X_diff,final_diff]=cell_diff(X,D) 
% X is a data sequence with 'cell' format in MATLAB, 'D' is the order of differencing
% X_dif is the d-order differencing result, 'final_dif' records the last one of each order differencing:
% final_dif{1} is the last one of 0-order differencing, final_dif{2} is the last one of 1-order differencing...
len_x=length(X);
if len_x<D+1
    error('the data sequence cannot difference with given order')
end
final_diff=cell(D+1,1);
final_diff{1}=X{end};
for d=1:D
    for j=1:len_x-d
        X{j}=X{j+1}-X{j};
    end
    final_diff{d+1}=X{len_x-d};
end
X_diff=X(1:len_x-D);
end