function X = slice_thr_dot(A,B,C)
    % A:n1*n2*n5 , B:n2*n3*n5, C:n3*n4*n5 are three tensors£¬returnX:n1*n4*n5
    % where X(:,:,i) = A(:,:,i)*B(:,:,i)*C(:,:,i)
    n = size(A,3);
    X = zeros(size(A,1), size(C,2), n);
    for i = 1: n
        X(:,:,i) = A(:,:,i)*B(:,:,i)*C(:,:,i);
    end
end

