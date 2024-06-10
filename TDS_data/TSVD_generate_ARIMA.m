function [data,fft_S,fft_U,fft_V]=TSVD_generate_ARIMA(T,scale_data,scale_core,AR_evolution,d)
% -----------TDS: abbr. of Tensor Data Sequence--------
% T                 - Length of Sequence
% scale_data        - scale of each tensor data
% scale_core        - scale of each tensor core, number
% AR_evolution      - all AR parameters with evolution from time max(length(alpha),length(beta))+1 
%                         to the end in a 'matrix' format, each row corresponds to the AR parameters in a time
% MA_evolution      - all MA parameters with evolution from time max(length(alpha),length(beta))+1 
%                         to the end in a 'matrix' format, each row corresponds to the MA parameters in a time
% sigma             - square root of invariance of i.i.d. random noises
% d                 - order of differencing
p=size(AR_evolution,2);
n1 = scale_data(1);
n2 = scale_data(2);
n3 = scale_data(3);
s=p+d;
M=3;

%%  Error hint
if length(scale_data)~=M
    error('The dimension of tensor data is not equal to 3')
end
if min(scale_data(1:2)-scale_core)<0
    error('The scale of tensor core can not large than the scale of data')
end
if T<=s
    error('The length of data sequence is not more than the regression length')
end

%%  Initialization
fft_data=cell(T,1);     data = fft_data;
S=cell(T,1);
fft_S = S;

% generate real tensor U:n1*r*n3, V:n2*r*n3
A = randn(n1,n2,n3);
fft_A = fft(A,[],3);
fft_U = zeros(n1,scale_core,n3);  fft_V = zeros(n2,scale_core,n3);
for i = 1: n3
    [U, ~, V] = svd(fft_A(:,:,i),'econ');
    fft_U(:,1:scale_core,i) = U(:,1:scale_core);
    fft_V(:,1:scale_core,i) = V(:,1:scale_core);
end

for t=1:s
    S{t} = 1000*rand(scale_core, n3);
    fft_S{t} = fft(S{t},[],2);
    % fft_S{t} = 1000*rand(scale_core, n3);
end

TDS_dif=cell(T-d,1);
[TDS_dif(1:s-d),TDS_final_dif]=cell_diff(fft_S(1:s),d);
%% Data Generation with ARIMA(p,d,q) model
for t=(s+1):T
    TDS_dif{t-d}=weighted_sum(AR_evolution(t-s,:),TDS_dif(t-d-p:t-d-1));
    TDS_final_dif=update_final_dif(TDS_final_dif,TDS_dif{t-d});
    fft_S{t}=TDS_final_dif{1};
end
for t = 1: T
    fft_data{t} = zeros(scale_data);
    GG = zeros(scale_core,scale_core,n3);
    for i = 1: n3
        GG(:,:,i) = diag(fft_S{t}(:,i));
        fft_data{t}(:,:,i) = fft_U(:,:,i) * GG(:,:,i) * fft_V(:,:,i)';
    end
    fft_S{t} = GG;
    data{t} = ifft(fft_data{t}, [], 3);
end


end


%% subfunction for generating core tensor
function result=weighted_sum(a,b) 
% a is a list, b is a cell, len(a)=len(b), 
% output: summation with reverse order
    result=zeros(size(b{1}));
    for i=1:length(a)
        result=result+a(i)*b{end+1-i};
    end
end