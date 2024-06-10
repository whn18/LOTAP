%% load data
clear
load tds;

%% write in .h5 file
name = "/data";
h5create("TDS.h5",name,[1000, 100, 100])
data = zeros(100, 100, 1000);
for t = 1: 100
    data(t, :, :) = reshape(X{t}(:,:,:), 100, 1000);
end
data = permute(data,[3 2 1]);
h5write("TDS.h5","/data",data)
h5disp("TDS.h5")