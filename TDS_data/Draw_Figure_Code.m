clear
clc
load time_vs_alg
%%
Time = 1:50;
alg_name={'AR_ERR', 'MCAR_ERR', 'MCAR_ERR2', 'thot_err2'};
legend_name={'ARMA', 'MCAR(small core)', 'MCAR(large core)', 'TSVD-AR'};
LINE={'-','--',':','--','-.',   '-',':','-','--','-.',   ':','--' };
MARK={'^','o','diamond','+','*',  'h','s','<','>','p',   'v','x'};


figure(1)
grid on
for i=1:length(alg_name)
    hold on
    plot(Time,eval(alg_name{i}),'Marker',MARK{i},'LineStyle',LINE{i},'linewidth',1.25)
end
xlim([Time(1)-1 Time(end)+1])
xlabel('Prediction Time','FontName','Times New Roman','FontSize',18);
ylabel('MSRE','FontName','Times New Roman','FontSize',18);
LEGEND=cell(length(alg_name),1);
for i=1:length(alg_name)
    LEGEND{i}=legend_name{i};
end
legend(LEGEND)
