%% global path
clear
Path =  './';
cd(Path)
rc_wd = [Path 'rawcount.csv'];
ginfo_wd = [Path 'genomeFile/NC000913.3.gff3'];
mkdir("QC_analysis")
ncRNA = [Path 'genomeFile/ncRNA.txt'];

%% loading data
% normalize rawcount to tpm and remove counts of ssrA & rRNA
[tpm,gene,variableName,rpkm,rc_table] = rcNormalize(rc_wd, ginfo_wd, 'Remove_ncRNA',ncRNA);
if numel(variableName) >= 2
   v_tmp = split(variableName,'_'); variableName = v_tmp(:,:,1);
else
   v_tmp = split(variableName,'_'); variableName = v_tmp(1);
end
t = array2table(tpm);
t.Properties.VariableNames = variableName;
t.Properties.RowNames = gene;
rc = table2array(rc_table);
writetable(t, [Path 'TPM_rm_ncRNA.csv'],'WriteRowNames',true,'WriteVariableNames',true);

%% calculate the mean and median of tpm
tpm_mean = mean(tpm);
tpm_median = median(tpm);
t_m = array2table([tpm_mean;tpm_median]);
t_m.Properties.VariableNames = variableName;
t_m.Properties.RowNames = {'Mean','Median'};
writetable(t_m, [Path, 'QC_analysis/TPM_mean_median.txt'],'Delimiter','tab','WriteRowNames',true,'WriteVariableNames',true)

%% normalize tpm to log tpm
tpmlog = log10(tpm + 1);

%% count remaining rRNA and ssrA
nc_ratio = plotRiboRatio(rc, rc_table.Properties.RowNames, variableName, ncRNA);
saveas(gca,[Path 'QC_analysis/rRNA&ssrA_depletion_rawcount.png'])
close all

%% plot gene coverage based on rawcount
% count number of genes with at least 1 read
coverage = sum(rc >= 1)./ numel(rc_table.Properties.RowNames) *100;
figure
set(gcf,'position', [100 100 1100 650])
b = bar(coverage,'BarWidth',0.6);
b.EdgeColor = 'white';
b.FaceColor = [0.5 0.5 0.5];
set(gca,'Xtick',1:numel(variableName),'Xticklabel',variableName,'FontSize',18,'YGrid','on')
xtickangle(330)
xlabel('Library','FontSize',20)
ylabel('Coverage of genes (%)', 'FontSize',25)
ylim([0 110])
set(gca, 'yTick',[0:20:100])
box on
rc_sum = round(sum(rc)./1000);
txt = compose('%d k',rc_sum);
text(1:numel(variableName), coverage, txt,'HorizontalAlignment','center','Fontsize',25, ...
    'VerticalAlignment', 'bottom')
saveas(gca,[Path 'QC_analysis/gene_coverage.png'])
close all

%% histogram of TPM
figure
set(gcf, 'position',[100 0 400 1000])
for x = 1:numel(variableName)
    subplot(length(variableName),1,x)
    h = histogram(tpmlog(:,x));
    Counts{x} = h.Values'; 
    FPKM_lim{x} = h.BinEdges(1:end-1) + h.BinWidth/2;
    title(variableName{x},'FontSize', 15)
    ylabel('Number of genes','FontSize', 12)
    xlabel('log_1_0(TPM+1)','FontSize',12)
end
saveas(gca,[Path 'QC_analysis/histogram_log10TPM_1.png'])
% plot histogram of all clusters in 1 figure
figure
cmp = colormap(jet(numel(variableName)));
for y = 1:numel(variableName)
    hold on
    plot(FPKM_lim{y}, Counts{y}, 'color',cmp(y,:), 'LineWidth',3)
    box on
end
legend(variableName)
xlabel('log_1_0(TPM+1)','FontSize',16)
ylabel('Number of genes','FontSize', 18)
hold off
saveas(gca,[Path 'QC_analysis/histogram_log10TPM_2.png'])
% plot Cumulative distribution 
figure
hold on
for m = 1:numel(variableName)
    h = cdfplot(tpmlog(:,m));
    h.LineWidth = 3;
    h.Color = cmp(m,:);
end
legend(variableName,'FontSize', 16,'Location','southeast')
ylabel('Cumulative distribution ( F(x))','FontSize',18)
xlabel('log_1_0(TPM+1)','FontSize',20)
title('')
box on
hold off
saveas(gca,[Path 'QC_analysis/cdf_log10TPM.png'])
close all

%% plot the pearson correllation of log10(tpm+1)
correlation = corr(tpmlog,'type', 'pearson');
figure
set(gcf,'position',[150 150 900 800])
colormap(jet)
h = heatmap(variableName, variableName, correlation,'ColorScaling','scaled','GridVisible','on','ColorLimits',[0 1]);
h.FontSize = 18;
saveas(gca,[Path 'QC_analysis/pearson_correllation_log10TPM.png'])
close all

%% view the gene counts map of log10(tpm+1)
figure
set(gcf,'position', [0 0 1200 1000])
n = 1;
for r = 1:numel(variableName)
    for c = 1:numel(variableName)
        subplot(numel(variableName),numel(variableName),n);
        n = n + 1;
        scatter_kde(tpmlog(:,r),tpmlog(:,c), 'filled','MarkerSize', 5);
        set(gca,'xtick',-inf:inf:inf);
        set(gca,'ytick',-inf:inf:inf);
        box on
        if c == 1
           ylabel(variableName{r},'Rotation',0,'FontSize',15)
        end
        if r == numel(variableName)
           xlabel(variableName{c},'FontSize',15)
        end
    end
end
axis equal
saveas(gca,[Path 'QC_analysis/scatter_plot_log10TPM.png'])
close all
