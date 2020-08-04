% demo_hpt() - A demo of 'Hyperbolic power transformation'
%
function demo_hpt

close all;
h=figure;
fig_height=1600;
set(gcf,'PaperPositionMode','auto');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), fig_height./2, fig_height]);

transformation('lognormal',h,1);
transformation('Gamma',h,2);
transformation('beta',h,3);
transformation('Laplacian',h,4);
transformation('uniform',h,5);
transformation('bimodal',h,6);

function transformation(distfamily,h,j)
figure(h);
M=6; % number of rows to plot
FontSize=10;

load([distfamily '.mat'],'x');

% ==== 1. original values =====================================================
subplot(M,3,3*j-2);
hist2(x,20);
if j==M, xlabel('original value'); end;

txtline0=distfamily;
txtline1=['skewness=' num2str(skewness(x),1), ', kurtosis=' num2str(kurtosis(x),2), ', jb=' num2str(jbtest(x),1)];
title(sprintf('%s\n%s',txtline0,txtline1),'Units', 'normalized','Position', [0 1], 'HorizontalAlignment', 'left');

% ==== 2. transformed values by method of percentile  ===========================
subplot(M,3,3*j-1);
[a, betaminus, lambdaminus, betaplus, lambdaplus, x] = hyperdistmop(x);
y=zeros(length(x),1);
y((x>=0))=(a./betaplus).*sinh(betaplus.*x(x>=0)) .* ( sech(betaplus.*x(x>=0)) ).^lambdaplus;
y((x<0))=(a./betaminus).*sinh( betaminus.*x(x<0)) .* ( sech(betaminus.*x(x<0)) ).^lambdaminus;

hist2(y, 20); hold on;

if j>1, txtline0=''; else, txtline0=['By Method of Percentiles.']; end;
txtline1=['skewness=' num2str(skewness(y),1), ', kurtosis=' num2str(kurtosis(y),2)];
txtline2=['\alpha=' num2str(a,2)...
' \beta_-=' num2str(betaminus,2) ' \lambda_-=' num2str(lambdaminus,2)...
' \beta_+=' num2str(betaplus,2)  ' \lambda_+=' num2str(lambdaplus,2)];
title(sprintf('%s\n%s',txtline0,txtline1),'Units', 'normalized','Position', [0 1], 'HorizontalAlignment', 'left');
xlabel(txtline2); set(gca,'xlim',[-4 4]); set(gca,'ylim',[0 1]);
h2 = pG(y, 'b-', 0); % plot an standard normal distribution.

% ==== 3. transformed values by MLE ===========================================
subplot(M,3,3*j);
[a, betaminus, lambdaminus, betaplus, lambdaplus] = hyperdistfminsearch(x, a, betaminus, lambdaminus, betaplus, lambdaplus);
y=zeros(length(x),1);
y((x>=0))=(a./betaplus).*sinh(betaplus.*x(x>=0)) .* ( sech(betaplus.*x(x>=0)) ).^lambdaplus;
y((x<0))=(a./betaminus).*sinh( betaminus.*x(x<0)) .* ( sech(betaminus.*x(x<0)) ).^lambdaminus;

[h1, count]=hist2(y, 20); hold on;

if j>1, txtline0=''; else, txtline0=['By MLE.']; end;
txtline1=['skewness=' num2str(skewness(y),1), ', kurtosis=' num2str(kurtosis(y),2), ', jb=' num2str(jbtest(y),1)];
txtline2=['\alpha=' num2str(a,2)...
' \beta_-=' num2str(betaminus,2) ' \lambda_-=' num2str(lambdaminus,2)...
' \beta_+=' num2str(betaplus,2)  ' \lambda_+=' num2str(lambdaplus,2)];
title(sprintf('%s\n%s',txtline0,txtline1),'Units', 'normalized','Position', [0 1], 'HorizontalAlignment', 'left');
xlabel(txtline2); set(gca,'xlim',[-4 4]); set(gca,'ylim',[0 1]);
hold on;
h2 = pG(y, 'b-', 0); % plot an standard normal distribution.

%DagosPtest(y); % D'Agostino-Pearson's test to assessing normality
fig_height=1600;
set(gcf,'PaperPositionMode','auto');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), fig_height, fig_height./1.8]);

print('-djpeg', ['demp_hpt.jpg']);
%print('-depsc', ['demp_hpt.eps']); % use -depsc instead of -deps, so you can have color with gray level

function [h counts] = hist2(data, default_counts)
    binWidth = (max(data)-min(data))/default_counts;
    binCtrs = min(data):binWidth:max(data);
    counts = hist(data,binCtrs); 
    n = length(data);
    prob = counts / (n * binWidth);
    bar(binCtrs,prob,'hist');
    h = get(gca,'child');
    set(h(1),'FaceColor',[.9 .9 .9]);

