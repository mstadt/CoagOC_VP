% Use this script to makes figures for VP results
clear all;

%% Load Data
% file where VP results are stored
fname = './results/22-Feb-2024_VP_factor_X_notes-final.mat'; 

load(fname); % load data


%% Plot data
figure(1)
if strcmp(factor,'VIII')
    yrange = [0, 5];
    xrange = [40, 200];
    w_bin = 5;
elseif strcmp(factor,'II')
    yrange = [0,8];
    xrange = [40, 200];
    w_bin = 2;
elseif strcmp(factor,'V')
    yrange = [0,12];
    xrange = [40,200];
    w_bin = 2;
elseif strcmp(factor,'VII')
    yrange = [0,10];
    xrange = [40,200];
    w_bin = 2;
elseif strcmp(factor,'X')
    yrange = [0,10];
    xrange = [50,200];
    w_bin =2;
else
    yrange = [0,28];
    xrange = [40, 200];
    w_bin = 5;
end
clf;
subplot(1,3,1)
hold on
histogram(F_noOC, ...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s before OC data', factor))
set(gca, 'fontsize', 16)

subplot(1,3,2)
histogram(F_lev, ...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s after Lev data', factor))
set(gca, 'fontsize', 16)

subplot(1,3,3)
histogram(F_dsg, ...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s after Dsg data', factor))
set(gca, 'fontsize', 16)

sgtitle(sprintf('Factor %s Data', factor))

%% Plot iCDF
lw = 3;

figure(2)
clf;
hold on
plot(pi, xiNoOC, 'color', c_noOC, 'linewidth', lw)
plot(pi, xiLev, 'color', c_lev, 'linewidth', lw)
plot(pi, xiDsg, 'color', c_dsg, 'linewidth', lw)
xrange = [0,1];
if strcmp(factor, 'VIII')
    yrange = [50,210];
elseif strcmp(factor, 'II')
    yrange = [50,200];
elseif strcmp(factor,'V')
    yrange = [50,200];
elseif strcmp(factor,'VII')
    yrange = [50,200];
elseif strcmp(factor,'X')
    yrange = [50,200];
else
    yrange = [50,200];
end
xlim(xrange)
ylim(yrange)
xlabel('p')
ylabel(sprintf('Factor %s', factor))
legend('No OC', 'Lev', 'Dsg', 'location', 'northwest')
title('Inverse CDF')
set(gca, 'fontsize', 18)

%% Plot PDF on one axis
figure(3);
if strcmp(factor,'VIII')
    xrange = [30, 200];
    yrange = [0,0.03];
elseif strcmp(factor, 'II')
    xrange = [30, 200];
    yrange = [0,0.04];
elseif strcmp(factor,'V')
    xrange = [40,200];
    yrange = [0,0.03];
elseif strcmp(factor,'VII')
    xrange = [30,200];
    yrange = [0,0.03];
elseif strcmp(factor, 'X')
    xrange = [30,200];
    yrange = [0,0.03];
else
    xrange = [30, 200];
    yrange = [0,0.03];
end
clf;
hold on;
plot(XF_noOC, PDF_noOC, 'color', c_noOC, 'linewidth',lw)
plot(XF_lev, PDF_lev, 'color', c_lev, 'linewidth',lw)
plot(XF_dsg, PDF_dsg, 'color', c_dsg, 'linewidth', lw)
xlabel(sprintf('Factor %s', factor))
ylim(yrange)
xlim(xrange)
legend('No OC', 'Lev', 'Dsg')
title('Kernel density functions')
set(gca, 'fontsize', 18)

%% Plot KDF and samples
figure(4);

if strcmp(factor,'VIII')
    xrange = [40, 220];
    yrange = [0,0.02];
    w_bin2 = 1;
elseif strcmp(factor,'II')
    xrange = [40,220];
    yrange = [0.0,0.035];
    w_bin2 = 1;
elseif strcmp(factor,'V')
    xrange = [40,220];
    yrange = [0,0.03];
    w_bin2 = 1;
elseif strcmp(factor,'VII')
    xrange = [50,200];
    yrange = [0,0.03];
    w_bin2 = 4;
elseif strcmp(factor, 'X')
    xrange = [50,200];
    yrange = [0,0.03];
    w_bin2 = 1;
else
    xrange = [40,220];
    yrange = [0,0.04];
    w_bin2 = 1;
end
lw = 4.5;
cmap2 = parula(6);
c_samp = cmap2(5,:);
c_kdf = 'red'; %cmap2(2,:);

clf;
subplot(1,3,1)
hold on;
histogram(samplesNoOC,'Normalization','pdf', ...
                'BinWidth', w_bin2, 'FaceColor', c_samp)
plot(XF_noOC, PDF_noOC, 'color', c_kdf, 'linewidth',lw)
xlabel(sprintf('Factor %s before OC', factor))
xlim(xrange)
ylim(yrange)
set(gca,'fontsize',18)

subplot(1,3,2)
hold on
histogram(samplesLev, 'Normalization', 'pdf',...
                'BinWidth', w_bin2, 'FaceColor', c_samp)
plot(XF_lev, PDF_lev, 'color', c_kdf, 'linewidth',lw)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s after Lev', factor))
set(gca,'fontsize',18)

subplot(1,3,3)
hold on
histogram(samplesDsg, 'Normalization', 'pdf', ...
                'BinWidth', w_bin2, 'FaceColor', c_samp)
plot(XF_dsg, PDF_dsg, 'color', c_kdf, 'linewidth', lw)
xlabel(sprintf('Factor %s after Dsg', factor))
xlim(xrange)
ylim(yrange)
set(gca,'fontsize',18)

legend({'VP', 'Kernel Density function'})
sgtitle({'VP and Kernel Density Function',sprintf('Factor %s', factor)})

%% Plot pairs
figure(5)
if strcmp(factor,'VIII')
    xrange = [40, 220];
elseif strcmp(factor,'II')
    xrange = [50,200];
elseif strcmp(factor,'V')
    xrange = [50,200];
elseif strcmp(factor,'VII')
    xrange = [50,200];
elseif strcmp(factor, 'X')
    xrange = [50,200];
else
    xrange = [50,200];
end
yrange = xrange;
clf;
hold on
x = linspace(xrange(1),xrange(2), 100);
y = x;
temp = gray(3);
cgray = temp(2,:);
plot(x,y, 'linewidth',1.0, 'color',cgray)
plot(samplesNoOC, samplesLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_lev)
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_dsg)
xlabel(sprintf('Factor %s before OC',factor))
ylabel(sprintf('Factor %s after OC',factor))
title({'VP pairs', ['Factor ', factor]})
legend('','Lev','Dsg') 
ylim(yrange)
xlim(xrange)
hold off
set(gca,'fontsize',18)

%% Plot probabilities pairs
figure(6)
clf;
hold on
x = linspace(0,1, 100);
y = x;
temp = gray(3);
cgray = temp(2,:);
plot(x,y, 'linewidth',4.0, 'color',cgray,'HandleVisibility','off')
plot(xqNoOC, xqLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_lev)
plot(xqNoOC, xqDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_dsg)
xlabel('xq before OC')
ylabel('xq after OC')
title({'probability pairs', ['Factor ', factor]})
legend('Lev','Dsg') 
ylim([0,1])
xlim([0,1])
set(gca,'fontsize',18)
hold off

%% Plot pair differences
figure(7)
clf;
hold on
diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;
if strcmp(factor,'VIII')
    yrange = [0,0.03];
elseif strcmp(factor, 'II')
    yrange = [0,0.12];
elseif strcmp(factor,'V')
    yrange = [0,0.08];
elseif strcmp(factor,'VII')
    yrange = [0,0.06];
else
    yrange = [0,0.1];
end
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', c_lev, ...
                'Normalization', 'pdf')
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', c_dsg,...
                'Normalization', 'pdf')

xlabel(sprintf('Factor %s level difference after OC',factor))
ylabel('density')
legend('Lev', 'Dsg')
ylim(yrange)
set(gca,'fontsize',18)
hold off

%% Plot xqNoOC, xqLev, xqDsg
figure(8)
wbin = 0.05;
clf;
subplot(1,2,1)
hold on
histogram(xqNoOC,...
            'BinWidth', wbin, 'FaceColor', c_noOC)
histogram(xqLev,...
            'BinWidth', wbin, 'FaceColor',c_lev)
xlabel('p')
ylabel('count')
legend('no oc', 'lev')

set(gca,'fontsize',18)
hold off

subplot(1,2,2)
hold on
histogram(xqNoOC,...
            'BinWidth', wbin, 'FaceColor', c_noOC)
histogram(xqDsg,...
            'BinWidth', wbin, 'FaceColor',c_dsg)
xlabel('p')
ylabel('count')
legend('no oc', 'dsg')

set(gca,'fontsize',18)
hold off


sgtitle(sprintf('Factor %s', factor), 'fontsize', 20)