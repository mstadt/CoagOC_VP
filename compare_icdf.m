% Compare icdf after changing parameters
clear all;

% file where data is stored
fname = './results/22-Feb-2024_VP_factor_II_notes-final.mat'; 
load(fname); % load data

%% Compute default KDE and iCDF
xiNoOC_def = ksdensity(F_noOC, pi, 'Function','icdf');
xiLev_def  = ksdensity(F_lev, pi, 'Function','icdf');
xiDsg_def  = ksdensity(F_dsg, pi, 'Function','icdf');

%% Plot icdf
figure(1)
clf;
hold on
ls2 = ':';
plot(pi, xiNoOC, 'color', c_noOC, 'linewidth', lw)
plot(pi, xiLev, 'color', c_lev, 'linewidth', lw)
plot(pi, xiDsg, 'color', c_dsg, 'linewidth', lw)

plot(pi, xiNoOC_def, 'color', c_noOC, 'linewidth', lw, 'linestyle', ls2)
plot(pi, xiLev_def, 'color', c_lev, 'linewidth', lw, 'linestyle', ls2)
plot(pi, xiDsg_def, 'color', c_dsg, 'linewidth', lw, 'linestyle', ls2)
xrange = [0,1];
yrange = [20,200];
xlim(xrange)
ylim(yrange)
xlabel('p')
ylabel(sprintf('Factor %s', factor))
legend('No OC', 'Lev', 'Dsg',...
    'No OC (default)', 'Lev (default)', 'Dsg (default)',...
    'location', 'northwest')
title('Inverse CDF')
set(gca, 'fontsize', 18)



%% Compute PDF using ksdensity with same pars
[PDF_noOC2, XF_noOC2] = ksdensity(F_noOC);
[PDF_lev2, XF_lev2] = ksdensity(F_lev);
[PDF_dsg2, XF_dsg2] = ksdensity(F_dsg);

% PDF all on one axes
figure(2);
bw = 10;%5; %10;
cbin = cmap(6,:);
lw = 5;
labs = {'data', 'KDF', 'default KDF'};
xrange = [20,200];
yrange = [0,0.04];
clf;
subplot(1,3,1)
hold on;
histogram(F_noOC, 'Normalization', 'pdf', 'FaceColor', cbin, 'BinWidth',bw)
plot(XF_noOC, PDF_noOC, 'color', c_noOC, 'linewidth',lw)
plot(XF_noOC2, PDF_noOC2, 'color', c_noOC, 'linewidth',lw, 'linestyle',ls2)
xlabel(sprintf('Factor %s', factor))
ylim(yrange)
xlim(xrange)
title('No OC')
legend(labs)
set(gca, 'fontsize', 18)

subplot(1,3,2)
hold on;
histogram(F_lev, 'Normalization', 'pdf', 'FaceColor', cbin, 'BinWidth',bw)
plot(XF_lev, PDF_lev, 'color', c_lev, 'linewidth',lw)
plot(XF_lev2, PDF_lev2, 'color', c_lev, 'linewidth',lw, 'linestyle',ls2)
xlabel(sprintf('Lev Factor %s', factor))
ylim(yrange)
xlim(xrange)
title('Lev')
legend(labs)
set(gca, 'fontsize', 18)

subplot(1,3,3)
hold on;
histogram(F_dsg, 'Normalization', 'pdf', 'FaceColor', cbin, 'BinWidth',bw)
plot(XF_dsg, PDF_dsg, 'color', c_dsg, 'linewidth', lw)
plot(XF_dsg2, PDF_dsg2, 'color', c_dsg, 'linewidth', lw, 'linestyle',ls2)
xlabel(sprintf('Dsg Factor %s', factor))
ylim(yrange)
xlim(xrange)
title('Dsg')
legend(labs)
%sgtitle('Kernel density functions')
set(gca, 'fontsize', 18)

sgtitle(sprintf('Factor %s', factor), 'fontsize', 24)