% Generate N VP for Factor X

% clear
clear all;

% set factor
factor = 'X';
note = 'final';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Data
temp = strcat('data/Factor', factor, '_noOC1.csv');
csvtable1 =  readtable(temp);
temp = strcat('data/Factor', factor, '_lev.csv');
csvtable2 =  readtable(temp);
temp = strcat('data/Factor', factor, '_noOC2.csv');
csvtable3 =  readtable(temp);
temp = strcat('data/Factor', factor, '_dsg.csv');
csvtable4 =  readtable(temp);
F_noOC1 = sort(csvtable1.Var2); 
F_lev   = sort(csvtable2.Var2);
F_noOC2 = sort(csvtable3.Var2);
F_dsg   = sort(csvtable4.Var2);

F_noOC = sort([F_noOC1; F_noOC2]); % merge noOC data together

clearvars -except F_noOC F_lev F_dsg factor note;

%Desired Mean Difference After Treatment (Lev - NoOC):
% Factor X
MEAN_lev = 22; %Table 1 taken from Middeldorp et al. 2000
MEAN_lev_range = [max(MEAN_lev - 0.5, 0.95*MEAN_lev),...
                        min(MEAN_lev +  0.5, 1.05*MEAN_lev)]; 
STD_lev  = 14;  %The Standard Deviation.
STD_lev_range = [max(STD_lev - 0.5, 0.95 * STD_lev),...
                        min(STD_lev + 0.5, 1.05*STD_lev);];

%Desired Mean Difference After Treatment (Dsg - NoOC):
MEAN_dsg = 25;
MEAN_dsg_range = [max(0.95 * MEAN_dsg, MEAN_dsg - 0.5), ...
                        min(1.05*MEAN_dsg, MEAN_dsg + 0.5)];
STD_dsg  = 12;
STD_dsg_range = [max(0.95 * STD_dsg,STD_dsg - 0.5),...
                        min(STD_dsg + 0.5, 1.05 * STD_dsg)];

% Hyperparameters for VP pairing
sigma_lev = 0.25;
sigma_dsg = 0.2;

N_vp = 1e4;

MAX_TRIALS = 25;

% set random seed
rng(25)

%% Plot Data
cmap = parula(7);
c_noOC = cmap(1,:); c_lev = cmap(3,:); c_dsg = cmap(5,:);
c_data = cmap(7,:);
figure(3)
w_bin = 5;
xrange = [40, 200];
yrange = [0, 0.06];
clf;
subplot(1,3,1)
hold on
histogram(F_noOC, 'Normalization','pdf',...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s before OC data', factor))
set(gca, 'fontsize', 16)

subplot(1,3,2)
histogram(F_lev, 'Normalization', 'pdf',...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s after Lev data', factor))
set(gca, 'fontsize', 16)

subplot(1,3,3)
histogram(F_dsg, 'Normalization', 'pdf',...
                    'BinWidth', w_bin, 'FaceColor', c_data)
xlim(xrange)
ylim(yrange)
xlabel(sprintf('Factor %s after Dsg data', factor))
set(gca, 'fontsize', 16)

sgtitle(sprintf('Factor %s Data', factor))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute kernel density estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel density parameters
fac_range = {[50, 175], [60, 175], [70, 175]}; % range allowed for factor values
bw = [0.35, 0.5, 0.50]; % bandwidth values for noOC, lev, dsg
kernels = {'normal', 'normal', 'normal'}; % kernels for noOC, lev, dsg

% weights for data
weight_noOC = ones(size(F_noOC));
p = 0.1; 
for ii = 1:length(F_noOC)
    noOC = F_noOC(ii);
    diff_mean = abs(noOC - mean(F_noOC));
    if diff_mean < p * mean(F_noOC) % weight values within p% of the mean
        weight_noOC(ii) = 1.2;
    end

end

weight_lev = ones(size(F_lev));
p = 0.1;
for ii = 1:length(F_lev)
    lev = F_lev(ii);
    diff_mean = abs(lev - mean(F_lev));
    if diff_mean < p * mean(F_lev) % weight values within p% of the mean
        weight_lev(ii) = 1.15;
    end

%     if lev < 100
%         weight_lev(ii) = 0.9;
%     end
end

weight_dsg = ones(size(F_dsg));
p = 0.1;
for ii = 1:length(F_dsg)
    dsg = F_dsg(ii);
    diff_mean = abs(dsg - mean(F_dsg));
    if diff_mean < p * mean(F_dsg) % weight values within p% of the mean
        if dsg > mean(F_dsg)
            weight_dsg(ii) = 1.3;
        else
            weight_dsg(ii) = 1.2;
        end
    end

    if dsg < 100
        weight_dsg(ii) = 0.9;
    end
    
end
weights = {weight_noOC, weight_lev, weight_dsg}; % weigths for kernel density
figure(5)
xrange = [50,200];
clf;
subplot(3,1,1)
plot(F_noOC, weight_noOC, 'marker', '.', 'markersize', 20)
ylabel('weight noOC')
xlabel('data')
xlim(xrange)
grid on
set(gca, 'fontsize', 16)
subplot(3,1,2)
plot(F_lev, weight_lev, 'marker', '.', 'markersize', 20)
ylabel('weight lev')
xlabel('data')
xlim(xrange)
grid on
set(gca, 'fontsize', 16)
subplot(3,1,3)
plot(F_dsg, weight_dsg, 'marker', '.', 'markersize', 20)
ylabel('weight dsg')
xlabel('data')
xlim(xrange)
grid on
set(gca, 'fontsize', 16)

sgtitle('KDF weights')

% Compute inverse kernel density functions
%Range 
pi = linspace(.0001,.9999,10000);
xiNoOC = ksdensity(F_noOC, pi, 'Function', 'icdf',...
                            'Support', fac_range{1},... % range of factor level
                            'Bandwidth', bw(1),... % bandwidth
                            'Kernel', kernels{1},... % kernel
                            'Weights', weights{1}... % weights
                            );
xiLev = ksdensity(F_lev, pi, 'Function', 'icdf',...
                            'Support', fac_range{2},... % range of factor level
                            'Bandwidth', bw(2),... % bandwidth
                            'Kernel', kernels{2},... % kernel
                            'Weights', weights{2}... % weights
                            );
xiDsg = ksdensity(F_dsg, pi, 'Function', 'icdf',...
                            'Support', fac_range{3},... % range of factor level
                            'Bandwidth', bw(3),... % bandwidth
                            'Kernel', kernels{3},... % kernel
                            'Weights', weights{3}... % weights
                            );

%% plot icdf
lw = 3;

figure(1)
clf;
hold on
plot(pi, xiNoOC, 'color', c_noOC, 'linewidth', lw)
plot(pi, xiLev, 'color', c_lev, 'linewidth', lw)
plot(pi, xiDsg, 'color', c_dsg, 'linewidth', lw)
xrange = [0,1];
yrange = [50,180];
xlim(xrange)
ylim(yrange)
xlabel('p')
ylabel(sprintf('Factor %s', factor))
legend('No OC', 'Lev', 'Dsg', 'location', 'northwest')
title('Inverse CDF')
set(gca, 'fontsize', 18)

%% Compute PDF using ksdensity with same pars
[PDF_noOC, XF_noOC] = ksdensity(F_noOC, ...
                            'Support', fac_range{1},... % range of factor level
                            'Bandwidth', bw(1),... % bandwidth
                            'Kernel', kernels{1},... % kernel
                            'Weights', weights{1}... % weights
                            );
[PDF_lev, XF_lev] = ksdensity(F_lev, ...
                            'Support', fac_range{2},... % range of factor level
                            'Bandwidth', bw(2),... % bandwidth
                            'Kernel', kernels{2},... % kernel
                            'Weights', weights{2}... % weights
                            );
[PDF_dsg, XF_dsg] = ksdensity(F_dsg, ...
                            'Support', fac_range{3},... % range of factor level
                            'Bandwidth', bw(3),... % bandwidth
                            'Kernel', kernels{3},... % kernel
                            'Weights', weights{3}... % weights
                            );
% PDF all on one axes
figure(2);
xrange = [30, 200];
yrange = [0,0.03];
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

%% Sample VPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: find NoOC, Lev, and Dsg samples that within tolerance of 
%    mean difference and std difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OBJ_mean = 0;
OBJ = 0; % Obj: both OBJ_mean and OBJ_std are met

NUM_TRIALS = 0;

fprintf('start algorithm \n')

% Figure for plotting MEAN and STD differences
figure(4);
clf;
subplot(2,2,1)
hold on
yline(MEAN_lev_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--', ...
                        'HandleVisibility', 'off')
yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-',...
                            'HandleVisibility', 'off')
yline(MEAN_lev_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--',...
                            'HandleVisibility', 'off')
xlabel('trial')
ylabel('mean difference (Lev - NoOC)')
title('MEAN difference')
set(gca, 'fontsize', 18)

subplot(2,2,2)
hold on
yline(MEAN_dsg_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--',...
                            'HandleVisibility', 'off')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-',...
                            'HandleVisibility', 'off')
yline(MEAN_dsg_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--',...
                            'HandleVisibility', 'off')
xlabel('trial')
ylabel('mean difference (Dsg - NoOC)')
title('MEAN difference')
set(gca, 'fontsize', 18)

subplot(2,2,3)
hold on
yline(STD_lev_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_lev_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference (Lev - NoOC)')
title('STD difference')
set(gca, 'fontsize', 18)

subplot(2,2,4)
hold on
yline(STD_dsg_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_dsg_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('STD difference (Dsg - NoOC)')
title('STD difference')
set(gca, 'fontsize', 18)

while ~OBJ
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS > MAX_TRIALS
        fprintf('MAX TRIALS reached. \n')
        break;
    end

    % Sample probabilities for No OC
    xqNoOC = rand(1,N_vp); % sample probabilities

    % Interpolate corresponding factor levels for xqNoOC using
    %   inverse CDF
    samplesNoOC = interp1(pi, xiNoOC, xqNoOC);

    % resamples xqNoOC if nan values
    while sum(isnan(samplesNoOC)) > 0
        % resample nan values
        nan_ids = find(isnan(samplesNoOC));
        if length(nan_ids) > 10
            error('nan ids longer than 10 \n')
        end

        for ii = 1:length(nan_ids)
            id = nan_ids(ii);
            xqNoOC(id) = rand(1); % new value
        end
        samplesNoOC = interp1(pi,xiNoOC,xqNoOC);
    end % while sum(isnan)>0
    
    % Add random noise to xqNoOC for xqLev
    DELTA_lev = sigma_lev * randn(size(xqNoOC)); % random noise (normal dist)
    xqLev = xqNoOC + DELTA_lev; % add noise to xqNoOC

    % Resample xqLev that are out of range
    inds = find(xqLev > max(pi) | xqLev < min(pi)); % out of range indices
    for ii = 1:length(inds)
        ind = inds(ii);
        xqNoOC_val = xqNoOC(ind);

        % Redraw from truncated probability distribution
        pd = makedist('Normal', 'mu', 0, 'sigma', sigma_lev); % normal distribution
        % truncate distribution
        max_inc = max(pi) - xqNoOC_val; % maximum increase
        max_dec = min(pi) - xqNoOC_val; % maximum decrease
        t = truncate(pd, max_dec, max_inc);

        DELTA_new = random(t, 1,1); % get new noise value
        xqLev_val_new = xqNoOC_val + DELTA_new;

        if or(xqLev_val_new > max(pi), xqLev_val_new < min(pi))
            fprintf('xqlev: %0.2f',xqLev_val_new)
            error('xqlev out of range!')
        end

        DELTA_lev(ind) = DELTA_new;
        xqLev(ind) = xqLev_val_new;
    end % for ii

    % Add random noise to xqNoOC for xqDsg
    DELTA_dsg = sigma_dsg * randn(size(xqNoOC)); % random noise (normal dist)
    xqDsg = xqNoOC + DELTA_dsg;

    % Resample xqDsg that are out of range
    inds = find(xqDsg > max(pi) | xqDsg < min(pi)); % out of range indices
    for ii = 1:length(inds)
        ind = inds(ii);
        xqNoOC_val = xqNoOC(ind);

        % Redraw from truncated probability distribution
        pd = makedist('Normal', 'mu', 0, 'sigma', sigma_dsg); % normal distribution
        
        max_inc = max(pi) - xqNoOC_val; % maximum increase
        max_dec = min(pi) - xqNoOC_val; % maximum decrease
        t = truncate(pd, max_dec, max_inc); % truncate distribution

        DELTA_new = random(t,1,1);% get new noise value
        xqDsg_val_new = xqNoOC_val + DELTA_new;

        if or(xqDsg_val_new > max(pi), xqDsg_val_new < min(pi))
            fprintf('xqdsg: %0.2f',xqDsg_val_new)
            error('xqdsg out of range!')
        end

        DELTA_dsg(ind) = DELTA_new;
        xqDsg(ind) = xqDsg_val_new;
    end % for ii

    % Interpolate LEV and DSG values from corresponding icdf
    samplesLev = interp1(pi,xiLev,xqLev);
    samplesDsg = interp1(pi,xiDsg,xqDsg);

    % Compute mean and std difference
    diff_lev = samplesLev - samplesNoOC;
    mean_diff_lev = mean(diff_lev);
    std_diff_lev = std(diff_lev);

    diff_dsg = samplesDsg - samplesNoOC;
    mean_diff_dsg = mean(diff_dsg);
    std_diff_dsg = std(diff_dsg);


    % Plot mean and std difference for trial
    figure(4);
    subplot(2,2,1)
    hold on
    plot(NUM_TRIALS,mean_diff_lev, 'linestyle', 'none', 'marker', 'square',...
                            'markersize', 10,...
                            'MarkerFaceColor', c_lev, 'color', c_lev)

    subplot(2,2,2)
    hold on
    plot(NUM_TRIALS,mean_diff_dsg, 'linestyle', 'none', 'marker', 'square',...
                            'markersize', 10,...
                            'MarkerFaceColor', c_dsg, 'color', c_dsg)
    subplot(2,2,3)
    hold on
    plot(NUM_TRIALS, std_diff_lev, 'linestyle','none','marker','o',...
                    'markersize', 10,...
                    'MarkerFaceColor', c_lev, 'color', c_lev)
    subplot(2,2,4)
    hold on
    plot(NUM_TRIALS, std_diff_dsg, 'linestyle', 'none', 'marker', 'o',...
                    'markersize', 10,...
                    'markerfacecolor', c_dsg, 'color', c_dsg)
    
    % Check objectives
    OBJ_mean = 0; OBJ_std = 0;

    % MEAN difference objective
    OBJ_lev = 0; OBJ_dsg = 0;
    if and(MEAN_lev_range(1) < mean_diff_lev , mean_diff_lev < MEAN_lev_range(2))
        OBJ_lev = 1;
    end
    if and(MEAN_dsg_range(1) < mean_diff_dsg , mean_diff_dsg < MEAN_dsg_range(2))
        OBJ_dsg = 1;
    end

    if and(OBJ_lev, OBJ_dsg)
        OBJ_mean = 1;
    end

    % STD difference objective
    OBJ_lev = 0; OBJ_dsg = 0;
    if and(STD_lev_range(1) < std_diff_lev, std_diff_lev < STD_lev_range(2))
        OBJ_lev = 1;
    end
    if and(STD_dsg_range(1) < std_diff_dsg, std_diff_dsg < STD_dsg_range(2))
        OBJ_dsg = 1;
    end
    if and(OBJ_lev, OBJ_dsg)
        OBJ_std = 1;
    end

    % OBJ
    if and(OBJ_mean, OBJ_std)
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
        fprintf('Factor %s samples complete! \n', factor)
    end
end % ~OBJ




%% Plot results
fprintf('plotting results \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples on kernel density function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
w_bin2 = 2;
xrange = [40, 190];
yrange = [0,0.03];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrange = [50, 200];
yrange = xrange;
figure(8)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot probabilities pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
clf;
hold on
plot(xqNoOC, xqLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_lev)
plot(xqNoOC, xqDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', c_dsg)
xlabel(sprintf('xq before OC',factor))
ylabel(sprintf('xq after OC',factor))
title({'probability pairs', ['Factor ', factor]})
legend('Lev','Dsg') 
ylim([0,1])
xlim([0,1])
set(gca,'fontsize',18)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot differences for pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
clf; 
hold on
diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;

yrange = [0,0.05];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xqNoOC, xqLev, xqDsg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
wbin = 0.05;
clf;
hold on
histogram(xqNoOC,...
            'BinWidth', wbin, 'FaceColor', c_noOC)
histogram(xqLev,...
            'BinWidth', wbin, 'FaceColor',c_lev)
histogram(xqDsg,...
            'BinWidth', wbin, 'FaceColor', c_dsg)
xlabel('p')
ylabel('count')
legend('no oc', 'lev', 'dsg')
title(sprintf('Factor %s', factor))
set(gca,'fontsize',18)
hold off

figure(14)
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
