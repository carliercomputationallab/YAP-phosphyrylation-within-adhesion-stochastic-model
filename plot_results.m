%##################
%Plot scripts  
%##################


%% fig 2a

clc;
clear;
close all;

load('20230908T222115_results\sspYAP.mat'); % N=1 small
results_1 = results; 
load('20230815T170105_results\sspYAP.mat'); % this is low diff rate (N=9 small , N=1 large)

figure('Position', [10 10 550 400]); 
line_style = {'-','--',':'};

plot(results(1).time{1,4}(:,1:8.5*10^6)', results(1).pYAPAll{1,4}(:,1:8.5*10^6)'/particleNumber,'LineWidth',1, 'Color','k','LineStyle',line_style{1});
set(gca,'Fontname', 'Times New Roman','fontsize',14);
hold on;

plot(results(1).time{1,1}(:,1:8.5*10^6)', results(1).pYAPAll{1,1}(:,1:8.5*10^6)'/particleNumber,'LineWidth',1, 'Color','k','LineStyle',line_style{2});
set(gca,'Fontname', 'Times New Roman','fontsize',14);
hold on;

plot(results_1(1).time{1,1}(:,1:8.5*10^6)', results_1(1).pYAPAll{1,1}(:,1:8.5*10^6)'/particleNumber,'LineWidth',1, 'Color','k','LineStyle',line_style{3});
set(gca,'Fontname', 'Times New Roman','fontsize',14);
xlabel('Time (s)'); ylabel('pYAP ratio');
hold on;

% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

axis([0 135 0 0.6])

legend('{\it{N}} = 9, small adhesion','','','{\it{N}} = 1, large adhesion','','','{\it{N}} = 1, small adhesion','','', ...
    'Position', [0.6, 0.3, 0.2, 0.2],'Box','off','FontSize', 12)

folder = 'figures\';
print([folder 'fig2a'],'-dpng','-r300');

%% fig 2c

clc
clear
close all

%diffusion_rate = [19 0.8];

load('20230815T170105_results\sspYAP.mat'); % this is low diff rate
copy_results1 = results;

load('20230823T005502_results\sspYAP.mat'); % this is 25% area
copy_results2 = results;

load('20230823T120041_results\sspYAP.mat'); % this 9% area


figure('Position', [10 10 350 500]);
name = string();
% name(1) = ['{\it{N}} = 1']; name(2) = ['{\it{N}} = 9,  r = 0.09'];
% name(3) = ['{\it{N}} = 9, r = 0.25']; name(4) = ['{\it{N}} = 9, r = 1'];
name(1) = 1; name(2) = 2;
name(3) = 3; name(4) = 4;

data_1 = [copy_results1(1).avg_iterations_pYAPAll{1,1};
          results(1).avg_iterations_pYAPAll{1,1};
          copy_results2(1).avg_iterations_pYAPAll{1,1};
          copy_results1(1).avg_iterations_pYAPAll{1,4}]/particleNumber;

data_1_error = [copy_results1(1).std_iterations_pYAPAll{1,1};
          results(1).std_iterations_pYAPAll{1,1};
          copy_results2(1).std_iterations_pYAPAll{1,1};
          copy_results1(1).std_iterations_pYAPAll{1,4}]/particleNumber;

%b = bar(data_1,'grouped','FaceColor', 'flat'); hold on;
b = bar(data_1,'grouped'); hold on;

bar_color = gray(4);
b(1).FaceColor = bar_color(3,:);

[ngroups,nbars] = size(data_1);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',data_1,data_1_error,'k','linestyle','none');
hold off;

ylabel('pYAP ratio');
xlabel('Case');

set(gca,'Fontname', 'Times New Roman','fontsize',22);
set(gca,'LineWidth',1);
xticklabels(name);
ax = gca;
ax.XAxis.FontSize = 18;

b.CData = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560];


% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

ylim([0.2, 0.55]);
% Set ticks every 0.1 on the y-axis
yticks(0.2:0.1:0.5);

folder = 'figures\';
print([folder 'fig2c'],'-dpng','-r300');

%% fig 2d

clc
clear
close all

diffusion_rate = [19 0.8];

load('20230815T170603_results\sspYAP.mat'); %this is high diff rate
copy_results = results;

load('20230815T170105_results\sspYAP.mat'); % this is low diff rate

figure('Position', [10 10 500 350]);
name = string();
name(1) = ['{\it{N}} = 1']; name(2) = ['{\it{N}} = 2'];
name(3) = ['{\it{N}} = 4']; name(4) = ['{\it{N}} = 9'];

data_1 = [results(1).avg_iterations_pYAPAll{1,1} copy_results(1).avg_iterations_pYAPAll{1,1}; ...
          results(1).avg_iterations_pYAPAll{1,2} copy_results(1).avg_iterations_pYAPAll{1,2}; ...
          results(1).avg_iterations_pYAPAll{1,3} copy_results(1).avg_iterations_pYAPAll{1,3}; ...
          results(1).avg_iterations_pYAPAll{1,4} copy_results(1).avg_iterations_pYAPAll{1,4}]/particleNumber;


data_1_error = [results(1).std_iterations_pYAPAll{1,1} copy_results(1).std_iterations_pYAPAll{1,1}; ...
                results(1).std_iterations_pYAPAll{1,2} copy_results(1).std_iterations_pYAPAll{1,2}; ...
                results(1).std_iterations_pYAPAll{1,3} copy_results(1).std_iterations_pYAPAll{1,3}; ...
                results(1).std_iterations_pYAPAll{1,4} copy_results(1).std_iterations_pYAPAll{1,4}]/particleNumber;

bar_color = gray(4);
b = bar(data_1,'grouped'); hold on;
b(1).FaceColor = bar_color(2,:);
b(2).FaceColor = bar_color(3,:);

[ngroups,nbars] = size(data_1);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',data_1,data_1_error,'k','linestyle','none');
hold off;

ylabel('pYAP ratio');
xticklabels(name);

legend(' {\it{D}} = ' + string(flip(diffusion_rate)) + ' {\mu}m^{2}/s','location','northoutside','orientation', 'horizontal','Box','off');

set(gca,'Fontname', 'Times New Roman','fontsize',18);
set(gca,'LineWidth',1);

% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

ylim([0.2, 0.8]);

folder = 'figures\';
print([folder 'fig2d'],'-dpng','-r300');



%% fig 3a 

clc;
clear;
close all;

load('20230819T164845_results\sspYAP.mat');
results_copy = results;

load('20230816T142738_results\sspYAP.mat');

results.avg_iterations_pYAPAll(4, :) = results_copy.avg_iterations_pYAPAll(1, :);
results.std_iterations_pYAPAll(4, :) = results_copy.std_iterations_pYAPAll(1, :);


figure('Position', [10 10 900 350]);
name = string();
name(1) = ['{\it{N}} = 1']; name(2) = ['{\it{N}} = 9'];

data_1 = [results(1).avg_iterations_pYAPAll{1,1} results(1).avg_iterations_pYAPAll{2,1} results(1).avg_iterations_pYAPAll{3,1} results(1).avg_iterations_pYAPAll{4,1}; ...
          results(1).avg_iterations_pYAPAll{1,2} results(1).avg_iterations_pYAPAll{2,2} results(1).avg_iterations_pYAPAll{3,2} results(1).avg_iterations_pYAPAll{4,2}]/particleNumber;

data_1_error = [results(1).std_iterations_pYAPAll{1,1} results(1).std_iterations_pYAPAll{2,1} results(1).std_iterations_pYAPAll{3,1} results(1).std_iterations_pYAPAll{4,1}; ...
                results(1).std_iterations_pYAPAll{1,2} results(1).std_iterations_pYAPAll{2,2} results(1).std_iterations_pYAPAll{3,2} results(1).std_iterations_pYAPAll{4,2}]/particleNumber;

bar_color = gray(6);
b = bar(data_1,'grouped'); hold on;
[ngroups,nbars] = size(data_1);

b(1).FaceColor = bar_color(2,:);
b(2).FaceColor = bar_color(3,:);
b(3).FaceColor = bar_color(4,:);
b(4).FaceColor = bar_color(5,:);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',data_1,data_1_error,'k','linestyle','none');
hold off;

ylabel('pYAP ratio');
xticklabels(name);

set(gca,'Fontname', 'Times New Roman','fontsize',18);
set(gca,'LineWidth',1);
legend('{\it{R}}_{b} = 50 s^{-1}, {\it{R}}_{u,YAP} = 0.1 s^{-1}', ...
    '{\it{R}}_{b} = 10 s^{-1}, {\it{R}}_{u,YAP} = 0.1 s^{-1}', ...
    '{\it{R}}_{b} = 50 s^{-1}, {\it{R}}_{u,YAP} = 0.5 s^{-1}', ...
    '{\it{R}}_{b} = 10 s^{-1}, {\it{R}}_{u,YAP} = 0.5 s^{-1}', ...
    'location','northeastoutside','Box','off','fontsize',14);
% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

folder = 'figures\';
print([folder 'fig3a'],'-dpng','-r300');

%% fig 3b 

clc
clear
close all

load('20230817T000936_results\sspYAP.mat');

for i = 1:5
    r(i) = rates{i,1}(1);
end

close all
figure('Position', [10 10 600 400]); 

marker_type = {'o','square','^','>'};
line_style = {'-','--',':','-.'};

errorbar(r, [results.avg_iterations_pYAPAll{:,1}]/particleNumber, [results.std_iterations_pYAPAll{:,1}]/particleNumber,marker_type{1},'MarkerFaceColor','k','MarkerSize',8,'Color','k');
hold on;
errorbar(r, [results.avg_iterations_pYAPAll{:,2}]/particleNumber, [results.std_iterations_pYAPAll{:,2}]/particleNumber,marker_type{2},'MarkerFaceColor','k','MarkerSize',8,'Color','k');



F_1 = @(x_1,xdata)x_1(1)*xdata./(x_1(2)+xdata);
x0 = [10 10];
[x_1,resnorm,~,exitflag,output] = lsqcurvefit(F_1,x0,r,[results.avg_iterations_pYAPAll{:,1}]);
hold on
% plot(1:0.1:100,F_1(x_1,1:0.1:100)/particleNumber,'LineWidth',1,'Color',[0 0.4470 0.7410]);
plot(1:0.1:100,F_1(x_1,1:0.1:100)/particleNumber,'LineWidth',1,'LineStyle',line_style{1},'Color','k');
fit_params(1,:) = x_1;

F_2 = @(x_2,xdata)x_2(1)*xdata./(x_2(2)+xdata);
x0 = [10 10];
[x_2,resnorm,~,exitflag,output] = lsqcurvefit(F_2,x0,r,[results.avg_iterations_pYAPAll{:,2}]);
hold on
% plot(1:0.1:100,F_2(x_2,1:0.1:100)/particleNumber,'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
plot(1:0.1:100,F_2(x_2,1:0.1:100)/particleNumber,'LineWidth',1,'LineStyle',line_style{2},'Color','k');

hold on
fit_params(2,:) = x_2

clear
marker_type = {'o','square','^','>'};
line_style = {'-','--',':','-.'};

load('20230818T203559_results\sspYAP.mat');
results_copy = results;

load('20230817T223348_results\sspYAP.mat');

results.avg_iterations_pYAPAll(2:5, 2) = results_copy.avg_iterations_pYAPAll(2:5, 2);
results.std_iterations_pYAPAll(2:5, 2) = results_copy.std_iterations_pYAPAll(2:5, 2);

for i = 1:5
    r(i) = rates{i,1}(1);
end
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
errorbar(r, [results.avg_iterations_pYAPAll{:,1}]/particleNumber, [results.std_iterations_pYAPAll{:,1}]/particleNumber,marker_type{3},'MarkerFaceColor','k','MarkerSize',8,'Color','k');
hold on;
errorbar(r, [results.avg_iterations_pYAPAll{:,2}]/particleNumber, [results.std_iterations_pYAPAll{:,2}]/particleNumber,marker_type{4},'MarkerFaceColor','k','MarkerSize',8,'Color','k');


F_1 = @(x_1,xdata)x_1(1)*xdata./(x_1(2)+xdata);
x0 = [10 10];
[x_1,resnorm,~,exitflag,output] = lsqcurvefit(F_1,x0,r,[results.avg_iterations_pYAPAll{:,1}]);
hold on
% plot(1:0.1:100,F_1(x_1,1:0.1:100)/particleNumber,'LineWidth',1,'Color',[0.9290    0.6940    0.1250]);
plot(1:0.1:100,F_1(x_1,1:0.1:100)/particleNumber,'LineWidth',1,'LineStyle',line_style{3},'Color','k');

fit_params(1,:) = x_1;

F_2 = @(x_2,xdata)x_2(1)*xdata./(x_2(2)+xdata);
x0 = [10 10];
[x_2,resnorm,~,exitflag,output] = lsqcurvefit(F_2,x0,r,[results.avg_iterations_pYAPAll{:,2}]);
hold on
plot(1:0.1:100,F_2(x_2,1:0.1:100)/particleNumber,'LineWidth',1,'LineStyle',line_style{4},'Color','k');
fit_params(2,:) = x_2


axis([0 110 0 0.82])

set(gca,'Fontname', 'Times New Roman','fontsize',18);
set(gca,'LineWidth',1);

xlabel('{\it{R}}_{b} (s^{-1})'); ylabel('pYAP ratio');
legend('{\it{N}} = 1, {\it{D}} = 0.8 {\mu}m^{2}/s','{\it{N}} = 9, {\it{D}} = 0.8 {\mu}m^{2}/s','','', ...
       '{\it{N}} = 1, {\it{D}} = 19 {\mu}m^{2}/s','{\it{N}} = 9, {\it{D}} = 19 {\mu}m^{2}/s', ...
       'location','southeast','Box','off','FontSize', 14);


% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

folder = 'figures\';
print([folder 'fig3b'],'-dpng','-r300');

    
%% fig 4

clc;
clear;
close all;

%##### low deph data

load('20230819T213159_results\sspYAP.mat'); %lifetime 60 
results_lt60_ld = results;

load('20230825T153857_results\sspYAP.mat'); % same range as previous no turnover
results_nolt_ld = results;

load('20230824T233621_results\sspYAP.mat'); %lifetime 135 
results_lt135_ld = results;

%##### high deph data

load('20230819T213159_results\sspYAP.mat'); %lifetime 60 
results_lt60_hd = results;

load('20230825T153857_results\sspYAP.mat'); % same range as previous no turnover
results_nolt_hd = results;

load('20230824T233621_results\sspYAP.mat'); %lifetime 135 
results_lt135_hd = results;

%
row_output_ld = 1:5;
row_output_hd = 6:10;

size_rate = size(rates);
for i = 1:length(row_output_ld)
    r_upyap(i) = rates{i,1}(4);
end

line_style = {'-','--',':'};
marker_type = {'o','square','^'};

color_def = [   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%pyap All low & high deph

figure('Position', [10 10 1700 500]); 
subplot(1,3,1)

errorbar(r_upyap, [results_nolt_ld.avg_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            [results_nolt_ld.std_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, [results_lt135_ld.avg_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            [results_lt135_ld.std_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, [results_lt60_ld.avg_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            [results_lt60_ld.std_iterations_pYAPAll{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;
 
%legend('No turnover', 'lifetime = 135 s', 'lifetime = 60 s','location','best','Box','off','FontSize', 12);

errorbar(r_upyap, [results_nolt_hd.avg_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
                [results_nolt_hd.std_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
                'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, [results_lt135_hd.avg_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
            [results_lt135_hd.std_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, [results_lt60_hd.avg_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
            [results_lt60_hd.std_iterations_pYAPAll{row_output_hd,1}]/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

legend('No turnover', 'lifetime = 135 s', 'lifetime = 60 s',...
       'No turnover', 'lifetime = 135 s', 'lifetime = 60 s','location','east','Box','off','FontSize', 12);

%  xlim([0.0008 0.25]);    %xlim([-0.01 0.22]);
%  set(gca, 'XScale', 'log');  
%  xticks([0.001, 0.01, 0.1, 0.2])
xlim([-0.01 0.22])
 set(gca,'Fontname', 'Times New Roman','fontsize',18);
 set(gca,'LineWidth',1);
 xlabel('{\it{R}}_{u,pYAP}'); ylabel('pYAP ratio');
% Get the current axis limits
ax = axis;

% Calculate the middle of the plot
mid_x = mean(ax(1:2));
mid_y = mean(ax(3:4));

% Add a text box in the middle of the plot
annotation('textbox', [mid_x + 0.06, mid_y + 0.08, 0.21, 0.21], 'String', ...
    '{\it{R}}_{deph}= 0.035 s^{-1}', 'FitBoxToText', 'on', 'EdgeColor','none','FontSize', 12,'Color', 'k', 'FontName', 'Times New Roman');

annotation('textbox', [mid_x + 0.06, mid_y - 0.03, 0.21, 0.21], 'String', ...
    '{\it{R}}_{deph}= 0.56 s^{-1}', 'FitBoxToText', 'on', 'EdgeColor','none','FontSize', 12,'Color', color_def(1,:), 'FontName', 'Times New Roman');

% pYAP in low and high deph

subplot(1,3,2);

errorbar(r_upyap, ([results_nolt_ld.avg_iterations_pYAPAll{row_output_ld,1}]-[results_nolt_ld.avg_iterations_pYAPOut{row_output_ld,1}])/particleNumber, ...
            sqrt([results_nolt_ld.std_iterations_pYAPAll{row_output_ld,1}].^2+[results_nolt_ld.std_iterations_pYAPOut{row_output_ld,1}].^2)/particleNumber, ...
            'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, ([results_lt135_ld.avg_iterations_pYAPAll{row_output_ld,1}]-[results_lt135_ld.avg_iterations_pYAPOut{row_output_ld,1}])/particleNumber, ...
            sqrt([results_lt135_ld.std_iterations_pYAPAll{row_output_ld,1}].^2+[results_lt135_ld.std_iterations_pYAPOut{row_output_ld,1}].^2)/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, ([results_lt60_ld.avg_iterations_pYAPAll{row_output_ld,1}]-[results_lt60_ld.avg_iterations_pYAPOut{row_output_ld,1}])/particleNumber, ...
            sqrt([results_lt60_ld.std_iterations_pYAPAll{row_output_ld,1}].^2+[results_lt60_ld.std_iterations_pYAPOut{row_output_ld,1}].^2)/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;


errorbar(r_upyap, ([results_nolt_hd.avg_iterations_pYAPAll{row_output_hd,1}]-[results_nolt_hd.avg_iterations_pYAPOut{row_output_hd,1}])/particleNumber, ...
            sqrt([results_nolt_hd.std_iterations_pYAPAll{row_output_hd,1}].^2+[results_nolt_hd.std_iterations_pYAPOut{row_output_hd,1}].^2)/particleNumber, ...
            'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, ([results_lt135_hd.avg_iterations_pYAPAll{row_output_hd,1}]-[results_lt135_hd.avg_iterations_pYAPOut{row_output_hd,1}])/particleNumber, ...
            sqrt([results_lt135_hd.std_iterations_pYAPAll{row_output_hd,1}].^2+[results_lt135_hd.std_iterations_pYAPOut{row_output_hd,1}].^2)/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, ([results_lt60_hd.avg_iterations_pYAPAll{row_output_hd,1}]-[results_lt60_hd.avg_iterations_pYAPOut{row_output_hd,1}])/particleNumber, ...
            sqrt([results_lt60_hd.std_iterations_pYAPAll{row_output_hd,1}].^2+[results_lt60_hd.std_iterations_pYAPOut{row_output_hd,1}].^2)/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;


%  xlim([0.0008 0.25]);    %xlim([-0.01 0.22]);
%  set(gca, 'XScale', 'log');  
%  xticks([0.001, 0.01, 0.1, 0.2])
xlim([-0.01 0.22])
 set(gca,'Fontname', 'Times New Roman','fontsize',18);
 set(gca,'LineWidth',1);
 xlabel('{\it{R}}_{u,pYAP}'); ylabel(sprintf('pYAP ratio\n inside the adhesions'));
% pYAP out low and high deph

subplot(1,3,3);

errorbar(r_upyap, [results_nolt_ld.avg_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            [results_nolt_ld.std_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, [results_lt135_ld.avg_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            [results_lt135_ld.std_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;

errorbar(r_upyap, [results_lt60_ld.avg_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            [results_lt60_ld.std_iterations_pYAPOut{row_output_ld,1}]/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor','k','MarkerSize',8, LineWidth=1,Color='k'); hold on;



errorbar(r_upyap, [results_nolt_hd.avg_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            [results_nolt_hd.std_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            'Marker',marker_type{1},'LineStyle',line_style{1},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, [results_lt135_hd.avg_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            [results_lt135_hd.std_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            'Marker',marker_type{3},'LineStyle',line_style{3},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;

errorbar(r_upyap, [results_lt60_hd.avg_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            [results_lt60_hd.std_iterations_pYAPOut{row_output_hd,1}]/particleNumber, ...
            'Marker',marker_type{2},'LineStyle',line_style{2},'MarkerFaceColor',color_def(1,:),'MarkerSize',8, LineWidth=1,Color=color_def(1,:)); hold on;




%  xlim([0.0008 0.25]);    %xlim([-0.01 0.22]);
%  set(gca, 'XScale', 'log');  
%  xticks([0.001, 0.01, 0.1, 0.2])
xlim([-0.01 0.22])
 set(gca,'Fontname', 'Times New Roman','fontsize',18);
 set(gca,'LineWidth',1);
 xlabel('{\it{R}}_{u,pYAP}'); ylabel(sprintf('pYAP ratio\n outside the adhesions'));

folder = 'figures\';
print([folder 'fig4'],'-dpng','-r300');

%% fig S3

clc;
clear;
close all;

load('20230824T233530_results\sspYAP.mat');

ns = [3 5 15];
figure('Position', [10 10 650 900]);

for i = 1:length(ns) % sets
    %figure('Position', [10 10 650 200]); 
    
    subplot(4,1,i)
    plot(results(1).time{1,4}(1:ns(i),:)', results(1).pYAPAll{1,4}(1:ns(i),:)'/particleNumber,'LineWidth',2);
    axis([0 425 0 0.6])
    xlabel('Time (s)'); ylabel('pYAP ratio');
    title([num2str(ns(i)) ' simulations']);
    set(gca,'Fontname', 'Times New Roman','fontsize',14);
    set(gca,'LineWidth',1);
    % Get the current axes handle
    ax = gca;
    % Remove the top and right axes by setting the 'box' property to 'off'
    ax.Box = 'off';


end


% plotting avg pYAP vs number of simulations
%figure('Position', [10 10 650 300]); 
subplot(4,1,4)
ns = [2 3 5 10 15];

for i = 1:length(ns)
errorbar(ns(i), mean(results.avgpYAPAll{1,4}(1:ns(i)))/particleNumber, std(results.avgpYAPAll{1,4}(1:ns(i)))/particleNumber, 'vertical', 'o', 'MarkerSize', 8, 'MarkerFaceColor',[ 0    0.4470    0.7410], 'Color',[ 0    0.4470    0.7410]);
axis([1 16 0.4 0.6]);
hold on;
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';
end

xlabel('No. simulations'); ylabel('pYAP ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);

folder = 'figures\'
print([folder 'figS3'],'-dpng','-r300');

%% fig. S4

clc;
clear;
close all;

load('20230829T122341_results\sspYAP.mat'); 

figure('Position', [10 10 800 400]); 

plot(results(1).time{1,1}(:,1:15*10^6)', results(1).pYAPAll{1,1}(:,1:15*10^6)'/particleNumber,'LineWidth',2,'Color',	[0  0.4470 0.7410]);
xlabel('Time (s)'); ylabel('pYAP ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
hold on;

plot(results(1).time{1,2}(:,1:15*10^6)', results(1).pYAPAll{1,2}(:,1:15*10^6)'/particleNumber,'LineWidth',2,'Color',	[0.8500  0.3250  0.0980]);
%xlabel('Time (s)'); ylabel('pYAP ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);

axis([0 300  0 0.4])
legend('{\it{N}} = 1','','','{\it{N}} = 9','','','location','northeastoutside','Box','off');

% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

folder = 'figures\';
print([folder 'figS4'],'-dpng','-r300');
%% fig S5

clc;
clear;
close all;

load('20230824T233530_results\sspYAP.mat');

% calculate pairwise distances
 
for i = 1:4 % N= 1,2,4,9
    for j = 1:15   % no. simulations    
        pd_integ{1,i}(j) = mean(pdist(results.integpos{1,i}{j}(:,1:2)));
        
    end
    pd_integ_avg(i) = mean(pd_integ{1,i})*0.2;
    pd_integ_std(i) = std(pd_integ{1,i})*0.2;
end

avg_pyapAll = cellfun(@mean, results.avgpYAPAll)/particleNumber;
std_pyapAll = cellfun(@std, results.avgpYAPAll)/particleNumber;

close all
figure('Position', [10 10 500 500]); 

color_def = [   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

marker_type = {'o','square','^','>'};

for i = 1:4 
    errorbar(pd_integ_avg(i), avg_pyapAll(i), std_pyapAll(i), 'vertical', marker_type{i}, 'MarkerSize', 8, 'MarkerFaceColor','k', 'Color', 'k', 'LineWidth', 1.5);
    hold on;
    errorbar(pd_integ_avg(i), avg_pyapAll(i), pd_integ_std(i), 'horizontal', marker_type{i}, 'MarkerSize', 8, 'MarkerFaceColor','k', 'Color', 'k', 'LineWidth', 1.5);
end

axis([0.5 4.5 0.3 0.55])

xlabel('pairwise distance (\mum)'); ylabel('pYAP ratio');

set(gca,'Fontname', 'Times New Roman','fontsize',14);

set(gca,'Fontname', 'Times New Roman','fontsize',18);
set(gca,'LineWidth',1);

hold on;

% Initialize variables to store fit results
best_degree = 0; % Degree with the best fit
best_rmse = Inf; % Initialize RMSE with a high value

for degree = 2:2
    % Fit a polynomial of the current degree
    coeffs = polyfit(pd_integ_avg, avg_pyapAll, degree);
    
    % Evaluate the polynomial on the x-values
    fitted_pyapAll = polyval(coeffs, pd_integ_avg);
    
    % Calculate the Root Mean Square Error (RMSE) as a measure of accuracy
    rmse = sqrt(mean((fitted_pyapAll - avg_pyapAll).^2));
    
    % Check if this degree provides a better fit
    if rmse < best_rmse
        best_rmse = rmse;
        best_degree = degree;
    end
    
    % Optionally, you can plot the fitted curve for visualization
%     figure;
%     plot(x, fitted_y);
    
end

hold on ;
plot(0.7:0.1:4.3,polyval(coeffs, 0.7:0.1:4.3),'LineWidth',1.5,'Color',[150 150 150]/250,'LineStyle','--');

 legend('{\it{N}} = 1','','{\it{N}} = 2','','{\it{N}} = 4','','{\it{N}} = 9','','location','northwest','Box','off');
% Get the current axes handle
ax = gca;

% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';


folder = 'figures\';
print([folder 'figS5'],'-dpng','-r300');



%% fig S6,S8

clc;
clear;
close all;

load('20230819T213159_results\sspYAP.mat'); %lifetime 60 factor 0.04
results_lt60 = results;

load('20230825T153857_results\sspYAP.mat'); % same range as previous no turnover
results_nolt = results;

load('20230824T233621_results\sspYAP.mat'); %lifetime 135 factor 0.04
results_lt135 = results;

color_def = [   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

close all
figure('Position', [10 10 1200 700]); 

 % nolt

subplot(3,3,1)
plot(results_nolt.time{2,1}', results_nolt.pYAPAll{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
%ylabel('pYAP/TAZ ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,2)
plot(results_nolt.time{2,1}', (results_nolt.pYAPAll{2,1}'-results_nolt.pYAPOut{2,1}')/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
%ylabel('pYAP/TAZ ratio inside');
title('No adhesion turnover');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,3)
plot(results_nolt.time{2,1}', results_nolt.pYAPOut{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
% ylabel('pYAP/TAZ ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

% lifetime 135

subplot(3,3,4)
plot(results_lt135.time{2,1}', results_lt135.pYAPAll{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
ylabel('pYAP ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,5)
plot(results_lt135.time{2,1}', (results_lt135.pYAPAll{2,1}'-results_lt135.pYAPOut{2,1}')/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
ylabel('pYAP ratio inside');
title('lifetime = 135 s');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,6)
plot(results_lt135.time{2,1}', results_lt135.pYAPOut{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
ylabel('pYAP ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

%lifetime 60

subplot(3,3,7)
plot(results_lt60.time{2,1}', results_lt60.pYAPAll{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
%ylabel('pYAP/TAZ ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,8)
plot(results_lt60.time{2,1}', (results_lt60.pYAPAll{2,1}'-results_lt60.pYAPOut{2,1}')/particleNumber,'LineWidth',1);
xlabel('Time (s)'); 
%ylabel('pYAP/TAZ ratio inside');
title('lifetime = 60 s');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,9)
plot(results_lt60.time{2,1}', results_lt60.pYAPOut{2,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
% ylabel('pYAP/TAZ ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

folder = 'figures\'
print([folder 'figS6' ],'-dpng','-r300');

%#################################### Ru = 0.001
figure('Position', [10 10 1200 700]); 
row = 5;
 % nolt

subplot(3,3,1)
plot(results_nolt.time{row,1}', results_nolt.pYAPAll{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
%ylabel('pYAP/TAZ ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,2)
plot(results_nolt.time{row,1}', (results_nolt.pYAPAll{row,1}'-results_nolt.pYAPOut{row,1}')/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
%ylabel('pYAP/TAZ ratio inside');
title('No adhesion turnover');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,3)
plot(results_nolt.time{row,1}', results_nolt.pYAPOut{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
% ylabel('pYAP/TAZ ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

% lifetime 135

subplot(3,3,4)
plot(results_lt135.time{row,1}', results_lt135.pYAPAll{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
ylabel('pYAP ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,5)
plot(results_lt135.time{row,1}', (results_lt135.pYAPAll{row,1}'-results_lt135.pYAPOut{row,1}')/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
ylabel('pYAP ratio inside');
title('lifeimte = 135 s');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,6)
plot(results_lt135.time{row,1}', results_lt135.pYAPOut{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
ylabel('pYAP ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

%lifetime 60

subplot(3,3,7)
plot(results_lt60.time{row,1}', results_lt60.pYAPAll{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)');
%ylabel('pYAP/TAZ ratio');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,8)
plot(results_lt60.time{row,1}', (results_lt60.pYAPAll{row,1}'-results_lt60.pYAPOut{row,1}')/particleNumber,'LineWidth',1);
xlabel('Time (s)'); 
%ylabel('pYAP/TAZ ratio inside');
title('lifetime = 60 s');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

subplot(3,3,9)
plot(results_lt60.time{row,1}', results_lt60.pYAPOut{row,1}'/particleNumber,'LineWidth',1);
%xlabel('Time (s)'); 
% ylabel('pYAP ratio outside');
set(gca,'Fontname', 'Times New Roman','fontsize',14);
set(gca,'LineWidth',1);
% Get the current axes handle
ax = gca;
% Remove the top and right axes by setting the 'box' property to 'off'
ax.Box = 'off';

print([folder 'figS8' ],'-dpng','-r300');



