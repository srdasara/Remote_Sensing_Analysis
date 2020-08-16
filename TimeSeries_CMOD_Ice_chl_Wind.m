

%% Here, I plot out my saved variables from 'Saving_out_CMOD_Ice_Wind.m' and plot them as monthly timeseries to view trends between CMOD, Ice, Chl-a, and Wind Speed



cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

%%

load CMOD_Monthly_avg_Vars.mat
load Depol_Ratio_Monthly_avg_Vars.mat
load amsrmf_Monthly_avg_Vars.mat

%%

cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/TimeSeries_Figures

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 
clear t1 t2
    
%%

% Line Plot Code...

black = rgb('black'); 
grass_green = rgb('leaf green');
ice_blue = rgb('lightish blue'); 
wind_blue = rgb('royal blue');

%% 
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY/
load('Master_chlor_a_monthly_full_res.mat') 
for i = 1:151 
Total_chl_monthly_2(i) = nanmean(Master_chl_a(:,:,i), [1 2]); 
end

Total_chl_a_monthly = Total_chl_monthly_2; 

cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/TimeSeries_Figures


%%
%%%%%% ----- Timeseries plots below ------ %%%%%%%%%
% doc addaxis If needed
% 
% 
% black = rgb('black'); 
% grass_green = rgb('true green');
% ice_blue = rgb('lightish blue'); 
% ice_blue = rgb('dodger blue');
% wind_blue = rgb('royal blue');
fig = figure; clf;
    
set(fig, 'defaultAxesColorOrder', [black; grass_green; ice_blue]);

   x = times;
   
   aa_splot(x, CMOD_Monthly_avg, '-',...
       'linewidth', 1.5, ...
       'Color', black);
% %    xticks(times(1) : calmonths(12) :times(end))
% calendarmonths = times(1) : calmonths(6) : times(end); 
% X = month(calendarmonths); 
set(gca, 'XTick', (times(1) : calmonths(6) : times(end)));

 xtickformat('MM-yyyy')
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = times(1) : calmonths(3) : times(end); % 
grid on
h.XMinorGrid = 'on';
% 


% set(gca, 'XTick', X);

xtickangle(35)
   
% CHLOROPHYLL
   addaxis(x,  Total_chl_a_monthly,'-',...
       'linewidth', 1.5,...
       'MarkerSize', 4,...
       'MarkerFaceColor', grass_green,...
       'MarkerEdgeColor', grass_green,...
       'Color', grass_green);

% ICE
   addaxis(x,  Depol_Ratio_Monthly_avg,'-',...
    'linewidth', 1.5,... 
    'MarkerSize', 4,...
    'MarkerFaceColor', ice_blue,...
    'MarkerEdgeColor', ice_blue,...
    'Color', ice_blue);
%    
% % WIND
%     
%    addaxis(x, amsrmf_Monthly_avg, 'x-',...
%     'linewidth', 1.5,... 
%     'MarkerSize', 4,...
%     'MarkerFaceColor', wind_blue,...
%     'MarkerEdgeColor', wind_blue,...
%     'Color', wind_blue);
%    
   
   addaxislabel(1,'CMOD Night');
   addaxislabel(2,'Chl-{\ita} (mg m^{-3})');
   addaxislabel(3,'Ice'); 
%    addaxislabel(4, 'Wind'); 
%    

   AX=findall(0,'type','axes'); 
   set(AX, 'FontSize', 25)

   grid on
% Make sure this is in the right order
 legend('CMOD', 'Chl-{\ita}', 'Ice'); 
%  title('Timeseries: CMOD, Chlorophyll-{\ita} concentration, & Ice Depolarization Ratio (\delta)')
 
 %%
 
 set(gcf,'PaperPositionMode','auto')
%  set(gcf,'PaperPosition','fillpage') 
%  orient(fig,'portrait')

print(fig,'Updated_TimeSeries_CMOD_Ice_CHL_for_publication_v5.png','-dpng','-r96');       %  *// 300 dpi
%  print(gcf, 'June_test.pdf', '-dpdf', '-r96', '-fillpage');

 %%
 
wind_blue = rgb('royal blue');
fig = figure; clf;
    
set(fig, 'defaultAxesColorOrder', [black; wind_blue]);

   x = times;
   
   aa_splot(x, CMOD_Monthly_avg, '-',...
       'linewidth', 1.5, ...
       'Color', black);
set(gca, 'XTick', (times(1) : calmonths(6) : times(end)));

 xtickformat('MM-yyyy')
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = times(1) : calmonths(3) : times(end); % 
grid on
h.XMinorGrid = 'on';

xtickangle(35)

% WIND
    
   addaxis(x, amsrmf_Monthly_avg, '-',...
    'linewidth', 1.5,... 
    'MarkerSize', 4,...
    'MarkerFaceColor', wind_blue,...
    'MarkerEdgeColor', wind_blue,...
    'Color', wind_blue);
   
   
   addaxislabel(1,'CMOD Night');
   addaxislabel(2,'Wind Speed (m s^{-1})');
%    
   AX=findall(0,'type','axes'); 
   set(AX, 'FontSize', 25)

   grid on
   
% Make sure this is in the right order
 legend('CMOD', 'Wind Speed'); 
%  title('Timeseries: CMOD & Wind Speed')
 

%%
 
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'Updated_TimeSeries_CMOD_Wind_for_publication_v5.png','-dpng','-r96');       %  *// 300 dpi


%%







%%
%%%%%% ----- Timeseries plots below ------ %%%%%%%%%
% doc addaxis If needed
% 
% 
% black = rgb('black'); 
% grass_green = rgb('true green');
% ice_blue = rgb('lightish blue'); 
% ice_blue = rgb('dodger blue');
% wind_blue = rgb('royal blue');
fig = figure; clf;

subplot(2,1,1)
set(fig, 'defaultAxesColorOrder', [black; grass_green; ice_blue]);

   x = times;
   
   aa_splot(x, CMOD_Monthly_avg, '-',...
       'linewidth', 1.5, ...
       'Color', black);
% %    xticks(times(1) : calmonths(12) :times(end))
% calendarmonths = times(1) : calmonths(6) : times(end); 
% X = month(calendarmonths); 
set(gca, 'XTick', (times(1) : calmonths(6) : times(end)));

 xtickformat('MM-yyyy')
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = times(1) : calmonths(3) : times(end); % 
grid on
h.XMinorGrid = 'on';
% 


% set(gca, 'XTick', X);

xtickangle(35)
   
% CHLOROPHYLL
   addaxis(x,  Total_chl_a_monthly,'-',...
       'linewidth', 1.5,...
       'MarkerSize', 4,...
       'MarkerFaceColor', grass_green,...
       'MarkerEdgeColor', grass_green,...
       'Color', grass_green);

% ICE
   addaxis(x,  Depol_Ratio_Monthly_avg,'-',...
    'linewidth', 1.5,... 
    'MarkerSize', 4,...
    'MarkerFaceColor', ice_blue,...
    'MarkerEdgeColor', ice_blue,...
    'Color', ice_blue);
%    
% % WIND
%     
%    addaxis(x, amsrmf_Monthly_avg, 'x-',...
%     'linewidth', 1.5,... 
%     'MarkerSize', 4,...
%     'MarkerFaceColor', wind_blue,...
%     'MarkerEdgeColor', wind_blue,...
%     'Color', wind_blue);
%    
   
   addaxislabel(1,'CMOD Night');
   addaxislabel(2,'Chl-{\ita} (mg m^{-3})');
   addaxislabel(3,'Ice'); 
%    addaxislabel(4, 'Wind'); 
%    

   AX=findall(0,'type','axes'); 
   set(AX, 'FontSize', 25)

   grid on
% Make sure this is in the right order
 legend('CMOD', 'Chl-{\ita}', 'Ice'); 
%  title('Timeseries: CMOD, Chlorophyll-{\ita} concentration, & Ice Depolarization Ratio (\delta)')
 
subplot(2,1,2)
 
 
wind_blue = rgb('royal blue');
% fig = figure; clf;
    
set(fig, 'defaultAxesColorOrder', [black; wind_blue]);

   x = times;
   
   aa_splot(x, CMOD_Monthly_avg, '-',...
       'linewidth', 1.5, ...
       'Color', black);
set(gca, 'XTick', (times(1) : calmonths(6) : times(end)));

 xtickformat('MM-yyyy')
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = times(1) : calmonths(3) : times(end); % 
grid on
h.XMinorGrid = 'on';

xtickangle(35)

% WIND
    
   addaxis(x, amsrmf_Monthly_avg, '-',...
    'linewidth', 1.5,... 
    'MarkerSize', 4,...
    'MarkerFaceColor', wind_blue,...
    'MarkerEdgeColor', wind_blue,...
    'Color', wind_blue);
   
   
   addaxislabel(1,'CMOD Night');
   addaxislabel(2,'Wind Speed (m s^{-1})');
%    
   AX=findall(0,'type','axes'); 
   set(AX, 'FontSize', 25)

   grid on
   
% Make sure this is in the right order
 legend('CMOD', 'Wind Speed'); 
%  title('Timeseries: CMOD & Wind Speed')
 

%%







