

% Here I plot out climatologies both spatially and temporally. Spatial
% climatologies will show seasonal patterns, and temporal climatologies
% will show monthly patterns. All variables of interest are analyzed
% within, including CMOD, chl-a, Ice, and wind speed. 

%%


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication
load('Total_timetable_amsrmf.mat')
load('Total_timetable_Depol_Ratio.mat')
load('Total_timetable_CMOD.mat')


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY

load('Aqua.mat') 
load('Aqua_Longitude_Subset_BellingshausenSea.mat') 
load('Aqua_Latitude_Subset_BellingshausenSea.mat')
load('Master_chlor_a_monthly_full_res.mat') 

cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 
clear t1 t2
OPTIONS.MaxIter = 1000;


    
%%

Aqua_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
Aqua_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 

Aqua_Lat = flip(Aqua_Lat) ; 

%%
% I can index out the chlorophyll master array which is 1104 x 360 x 151
% into separate seasons for Winter, Spring, Summer, & Fall 

% After that, the climatology would just be the average of those seasons. 

% For Winter: 

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 

times = times';

[winter_x, winter_y] = find(times.Month >= 6 & times.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times(winter_x)
Chl_a_winter = Master_chl_a; 
Chl_a_winter = Chl_a_winter(:,:,winter_x);

% For Spring: 
[spring_x, spring_y] = find(times.Month >= 9 & times.Month <= 11);
times(spring_x)
Chl_a_spring = Master_chl_a; 
Chl_a_spring = Chl_a_spring(:,:, spring_x); 


% For Summer: 

[summer_x, summer_y] = find(times.Month >= 1 & times.Month <=2 | times.Month == 12);

times(summer_x)
Chl_a_summer = Master_chl_a; 
Chl_a_summer = Chl_a_summer(:,:, summer_x); 

% For Fall: 

[fall_x, fall_y] = find(times.Month >= 3 & times.Month <= 5); 
times(fall_x) 
Chl_a_fall = Master_chl_a; 
Chl_a_fall = Chl_a_fall(:,:, fall_x); 

% So now for the raw and smoothn climatologies, you can average all of
% these and plot... (for raw) or average all of them and then smoothn on
% the 2D. 

Chl_a_winter_mean = mean(Chl_a_winter,3 ,'omitnan');
Chl_a_spring_mean = mean(Chl_a_spring, 3 , 'omitnan'); 
Chl_a_summer_mean = mean(Chl_a_summer, 3, 'omitnan'); 
Chl_a_fall_mean = mean(Chl_a_fall, 3, 'omitnan'); 

%% Adjust resolution 

step_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
step_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 
step_Lat = flip(step_Lat); 

 Chl_a_winter_mean_one_degree_res = BlockMean(Chl_a_winter_mean,...
        (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat)); 

 Chl_a_spring_mean_one_degree_res = BlockMean(Chl_a_spring_mean,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));
    
 Chl_a_summer_mean_one_degree_res = BlockMean(Chl_a_summer_mean,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));

 Chl_a_fall_mean_one_degree_res = BlockMean(Chl_a_fall_mean,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));

% climatology chl raw full degree Res 

Chl_a_summer_mean_rot = rot90(Chl_a_summer_mean); 
Chl_a_spring_mean_rot = rot90(Chl_a_spring_mean); 
Chl_a_fall_mean_rot   = rot90(Chl_a_fall_mean); 
Chl_a_winter_mean_rot = rot90(Chl_a_winter_mean); 

%% plot climatology chl raw figures full resolution:

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(2,2,3)
plot_BSea_figure_subplot(Chl_a_summer_mean_rot, ...
    Longitude_subset,...
    flip(Latitude_subset),...
    'algae', ...
    [0 1.5],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(Chl_a_spring_mean_rot, ...
    Longitude_subset,...
    flip(Latitude_subset),...
    'algae', ...
    [0 1.5],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(Chl_a_fall_mean_rot, ...
    Longitude_subset,...
    flip(Latitude_subset),...
    'algae', ...
    [0 1.5],...
    'Fall'); 
    
subplot(2,2,1)
plot_BSea_figure_subplot(Chl_a_winter_mean_rot, ...
    Longitude_subset,...
    flip(Latitude_subset),...
    'algae', ...
    [0 1.5],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, sprintf('mg m^{-3}'), 'FontSize', 18);
% hy.FontSize = 18;



%% trying to subplot seasonal climatologies so that they are all in the same figure. 
   
% now smoothn full res files, and  adjust the resolution 

Chl_a_winter_mean_smoothn = smoothn(Chl_a_winter_mean, 'robust', OPTIONS); 
Chl_a_spring_mean_smoothn = smoothn(Chl_a_spring_mean, 'robust', OPTIONS); 
Chl_a_summer_mean_smoothn = smoothn(Chl_a_summer_mean, 'robust', OPTIONS); 
Chl_a_fall_mean_smoothn = smoothn(Chl_a_fall_mean, 'robust', OPTIONS); 

 Chl_a_winter_mean_one_degree_smoothn = BlockMean(Chl_a_winter_mean_smoothn,...
        (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat)); 

 Chl_a_spring_mean_one_degree_smoothn = BlockMean(Chl_a_spring_mean_smoothn,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));
    
 Chl_a_summer_mean_one_degree_smoothn = BlockMean(Chl_a_summer_mean_smoothn,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));

 Chl_a_fall_mean_one_degree_smoothn = BlockMean(Chl_a_fall_mean_smoothn,...
     (length(Longitude_subset) ./ length(step_Lon)) ,...
        length(Latitude_subset)  ./ length(step_Lat));

% chl vars: gaps from full res interpolated before constructing one degree
% res figures. 

% one degree Res vars: 
Chl_a_winter_mean_smoothn_rot = rot90(Chl_a_winter_mean_one_degree_smoothn); 
Chl_a_spring_mean_smoothn_rot = rot90(Chl_a_spring_mean_one_degree_smoothn); 
Chl_a_summer_mean_smoothn_rot = rot90(Chl_a_summer_mean_one_degree_smoothn); 
Chl_a_fall_mean_smoothn_rot   = rot90(Chl_a_fall_mean_one_degree_smoothn); 


%%
% If you want the image to have land pixels masked, load 'WAP_landmask_updated.mat', below. I felt this was
% unecessary for the purpose of visualizaiton. 

% cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/LandPoints
% load('WAP_landmask_updated.mat') 

%%
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Climatology_Figures
%%
% winter = WAP_landmask_updated .* Chl_a_winter_mean_smoothn_rot; 
% spring = WAP_landmask_updated .* Chl_a_spring_mean_smoothn_rot; 
% summer = WAP_landmask_updated.*Chl_a_summer_mean_smoothn_rot;
% fall = WAP_landmask_updated .* Chl_a_fall_mean_smoothn_rot;

% one degree res... 
winter = Chl_a_winter_mean_smoothn_rot; 
spring = Chl_a_spring_mean_smoothn_rot; 
summer = Chl_a_summer_mean_smoothn_rot;
fall =  Chl_a_fall_mean_smoothn_rot;

%%

% here i'm plotting out the full res chl-climatologies that have also been
% interpolated 

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure,clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    'algae', ...
    [0 1.5],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    'algae', ...
    [0 1.5],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    'algae', ...
    [0 1.5],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    'algae', ...
    [0 1.5],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, sprintf('mg m^{-3}'), 'FontSize', 18);
% hy.FontSize = 18;





%% CMOD, Winds, & Ice Now 

   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
   
    TR_winter_2006 = timerange('2006-06-01', '2006-09-01'); % Winter: June 1st to September 1st
    TR_winter_2007 = timerange('2007-06-01', '2007-09-01'); 
    TR_winter_2008 = timerange('2008-06-01', '2008-09-01');
    TR_winter_2009 = timerange('2009-06-01', '2009-09-01');
    TR_winter_2010 = timerange('2010-06-01', '2010-09-01'); 
    TR_winter_2011 = timerange('2011-06-01', '2011-09-01');
    TR_winter_2012 = timerange('2012-06-01', '2012-09-01');
    TR_winter_2013 = timerange('2013-06-01', '2013-09-01');
    TR_winter_2014 = timerange('2014-06-01', '2014-09-01');
    TR_winter_2015 = timerange('2015-06-01', '2015-09-01');
    TR_winter_2016 = timerange('2016-06-01', '2016-09-01');
    TR_winter_2017 = timerange('2017-06-01', '2017-09-01');
    TR_winter_2018 = timerange('2018-06-01', '2018-09-01'); 
    
    TR_spring_2006 = timerange('2006-09-01', '2006-12-01'); % Spring: Sept 1st to Dec 1st
    TR_spring_2007 = timerange('2007-09-01', '2007-12-01');
    TR_spring_2008 = timerange('2008-09-01', '2008-12-01'); 
    TR_spring_2009 = timerange('2009-09-01', '2009-12-01'); 
    TR_spring_2010 = timerange('2010-09-01', '2010-12-01'); 
    TR_spring_2011 = timerange('2011-09-01', '2011-12-01');
    TR_spring_2012 = timerange('2012-09-01', '2012-12-01'); 
    TR_spring_2013 = timerange('2013-09-01', '2013-12-01'); 
    TR_spring_2014 = timerange('2014-09-01', '2014-12-01');
    TR_spring_2015 = timerange('2015-09-01', '2015-12-01'); 
    TR_spring_2016 = timerange('2016-09-01', '2016-12-01'); 
    TR_spring_2017 = timerange('2017-09-01', '2017-12-01');
    TR_spring_2018 = timerange('2018-09-01', '2018-12-01'); 
    
    TR_summer_2006 = timerange('2006-12-01', '2007-03-01'); % Summer: Dec 1st to March 1st
    TR_summer_2007 = timerange('2007-12-01', '2008-03-01');
    TR_summer_2008 = timerange('2008-12-01', '2009-03-01'); 
    TR_summer_2009 = timerange('2009-12-01', '2010-03-01'); 
    TR_summer_2010 = timerange('2010-12-01', '2011-03-01'); 
    TR_summer_2011 = timerange('2011-12-01', '2012-03-01');
    TR_summer_2012 = timerange('2012-12-01', '2013-03-01'); 
    TR_summer_2013 = timerange('2013-12-01', '2014-03-01'); 
    TR_summer_2014 = timerange('2014-12-01', '2015-03-01');
    TR_summer_2015 = timerange('2015-12-01', '2016-03-01'); 
    TR_summer_2016 = timerange('2016-12-01', '2017-03-01'); 
    TR_summer_2017 = timerange('2017-12-01', '2018-03-01');
    TR_summer_2018 = timerange('2018-12-01', '2019-01-01'); 
    
    %     TR_fall_2006   = timerange('2006-03-01', '2006-06-01'); % there
    %     is no CALIPSO DATA in fall of 2006
    TR_fall_2007   = timerange('2007-03-01', '2007-06-01'); %% Fall: March 1st to June 1st
    TR_fall_2008   = timerange('2008-03-01', '2008-06-01');
    TR_fall_2009   = timerange('2009-03-01', '2009-06-01');
    TR_fall_2010   = timerange('2010-03-01', '2010-06-01');
    TR_fall_2011   = timerange('2011-03-01', '2011-06-01');
    TR_fall_2012   = timerange('2012-03-01', '2012-06-01');
    TR_fall_2013   = timerange('2013-03-01', '2013-06-01');
    TR_fall_2014   = timerange('2014-03-01', '2014-06-01');
    TR_fall_2015   = timerange('2015-03-01', '2015-06-01');
    TR_fall_2016   = timerange('2016-03-01', '2016-06-01');
    TR_fall_2017   = timerange('2017-03-01', '2017-06-01');
    TR_fall_2018   = timerange('2018-03-01', '2018-06-01');

%%
% Initializing of my variables for the loop below. 

Total_CMOD_winter = [];
Total_Lat_CMOD_winter = [];
Total_Lon_CMOD_winter = [];

Total_CMOD_spring = [];
Total_Lat_CMOD_spring = [];
Total_Lon_CMOD_spring = [];

Total_CMOD_summer = [];
Total_Lat_CMOD_summer = [];
Total_Lon_CMOD_summer = [];

Total_CMOD_fall = [];
Total_Lat_CMOD_fall = [];
Total_Lon_CMOD_fall = [];


Total_Ice_winter = [];
Total_Lat_Ice_winter = [];
Total_Lon_Ice_winter = [];

Total_Ice_spring = [];
Total_Lat_Ice_spring = [];
Total_Lon_Ice_spring = [];

Total_Ice_summer = [];
Total_Lat_Ice_summer = [];
Total_Lon_Ice_summer = [];

Total_Ice_fall = [];
Total_Lat_Ice_fall = [];
Total_Lon_Ice_fall = [];

Total_Wind_winter = [];
Total_Lat_Wind_winter = [];
Total_Lon_Wind_winter = [];

Total_Wind_spring = [];
Total_Lat_Wind_spring = [];
Total_Lon_Wind_spring = [];

Total_Wind_summer = [];
Total_Lat_Wind_summer = [];
Total_Lon_Wind_summer = [];

Total_Wind_fall = [];
Total_Lat_Wind_fall = [];
Total_Lon_Wind_fall = [];

%%
for i = 2006:2018
    
    disp(i)
    
    SEASON = {'winter', 'spring', 'summer', 'fall'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here
            
            eval(sprintf('Lat_CMOD = Total_timetable_CMOD(TR_%s_%d,:).Total_Latitude_Night_Cloud_Free;', SEASON{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_CMOD(TR_%s_%d, :).Total_Longitude_Night_Cloud_Free;',SEASON{j}, i))
            eval(sprintf('OD       = Total_timetable_CMOD(TR_%s_%d, :).CMOD;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
            
            % I dont need to filter this out anymore since the
            % 'Saving_out_CMOD_Ice_Wind.m' file took care of this. 
%             
%             bad_Ice_values = Ice <= -0.2 | Ice > 1.2; 
%             Ice(bad_Ice_values) = NaN;
%             
%             nan_ice        = isnan(Ice(:,1)); 
%             Ice      = Ice(~nan_ice) ;
%             Lat_Ice  = Lat_Ice(~nan_ice); 
%             Lon_Ice  = Lon_Ice(~nan_ice); 
%             
            
            
            
            eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
            eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
            eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
         
            
            if j == 1
                
                Total_CMOD_winter = vertcat(Total_CMOD_winter, OD);                 
                Total_Lat_CMOD_winter = vertcat(Total_Lat_CMOD_winter, Lat_CMOD);
                Total_Lon_CMOD_winter = vertcat(Total_Lon_CMOD_winter, Lon_CMOD);
                
                Total_Ice_winter = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice); 
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice); 
                                
                Total_Wind_winter = vertcat(Total_Wind_winter, Wind); 
                Total_Lat_Wind_winter = vertcat(Total_Lat_Wind_winter, Lat_Wind); 
                Total_Lon_Wind_winter = vertcat(Total_Lon_Wind_winter, Lon_Wind); 
                
            elseif j == 2
                
                Total_CMOD_spring = vertcat(Total_CMOD_spring, OD);
                Total_Lat_CMOD_spring = vertcat(Total_Lat_CMOD_spring, Lat_CMOD);
                Total_Lon_CMOD_spring = vertcat(Total_Lon_CMOD_spring, Lon_CMOD);
                
                Total_Ice_spring = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice); 
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice); 
                                
                Total_Wind_spring = vertcat(Total_Wind_spring, Wind); 
                Total_Lat_Wind_spring = vertcat(Total_Lat_Wind_spring, Lat_Wind); 
                Total_Lon_Wind_spring = vertcat(Total_Lon_Wind_spring, Lon_Wind); 
                
            elseif j == 3
                
                Total_CMOD_summer = vertcat(Total_CMOD_summer, OD);               
                Total_Lat_CMOD_summer = vertcat(Total_Lat_CMOD_summer, Lat_CMOD);
                Total_Lon_CMOD_summer = vertcat(Total_Lon_CMOD_summer, Lon_CMOD);
                
                Total_Ice_summer = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice); 
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice); 
                                
                Total_Wind_summer = vertcat(Total_Wind_summer, Wind); 
                Total_Lat_Wind_summer = vertcat(Total_Lat_Wind_summer, Lat_Wind); 
                Total_Lon_Wind_summer = vertcat(Total_Lon_Wind_summer, Lon_Wind); 
                
            elseif j == 4
                
                Total_CMOD_fall = vertcat(Total_CMOD_fall, OD);                
                Total_Lat_CMOD_fall = vertcat(Total_Lat_CMOD_fall, Lat_CMOD);
                Total_Lon_CMOD_fall = vertcat(Total_Lon_CMOD_fall, Lon_CMOD);
                
                Total_Ice_fall = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice); 
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice); 
                                
                Total_Wind_fall = vertcat(Total_Wind_fall, Wind); 
                Total_Lat_Wind_fall = vertcat(Total_Lat_Wind_fall, Lat_Wind); 
                Total_Lon_Wind_fall = vertcat(Total_Lon_Wind_fall, Lon_Wind); 
                
            end
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR 
        
    end
end

%%

[CMOD_one_degree_winter, CMOD_OCC_winter, CMOD_STD_winter, CMOD_ERR_winter]  = hist_wt_occ_tot(Total_Lat_CMOD_winter, Total_Lon_CMOD_winter, Total_CMOD_winter, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_spring, CMOD_OCC_spring, CMOD_STD_spring, CMOD_ERR_spring]  = hist_wt_occ_tot(Total_Lat_CMOD_spring, Total_Lon_CMOD_spring, Total_CMOD_spring, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_summer, CMOD_OCC_summer, CMOD_STD_summer, CMOD_ERR_summer]  = hist_wt_occ_tot(Total_Lat_CMOD_summer, Total_Lon_CMOD_summer, Total_CMOD_summer, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_fall, CMOD_OCC_fall, CMOD_STD_fall, CMOD_ERR_fall]          = hist_wt_occ_tot(Total_Lat_CMOD_fall, Total_Lon_CMOD_fall, Total_CMOD_fall, Aqua_Lat', Aqua_Lon');



[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, Aqua_Lat', Aqua_Lon');




[Wind_one_degree_winter, Wind_OCC_winter, Wind_STD_winter, Wind_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Wind_winter, Total_Lon_Wind_winter, Total_Wind_winter, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_spring, Wind_OCC_spring, Wind_STD_spring, Wind_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Wind_spring, Total_Lon_Wind_spring, Total_Wind_spring, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_summer, Wind_OCC_summer, Wind_STD_summer, Wind_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Wind_summer, Total_Lon_Wind_summer, Total_Wind_summer, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_fall, Wind_OCC_fall, Wind_STD_fall, Wind_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Wind_fall, Total_Lon_Wind_fall, Total_Wind_fall, Aqua_Lat', Aqua_Lon');



%%
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Climatology_Figures
% here is where I'm saving my figures.

%%
winter = Wind_one_degree_winter; 
spring = Wind_one_degree_spring; 
summer = Wind_one_degree_summer;
fall   =  Wind_one_degree_fall;

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure,clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Fall'); 
    
subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, sprintf('m s^{-1}'), 'FontSize', 18);
% hy.FontSize = 18;

%%

winter = Ice_one_degree_winter; 
spring = Ice_one_degree_spring; 
summer = Ice_one_degree_summer;
fall   = Ice_one_degree_fall;


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure' clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    step_Lon,...
    step_Lat,...
    'Ice', ...
    [0 0.8],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    step_Lon,...
    step_Lat,...
    'Ice', ...
    [0 0.8],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    step_Lon,...
    step_Lat,...
    'Ice', ...
    [0 0.8],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    step_Lon,...
    step_Lat,...
    'Ice', ...
    [0 0.8],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'Depolarization Ratio \delta', 'FontSize', 16);
% hy.FontSize = 18;


%%
% Smoothn Climatologies 

CMOD_winter_smoothn = smoothn(CMOD_one_degree_winter);
CMOD_spring_smoothn = smoothn(CMOD_one_degree_spring); 
CMOD_summer_smoothn = smoothn(CMOD_one_degree_summer); 
CMOD_fall_smoothn = smoothn(CMOD_one_degree_fall); 


%%
winter = CMOD_winter_smoothn; 
spring = CMOD_spring_smoothn; 
summer = CMOD_summer_smoothn;
fall   = CMOD_fall_smoothn;

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    step_Lon,...
    step_Lat,...
    'tempo', ...
    [0 .08],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    step_Lon,...
    step_Lat,...
    'tempo', ...
    [0 .08],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    step_Lon,...
    step_Lat,...
    'tempo', ...
    [0 .08],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    step_Lon,...
    step_Lat,...
    'tempo', ...
    [0 .08],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'CMOD', 'FontSize', 18);
% hy.FontSize = 18;




%%
Wind_winter_smoothn = smoothn(Wind_one_degree_winter);
Wind_spring_smoothn = smoothn(Wind_one_degree_spring); 
Wind_summer_smoothn = smoothn(Wind_one_degree_summer); 
Wind_fall_smoothn = smoothn(Wind_one_degree_fall); 

%%
winter = Wind_winter_smoothn; 
spring = Wind_spring_smoothn; 
summer = Wind_summer_smoothn;
fall   = Wind_fall_smoothn;
% 

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure,clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    step_Lon,...
    step_Lat,...
    'amp', ...
    [0 20],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, sprintf('m s^{-1}'), 'FontSize', 18);
% hy.FontSize = 18;

%% Time Series Climatologies...


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

load('Depol_Ratio_Monthly_avg_Vars.mat')
load('amsrmf_Monthly_avg_Vars.mat')
load('CMOD_Monthly_avg_Vars.mat') 


cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Climatology_TimeSeries

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 



t1 = datetime(2007,01,01);
t2 = datetime(2007,12,31);
times_months = t1:calmonths(1):t2; 
times_months_num = month(times_months);

times = times'; 
Ice_Vector = double(Depol_Ratio_Monthly_avg(:,1));
Ice_Vector(117,:) = [];
CMOD_Vector = CMOD_Monthly_avg(:,1);
CMOD_Vector(117,:) = [];
Wind_Vector = amsrmf_Monthly_avg;
times_without_feb = times';
times_without_feb(117,:) = []; %
times_for_wind = times';
times_for_wind(117,:) =[];
times_for_wind(64:73,:) = [];
Wind_Vector(117,:) = [];
Wind_Vector(64:73,:)=[];



[Ac_CMOD, tc_CMOD] = climatology(CMOD_Vector,times_without_feb, 'monthly'); 
[Ac_Ice, tc_Ice] = climatology(Ice_Vector, times_without_feb, 'monthly'); 
[Ac_Wind, tc_Wind] = climatology(Wind_Vector, times_for_wind, 'monthly'); 

Master_chl_a(isnan(Master_chl_a)) = 0; 
[Ac_Chl, tc_Chl] = climatology(Master_chl_a, times, 'monthly'); 

for i = 1:12
Ac_Chl_monthly(i) = mean2(Ac_Chl(:,:,i));
end

Ac_Chl_monthly(Ac_Chl_monthly <= 0) = NaN;

%%
sienna = rgb('sienna'); 
gray = rgb('gray'); 
black = rgb('black'); 
grass_green = rgb('grass green');
ice_blue = rgb('sky blue'); 
wind_blue = rgb('royal blue');

% Create textbox
% pos = [x-start y-start x-width y-height]


%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.05], [0.1 0.05], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end


figure(1), clf;

ax(1) = subplot(4,1,1);
grid on
aa_splot(1:12, Ac_CMOD, 'x-',...
    'linewidth', 1.5, ...
    'Color', black,...
    'MarkerSize', 9,...
    'MarkerFaceColor', black,...
    'MarkerEdgeColor', black);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% 
% title('CMOD Climatology')
ylabel('CMOD')
%         xlabel('Months')
xlim([1,12])
ylim([0.01 0.1])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 16)
%                   xtickangle(20)


ax(2) = subplot(4,1,2);
grid on
aa_splot(1:12, Ac_Chl_monthly, 'v-',...
    'linewidth', 1.5, ...
    'Color', grass_green,...
    'MarkerSize', 7,...
    'MarkerFaceColor', grass_green,...
    'MarkerEdgeColor', grass_green);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Chl-{\ita} Climatology')
ylabel(sprintf('mg m^{-3}')),...
    %         xlabel('Months')
xlim([1,12])
ylim([0 0.35])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 16)
% xtickangle(45)




ax(3) = subplot(4,1,3); 
grid on
aa_splot(1:12, Ac_Ice, 'o-',...
    'linewidth', 1.5, ...
    'Color', ice_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', ice_blue,...
    'MarkerEdgeColor', ice_blue);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Ice Climatology')
ylabel('Depolarization Ratio \delta',...
   'FontName','Helvetica Neue');%         xlabel('Months')
xlim([1,12])
ylim([0.4 0.7])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 16)
% xtickangle(45)




ax(4) = subplot(4,1,4);
grid on
aa_splot(1:12, Ac_Wind, 'd-',...
    'linewidth', 1.5, ...
    'Color', wind_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', wind_blue,...
    'MarkerEdgeColor', wind_blue);
set(gca, 'xtick', 1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Wind Speed Climatology')
ylabel(sprintf('m s^{-1}'))
xlim([1,12])
%          ylim([0 20])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 16)
 xtickangle(45)


[ax(1:3).XTickLabel] = deal([]);
ax(1).YTick(1) = []; 
ax(2).YTick(1) = []; 
ax(3).YTick(1) = []; 
ax(4).YTick(1) = [];
ax(4).FontSize = 18;
% ax(1) to see properties 

% from below: an(1).Position(2) = 0.87
% ax(1).Position(2) = 0.7673

% (0.7673 - should give the correct position for first subplot 

an(1) = annotation(gcf,'textbox',... % I drew this on in the figure 
  [0.11729695024077 0.900712682379355 0.16 0.0349],...
    'String','CMOD Climatology',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(2) = annotation(gcf,'textbox',...
    [an(1).Position(1),...
    (ax(2).Position(2) + (an(1).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
    0.16 0.0349],...
    'String',{'Chl-{\ita} Climatology'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(3) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(3).Position(2) + (an(2).Position(2) - ax(2).Position(2))),...
    0.16 0.0349],...
    'String',{'Ice Climatology'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(4) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(4).Position(2) + (an(3).Position(2) - ax(3).Position(2))),...
    0.16 0.0349],...
    'String',{'Wind Climatology'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end
%% Upper and Lower Subplots with Titles
income = [3.2,4.1,5.0,5.6];
outgo = [2.5,4.0,3.35,4.9];
subplot(2,1,1); plot(income)
title('Income')
subplot(2,1,2); plot(outgo)
title('Outgo')
%% Subplots in Quadrants
figure
subplot(2,2,1)
text(.5,.5,{'subplot(2,2,1)';'or subplot 221'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,2)
text(.5,.5,{'subplot(2,2,2)';'or subplot 222'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,3)
text(.5,.5,{'subplot(2,2,3)';'or subplot 223'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,{'subplot(2,2,4)';'or subplot 224'},...
    'FontSize',14,'HorizontalAlignment','center')



%%






%%%%%%%%%%%%%% The rest of this script is not needed for now... Disregard all lines below %%%%%%%%%%%%%%%%%%%%%%









%%

%%%%%%%%% Saving CMOD, Ice, and Wind Data as 3D arrays with third dimension
%%%%%%%%% as each individual month from 2006 to 2018 

cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY

load('Aqua_Longitude_Subset_BellingshausenSea.mat')
load('Aqua_Latitude_Subset_BellingshausenSea.mat')


%%

cd '/Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Spatial_Correlation_Analysis_matfiles_2_03_2020/Across_fullrange'

load('Total_timetable_winds.mat')
load('Total_timetable_SeaIce.mat')
load('Total_timetable_CMOD_Night.mat')


%%


%%
   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
    
    %TR_JAN_2006 = timerange('2006-01-01', '2006-02-01'); there is no CALIPSO DATA until June of 2006
    TR_JAN_2007 = timerange('2007-01-01', '2007-02-01'); 
    TR_JAN_2008 = timerange('2008-01-01', '2008-02-01');
    TR_JAN_2009 = timerange('2009-01-01', '2009-02-01');
    TR_JAN_2010 = timerange('2010-01-01', '2010-02-01'); 
    TR_JAN_2011 = timerange('2011-01-01', '2011-02-01');
    TR_JAN_2012 = timerange('2012-01-01', '2012-02-01');
    TR_JAN_2013 = timerange('2013-01-01', '2013-02-01');
    TR_JAN_2014 = timerange('2014-01-01', '2014-02-01');
    TR_JAN_2015 = timerange('2015-01-01', '2015-02-01');
    TR_JAN_2016 = timerange('2016-01-01', '2016-02-01');
    TR_JAN_2017 = timerange('2017-01-01', '2017-02-01');
    TR_JAN_2018 = timerange('2018-01-01', '2018-02-01'); 
    
    %TR_FEB_2006 = timerange('2006-02-01', '2006-03-01'); 
    TR_FEB_2007 = timerange('2007-02-01', '2007-03-01');
    TR_FEB_2008 = timerange('2008-02-01', '2008-03-01'); 
    TR_FEB_2009 = timerange('2009-02-01', '2009-03-01'); 
    TR_FEB_2010 = timerange('2010-02-01', '2010-03-01'); 
    TR_FEB_2011 = timerange('2011-02-01', '2011-03-01');
    TR_FEB_2012 = timerange('2012-02-01', '2012-03-01'); 
    TR_FEB_2013 = timerange('2013-02-01', '2013-03-01'); 
    TR_FEB_2014 = timerange('2014-02-01', '2014-03-01');
    TR_FEB_2015 = timerange('2015-02-01', '2015-03-01'); 
    TR_FEB_2016 = timerange('2016-02-01', '2016-03-01'); 
    TR_FEB_2017 = timerange('2017-02-01', '2017-03-01');
    TR_FEB_2018 = timerange('2018-02-01', '2018-03-01'); 
    
    %TR_MAR_2006 = timerange('2006-03-01', '2006-04-01'); 
    TR_MAR_2007 = timerange('2007-03-01', '2007-04-01');
    TR_MAR_2008 = timerange('2008-03-01', '2008-04-01'); 
    TR_MAR_2009 = timerange('2009-03-01', '2009-04-01'); 
    TR_MAR_2010 = timerange('2010-03-01', '2010-04-01'); 
    TR_MAR_2011 = timerange('2011-03-01', '2011-04-01');
    TR_MAR_2012 = timerange('2012-03-01', '2012-04-01'); 
    TR_MAR_2013 = timerange('2013-03-01', '2013-04-01'); 
    TR_MAR_2014 = timerange('2014-03-01', '2014-04-01');
    TR_MAR_2015 = timerange('2015-03-01', '2015-04-01'); 
    TR_MAR_2016 = timerange('2016-03-01', '2016-04-01'); 
    TR_MAR_2017 = timerange('2017-03-01', '2017-04-01');
    TR_MAR_2018 = timerange('2018-03-01', '2018-04-01'); 
    
    %TR_APR_2006   = timerange('2006-03-01', '2006-01-01');
    TR_APR_2007   = timerange('2007-04-01', '2007-05-01'); 
    TR_APR_2008   = timerange('2008-04-01', '2008-05-01');
    TR_APR_2009   = timerange('2009-04-01', '2009-05-01');
    TR_APR_2010   = timerange('2010-04-01', '2010-05-01');
    TR_APR_2011   = timerange('2011-04-01', '2011-05-01');
    TR_APR_2012   = timerange('2012-04-01', '2012-05-01');
    TR_APR_2013   = timerange('2013-04-01', '2013-05-01');
    TR_APR_2014   = timerange('2014-04-01', '2014-05-01');
    TR_APR_2015   = timerange('2015-04-01', '2015-05-01');
    TR_APR_2016   = timerange('2016-04-01', '2016-05-01');
    TR_APR_2017   = timerange('2017-04-01', '2017-05-01');
    TR_APR_2018   = timerange('2018-04-01', '2018-05-01');
    
        %TR_MAY_2006 = timerange('2006-05-01', '2006-06-01'); 
    TR_MAY_2007 = timerange('2007-05-01', '2007-06-01'); 
    TR_MAY_2008 = timerange('2008-05-01', '2008-06-01');
    TR_MAY_2009 = timerange('2009-05-01', '2009-06-01');
    TR_MAY_2010 = timerange('2010-05-01', '2010-06-01'); 
    TR_MAY_2011 = timerange('2011-05-01', '2011-06-01');
    TR_MAY_2012 = timerange('2012-05-01', '2012-06-01');
    TR_MAY_2013 = timerange('2013-05-01', '2013-06-01');
    TR_MAY_2014 = timerange('2014-05-01', '2014-06-01');
    TR_MAY_2015 = timerange('2015-05-01', '2015-06-01');
    TR_MAY_2016 = timerange('2016-05-01', '2016-06-01');
    TR_MAY_2017 = timerange('2017-05-01', '2017-06-01');
    TR_MAY_2018 = timerange('2018-05-01', '2018-06-01'); 
    
    TR_JUN_2006 = timerange('2006-06-01', '2006-07-01'); 
    TR_JUN_2007 = timerange('2007-06-01', '2007-07-01'); 
    TR_JUN_2008 = timerange('2008-06-01', '2008-07-01');
    TR_JUN_2009 = timerange('2009-06-01', '2009-07-01');
    TR_JUN_2010 = timerange('2010-06-01', '2010-07-01'); 
    TR_JUN_2011 = timerange('2011-06-01', '2011-07-01');
    TR_JUN_2012 = timerange('2012-06-01', '2012-07-01');
    TR_JUN_2013 = timerange('2013-06-01', '2013-07-01');
    TR_JUN_2014 = timerange('2014-06-01', '2014-07-01');
    TR_JUN_2015 = timerange('2015-06-01', '2015-07-01');
    TR_JUN_2016 = timerange('2016-06-01', '2016-07-01');
    TR_JUN_2017 = timerange('2017-06-01', '2017-07-01');
    TR_JUN_2018 = timerange('2018-06-01', '2018-07-01'); 
    
    TR_JUL_2006 = timerange('2006-07-01', '2006-08-01'); 
    TR_JUL_2007 = timerange('2007-07-01', '2007-08-01'); 
    TR_JUL_2008 = timerange('2008-07-01', '2008-08-01');
    TR_JUL_2009 = timerange('2009-07-01', '2009-08-01');
    TR_JUL_2010 = timerange('2010-07-01', '2010-08-01'); 
    TR_JUL_2011 = timerange('2011-07-01', '2011-08-01');
    TR_JUL_2012 = timerange('2012-07-01', '2012-08-01');
    TR_JUL_2013 = timerange('2013-07-01', '2013-08-01');
    TR_JUL_2014 = timerange('2014-07-01', '2014-08-01');
    TR_JUL_2015 = timerange('2015-07-01', '2015-08-01');
    TR_JUL_2016 = timerange('2016-07-01', '2016-08-01');
    TR_JUL_2017 = timerange('2017-07-01', '2017-08-01');
    TR_JUL_2018 = timerange('2018-07-01', '2018-08-01'); 
    
    TR_AUG_2006 = timerange('2006-08-01', '2006-09-01'); 
    TR_AUG_2007 = timerange('2007-08-01', '2007-09-01'); 
    TR_AUG_2008 = timerange('2008-08-01', '2008-09-01');
    TR_AUG_2009 = timerange('2009-08-01', '2009-09-01');
    TR_AUG_2010 = timerange('2010-08-01', '2010-09-01'); 
    TR_AUG_2011 = timerange('2011-08-01', '2011-09-01');
    TR_AUG_2012 = timerange('2012-08-01', '2012-09-01');
    TR_AUG_2013 = timerange('2013-08-01', '2013-09-01');
    TR_AUG_2014 = timerange('2014-08-01', '2014-09-01');
    TR_AUG_2015 = timerange('2015-08-01', '2015-09-01');
    TR_AUG_2016 = timerange('2016-08-01', '2016-09-01');
    TR_AUG_2017 = timerange('2017-08-01', '2017-09-01');
    TR_AUG_2018 = timerange('2018-08-01', '2018-09-01'); 
    
    TR_SEP_2006 = timerange('2006-09-01', '2006-10-01'); 
    TR_SEP_2007 = timerange('2007-09-01', '2007-10-01'); 
    TR_SEP_2008 = timerange('2008-09-01', '2008-10-01');
    TR_SEP_2009 = timerange('2009-09-01', '2009-10-01');
    TR_SEP_2010 = timerange('2010-09-01', '2010-10-01'); 
    TR_SEP_2011 = timerange('2011-09-01', '2011-10-01');
    TR_SEP_2012 = timerange('2012-09-01', '2012-10-01');
    TR_SEP_2013 = timerange('2013-09-01', '2013-10-01');
    TR_SEP_2014 = timerange('2014-09-01', '2014-10-01');
    TR_SEP_2015 = timerange('2015-09-01', '2015-10-01');
    TR_SEP_2016 = timerange('2016-09-01', '2016-10-01');
    TR_SEP_2017 = timerange('2017-09-01', '2017-10-01');
    TR_SEP_2018 = timerange('2018-09-01', '2018-10-01'); 
    
    TR_OCT_2006 = timerange('2006-10-01', '2006-11-01'); 
    TR_OCT_2007 = timerange('2007-10-01', '2007-11-01'); 
    TR_OCT_2008 = timerange('2008-10-01', '2008-11-01');
    TR_OCT_2009 = timerange('2009-10-01', '2009-11-01');
    TR_OCT_2010 = timerange('2010-10-01', '2010-11-01'); 
    TR_OCT_2011 = timerange('2011-10-01', '2011-11-01');
    TR_OCT_2012 = timerange('2012-10-01', '2012-11-01');
    TR_OCT_2013 = timerange('2013-10-01', '2013-11-01');
    TR_OCT_2014 = timerange('2014-10-01', '2014-11-01');
    TR_OCT_2015 = timerange('2015-10-01', '2015-11-01');
    TR_OCT_2016 = timerange('2016-10-01', '2016-11-01');
    TR_OCT_2017 = timerange('2017-10-01', '2017-11-01');
    TR_OCT_2018 = timerange('2018-10-01', '2018-11-01'); 
    
    TR_NOV_2006 = timerange('2006-11-01', '2006-12-01'); 
    TR_NOV_2007 = timerange('2007-11-01', '2007-12-01'); 
    TR_NOV_2008 = timerange('2008-11-01', '2008-12-01');
    TR_NOV_2009 = timerange('2009-11-01', '2009-12-01');
    TR_NOV_2010 = timerange('2010-11-01', '2010-12-01'); 
    TR_NOV_2011 = timerange('2011-11-01', '2011-12-01');
    TR_NOV_2012 = timerange('2012-11-01', '2012-12-01');
    TR_NOV_2013 = timerange('2013-11-01', '2013-12-01');
    TR_NOV_2014 = timerange('2014-11-01', '2014-12-01');
    TR_NOV_2015 = timerange('2015-11-01', '2015-12-01');
    TR_NOV_2016 = timerange('2016-11-01', '2016-12-01');
    TR_NOV_2017 = timerange('2017-11-01', '2017-12-01');
    TR_NOV_2018 = timerange('2018-11-01', '2018-12-01'); 
    
    TR_DEC_2006 = timerange('2006-12-01', '2007-01-01'); 
    TR_DEC_2007 = timerange('2007-12-01', '2008-01-01'); 
    TR_DEC_2008 = timerange('2008-12-01', '2009-01-01');
    TR_DEC_2009 = timerange('2009-12-01', '2010-01-01');
    TR_DEC_2010 = timerange('2010-12-01', '2011-01-01'); 
    TR_DEC_2011 = timerange('2011-12-01', '2012-01-01');
    TR_DEC_2012 = timerange('2012-12-01', '2013-01-01');
    TR_DEC_2013 = timerange('2013-12-01', '2014-01-01');
    TR_DEC_2014 = timerange('2014-12-01', '2015-01-01');
    TR_DEC_2015 = timerange('2015-12-01', '2016-01-01');
    TR_DEC_2016 = timerange('2016-12-01', '2017-01-01');
    TR_DEC_2017 = timerange('2017-12-01', '2018-01-01');
    TR_DEC_2018 = timerange('2018-12-01', '2019-01-01'); 
   


%%
    
% cmap_night = cmocean('gray');

% night Winter

%%
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Interpolation_of_Data/New_Resolution_1degree_Interpolated_Monthly/Smoothn_on_2D_Slices/Smoothn_Climatologies




%%
clear CMOD_MONTHS
clear Winds_MONTHS
clear SeaIce_MONTHS
clear Ice 


Aqua_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
Aqua_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 

Aqua_Lat = flip(Aqua_Lat) ; 


% CMOD_MONTHS = zeros(15,23,151); 
% Winds_MONTHS = zeros(15,23,151); 
% SeaIce_MONTHS = zeros(15,23,151);

% 7 months in 2006, 12 months in the other, so third dimension on monthly
% subset should be 7 + 12(12) = 151

%%

% 2006 : 2010, third dimension should be 7 + 12(4) = 55

for i = 2006:2018
    
    disp(i)
    
    MONTH = {'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'};
    
    for j = 1:length(MONTH)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_JAN_occ', i)),...
        %     eval(sprintf('night_%d_JAN_std')), eval(sprintf('night_%d_JAN_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', MONTH{j}, i),'var') % only keep looping if the variable actually exists
            % TR_APR_2006 does not exist, which is why I need this
            % statement in here
            
            eval(sprintf('Lat_CMOD = Total_timetable_CMOD_Night(TR_%s_%d,:).Total_Latitude_Night_Cloud_Free;', MONTH{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_CMOD_Night(TR_%s_%d, :).Total_Longitude_Night_Cloud_Free;',MONTH{j}, i))
            eval(sprintf('OD       = Total_timetable_CMOD_Night(TR_%s_%d, :).Total_OD_Night_Cloud_Free;', MONTH{j}, i))
            eval(sprintf('Lat_Ice  = Total_timetable_SeaIce(TR_%s_%d,:).Total_Latitude;', MONTH{j}, i))
            eval(sprintf('Lon_Ice  = Total_timetable_SeaIce(TR_%s_%d,:).Total_Longitude;', MONTH{j}, i))
            eval(sprintf('Ice      = Total_timetable_SeaIce(TR_%s_%d,:).Total_SeaIce;', MONTH{j},i))
            eval(sprintf('Lat_Wind = Total_timetable_winds(TR_%s_%d,:).Total_Latitude_wind;', MONTH{j}, i))
            eval(sprintf('Lon_Wind = Total_timetable_winds(TR_%s_%d, :).Total_Longitude_wind;', MONTH{j},i))
            eval(sprintf('amsrMF   = Total_timetable_winds(TR_%s_%d,:).Total_windamsrMF_nonan;', MONTH{j},i))
            
           
            bad_Ice_values = Ice <= -0.2 | Ice > 1.2; 
            Ice(bad_Ice_values) = NaN;
            
            nan_ice        = isnan(Ice(:,1)); 
            Ice_nonan      = Ice(~nan_ice) ;
            Lat_Ice_nonan  = Lat_Ice(~nan_ice); 
            Lon_Ice_nonan  = Lon_Ice(~nan_ice); 
            % need to start with SEASON{j} because you start with %s
            % above
            
            [CMOD, CMOD_OCC, CMOD_STD, CMOD_ERR]         = hist_wt_occ_tot(Lat_CMOD, Lon_CMOD, OD, Aqua_Lat', Aqua_Lon');
            [SeaIce, SeaIce_OCC, SeaIce_STD, SeaIce_ERR] = hist_wt_occ_tot(Lat_Ice_nonan, Lon_Ice_nonan, Ice_nonan, Aqua_Lat', Aqua_Lon');
            [Winds, Winds_OCC, Winds_STD, Winds_ERR]     = hist_wt_occ_tot(Lat_Wind, Lon_Wind, amsrMF, Aqua_Lat', Aqua_Lon');
            
            % CMOD 3D Monthly Array Now
            if ~exist('CMOD_MONTHS', 'var')
                
                CMOD_MONTHS = CMOD;
                SeaIce_MONTHS = SeaIce;
                Winds_MONTHS = Winds;
                
            else
                
                CMOD_MONTHS(:,:, end + 1) = cat(3, CMOD);
                SeaIce_MONTHS(:,:, end + 1) = cat(3, SeaIce);
                Winds_MONTHS(:,:, end + 1) = cat(3, Winds);
                
                
            end
            
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR SeaIce SeaIce_OCC SeaIce_STD SeaIce_ERR Winds Winds_OCC Winds_STD Winds_ERR Lat_Ice_nonan Lon_Ice_nonan Ice_nonan
        
    end
end

%% 

CMOD_MONTHS_one_degree_res = CMOD_MONTHS; 
SeaIce_MONTHS_one_degree_res = SeaIce_MONTHS;
Winds_MONTHS_one_degree_res = Winds_MONTHS;

%%

% save('CMOD_MONTHS_one_degree_res.mat', 'CMOD_MONTHS_one_degree_res') 
% save('SeaIce_MONTHS_one_degree_res.mat', 'SeaIce_MONTHS_one_degree_res') 
% save('Winds_MONTHS_one_degree_res.mat', 'Winds_MONTHS_one_degree_res')

%%

CMOD_one_degree_MONTHS_interpolated_inpaintn = inpaintn(CMOD_MONTHS_one_degree_res);
CMOD_one_degree_MONTHS_interpolated_smoothn = smoothn(CMOD_MONTHS_one_degree_res, 'robust'); 

SeaIce_one_degree_interpolated_inpaintn = inpaintn(SeaIce_MONTHS_one_degree_res); 
SeaIce_one_degree_interpolated_smoothn = smoothn(SeaIce_MONTHS_one_degree_res); 

Winds_one_degree_interpolated_inpaintn = inpaintn(Winds_MONTHS_one_degree_res); 
Winds_one_degree_interpolated_smoothn = smoothn(Winds_MONTHS_one_degree_res); 


%%

save('CMOD_one_degree_MONTHS_interpolated_inpaintn.mat', 'CMOD_one_degree_MONTHS_interpolated_inpaintn')
save('CMOD_one_degree_MONTHS_interpolated_smoothn.mat', 'CMOD_one_degree_MONTHS_interpolated_smoothn') 

save('SeaIce_one_degree_interpolated_inpaintn.mat', 'SeaIce_one_degree_interpolated_inpaintn') 
save('SeaIce_one_degree_interpolated_smoothn.mat', 'SeaIce_one_degree_interpolated_smoothn') 

save('Winds_one_degree_interpolated_inpaintn.mat', 'Winds_one_degree_interpolated_inpaintn') 
save('Winds_one_degree_interpolated_smoothn.mat', 'Winds_one_degree_interpolated_smoothn') 

%% 

times = datetime(2006, 06, 01) : calmonths(1) : datetime(2018, 12, 31);

% Now let's make those movies... 

figure(1), clf

vidfile = VideoWriter('CMOD_Months_one_degree_rawdata','MPEG-4');
vidfile.FrameRate = 2;
open(vidfile);

for i = 1 : 2 : 151
   
    disp(i)
    
    plot_BSea_figure_climatology(CMOD_MONTHS_one_degree_res(:,:,i),...
        Aqua_Lon,...
        Aqua_Lat,...
        'tempo',...
        [0 0.1],...
        'CMOD',...
        (['CMOD Monthly Scale, one degree resolution, ', datestr(times(i))]))
   
    drawnow 
    
    
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));

end

close(vidfile)



%%

figure(2), clf

vidfile = VideoWriter('CMOD_Months_one_degree_inpaintn','MPEG-4');
vidfile.FrameRate = 2;
open(vidfile);

for i = 1 : 2 : 151
   
    disp(i)
    
    plot_BSea_figure_climatology(CMOD_one_degree_MONTHS_interpolated_inpaintn(:,:,i),...
        Aqua_Lon,...
        Aqua_Lat,...
        'tempo',...
        [0 0.1],...
        'CMOD',...
        (['CMOD Monthly Scale, inpaintn interpolation, one degree resolution, ', datestr(times(i))]))
   
    drawnow 
    
    
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));

end

close(vidfile)

%%

figure(3), clf

vidfile = VideoWriter('CMOD_Months_one_degree_smoothn','MPEG-4');
vidfile.FrameRate = 2;
open(vidfile);

for i = 1 : 2 : 151
   
    disp(i)
    
    plot_BSea_figure_climatology(CMOD_one_degree_MONTHS_interpolated_smoothn(:,:,i),...
        Aqua_Lon,...
        Aqua_Lat,...
        'tempo',...
        [0 0.1],...
        'CMOD',...
        (['CMOD Monthly Scale, inpaintn interpolation, one degree resolution, ', datestr(times(i))]))
   
    drawnow 
    
    
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));

end

close(vidfile)


%%









