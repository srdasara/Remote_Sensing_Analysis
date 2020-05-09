%% Read CMOD, Construct into timetable, save out in 
% /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/CALIPSO_LID_L2_05kmAPro-Standard-V4-20/UPDATED_Calipso_2006_2018_CMOD_Night_CloudFree_AdjustedAltRange

% can these push these into monthly averaged time series
% climatologies
% spearman's rho 

%% 0.0977KM TO 2.0137KM 
% these correspond to altitude bins 358 to 390. 

% first load these
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/CALIPSO_LID_L2_05kmAPro-Standard-V4-20/Original_Variables_BellingshausenSea

load('Total_EC_532.mat')
load('Total_Day_Night_Flag.mat')
load('Total_COD_Cloud.mat')
load('Total_altitudes.mat')
load('Total_Longitude.mat')
load('Total_Latitude.mat') 
load('Total_Profile_Time_New.mat')
load('Total_Surface_532_Integrated_Depolarization_Ratio.mat')
load('Total_Profile_Time.mat') 
load('Total_windamsrMF.mat') 
   
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication


%% Filter to only include Extinction Coefficients that are cloud free and nighttime only. 
% Total_COD_Cloud == 0 , DayNightFlag == 1

Cloud_Free                             = Total_COD_Cloud(:, 1) == 0; 
Total_Latitude_Cloud_Free              = Total_Latitude(Cloud_Free); 
Total_Longitude_Cloud_Free             = Total_Longitude(Cloud_Free); 
Total_Profile_Time_New_Cloud_Free      = Total_Profile_Time_New(Cloud_Free);
Total_EC_532_Cloud_Free                = Total_EC_532(Cloud_Free, :) ; 
Total_Day_Night_Flag_Cloud_Free        = Total_Day_Night_Flag(Cloud_Free,:); 
Total_EC_532_Cloud_Free_adjusted_alt   = Total_EC_532_Cloud_Free(:, 358:390) ; 
Total_adjusted_alt                     = Total_altitudes(358:390, :) ; 



Night_Cloud_Free                        = Total_Day_Night_Flag_Cloud_Free(:,1) == 1 ; 
% Total_Night_Flag_Cloud_Free           = Total_Day_Night_Flag_Cloud_Free(Night_Cloud_Free); 
Total_Latitude_Night_Cloud_Free         = Total_Latitude_Cloud_Free(Night_Cloud_Free); 
Total_Longitude_Night_Cloud_Free        = Total_Longitude_Cloud_Free(Night_Cloud_Free); 
Total_Profile_Time_New_Night_Cloud_Free = Total_Profile_Time_New_Cloud_Free(Night_Cloud_Free); 
Total_EC_532_Night_Cloud_Free           = Total_EC_532_Cloud_Free_adjusted_alt(Night_Cloud_Free, :); 

     
%%                              
    
% Here, my CMOD variable is cloud-free with only nighttime profiles used. 

 Total_EC_532_Night_Cloud_Free(isnan(Total_EC_532_Night_Cloud_Free)) = 0 ; % 0 is clear air, NaN has been filtered out by quality screening
%  Total_EC_532_Night_Cloud_Free(Total_EC_532_Night_Cloud_Free == 0) = NaN ; % Converting all zeros in sigma to NaNs. 

 % to keep NaN or not keep NaN?
 %%
clear CMOD

CMOD = zeros(length(Total_EC_532_Night_Cloud_Free(:,1)), 1); 

for i = 1:length(Total_EC_532_Night_Cloud_Free(:,1))
%     disp(i)
    CMOD(i) = -1 .* (trapz(Total_alt_cloud_free, Total_EC_532_Night_Cloud_Free(i,:))) ; 
    % -1 in equation above was to flip in consideration of the fact that altitudes start from 2.0137 km and end at 0.0977km
end

save('CMOD.mat', 'CMOD') ; 

    %% Make a time table with all of these values, 3 separate ones: CMOD, Ice, and Winds
    
    Total_table_CMOD = table(Total_Profile_Time_New_Night_Cloud_Free,...
        Total_Latitude_Night_Cloud_Free,...
        Total_Longitude_Night_Cloud_Free,...
        CMOD); 

    Total_table_CMOD             = sortrows(Total_table_CMOD,'Total_Profile_Time_New_Night_Cloud_Free','ascend'); % sort values with increasing time duration
    Total_timetable_CMOD         = table2timetable(Total_table_CMOD); % make table into a timetable
save('Total_timetable_CMOD.mat', 'Total_timetable_CMOD', '-v7.3')
    %% I first have to filter for only good values of Ice 
    
    bad_Ice_values = Total_Surface_532_Integrated_Depolarization_Ratio <= -0.2 | Total_Surface_532_Integrated_Depolarization_Ratio > 1.2;
    Total_Surface_532_Integrated_Depolarization_Ratio(bad_Ice_values) = NaN; % I set these bad values to NaNs so I can easily index and remove them
    
    nan_ice        = isnan(Total_Surface_532_Integrated_Depolarization_Ratio(:,1));
    Total_Surface_532_Integrated_Depolarization_Ratio      = Total_Surface_532_Integrated_Depolarization_Ratio(~nan_ice) ;
    Total_Latitude_Ice  = Total_Latitude(~nan_ice);
    Total_Longitude_Ice  = Total_Longitude(~nan_ice);
    Total_Profile_Time_New_Ice = Total_Profile_Time_New(~nan_ice); 
    
    Total_table_Depol_Ratio = table(Total_Profile_Time_New_Ice,...
        Total_Latitude_Ice,...
        Total_Longitude_Ice,...
        Total_Surface_532_Integrated_Depolarization_Ratio);
    
    Total_table_Depol_Ratio = sortrows(Total_table_Depol_Ratio, 'Total_Profile_Time_New_Ice', 'ascend');
    
    Total_timetable_Depol_Ratio = table2timetable(Total_table_Depol_Ratio);
    save('Total_timetable_Depol_Ratio.mat', 'Total_timetable_Depol_Ratio', '-v7.3') 
    
    %%
    
    bad_Wind_values = Total_windamsrMF <= 0| Total_windamsrMF > 50;
    Total_windamsrMF(bad_Wind_values) = NaN; % I set these bad values to NaNs so I can easily index and remove them

    nan_wind        = isnan(Total_windamsrMF(:,1));
    Total_windamsrMF      = Total_windamsrMF(~nan_wind) ;
    Total_Latitude_Wind  = Total_Latitude(~nan_wind);
    Total_Longitude_Wind  = Total_Longitude(~nan_wind);
    Total_Profile_Time_New_Wind = Total_Profile_Time_New(~nan_wind); 
    
    
    Total_table_amsrmf = table(Total_Profile_Time_New_Wind,... 
        Total_Latitude_Wind,...
        Total_Longitude_Wind,... 
        Total_windamsrMF);
    
    Total_table_amsrmf = sortrows(Total_table_amsrmf, 'Total_Profile_Time_New_Wind', 'ascend'); 
    Total_timetable_amsrmf = table2timetable(Total_table_amsrmf); 
    save('Total_timetable_amsrmf.mat', 'Total_timetable_amsrmf', '-v7.3') 
    
    %%
    
    timetable_CMOD_monthly_avg    = retime(Total_timetable_CMOD, 'monthly', @nanmean); 
    CMOD_Monthly_avg = timetable_CMOD_monthly_avg.CMOD; 
    CMOD_Time_Months = timetable_CMOD_monthly_avg.Total_Profile_Time_New_Night_Cloud_Free; 
    CMOD_Lat_Months  = timetable_CMOD_monthly_avg.Total_Latitude_Night_Cloud_Free;
    CMOD_Lon_Months  = timetable_CMOD_monthly_avg.Total_Longitude_Night_Cloud_Free; 
    
    timetable_Depol_Ratio_monthly_avg = retime(Total_timetable_Depol_Ratio, 'monthly', @nanmean); 
    Depol_Ratio_Monthly_avg = timetable_Depol_Ratio_monthly_avg.Total_Surface_532_Integrated_Depolarization_Ratio;
    Depol_Ratio_Time_Months = timetable_Depol_Ratio_monthly_avg.Total_Profile_Time_New_Ice; 
    Depol_Ratio_Lat_Months = timetable_Depol_Ratio_monthly_avg.Total_Latitude_Ice; 
    Depol_Ratio_Lon_Months = timetable_Depol_Ratio_monthly_avg.Total_Longitude_Ice; 
    
    timetable_amsrmf_monthly_avg = retime(Total_timetable_amsrmf, 'monthly', @nanmean); 
    amsrmf_Monthly_avg = timetable_amsrmf_monthly_avg.Total_windamsrMF; 
    amsrmf_Time_Months = timetable_amsrmf_monthly_avg.Total_Profile_Time_New_Wind; 
    amsrmf_Lat_Months = timetable_amsrmf_monthly_avg.Total_Latitude_Wind;
    amsrmf_Lon_Months = timetable_amsrmf_monthly_avg.Total_Longitude_Wind; 
    
    save('CMOD_Monthly_avg_Vars.mat', ...
        'CMOD_Monthly_avg',...
        'CMOD_Time_Months',...
        'CMOD_Lat_Months',...
        'CMOD_Lon_Months',...
        '-v7.3') 
    
    save('Depol_Ratio_Monthly_avg_Vars.mat',...
        'Depol_Ratio_Monthly_avg',...
        'Depol_Ratio_Time_Months',...
        'Depol_Ratio_Lat_Months',...
        'Depol_Ratio_Lon_Months',...
        '-v7.3') 
    
    save('amsrmf_Monthly_avg_Vars.mat',...
        'amsrmf_Monthly_avg',...
        'amsrmf_Time_Months',...
        'amsrmf_Lat_Months',...
        'amsrmf_Lon_Months',...
        '-v7.3')
    
    %%
    
    
    
    
    
    
    
    
    