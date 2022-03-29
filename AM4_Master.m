%% Full Analysis runner
% Created by Marshall Borrus on June 7th, 2021
%
% The purpose of this document is to do EGR calculations and merge all the
% plotting intruments. The order of events is:
% 2. Variable Checking
    % 2.1 Mean, STD, and RMS Variable Creation
    % 2.2 Mean+STD and RMS ploting
    % 2.3 UTP Scatter Plots
% 3. EGR Calculation
    % 3.1 dry and moist values + N^2
    % 3.2 Mean values
% 4. EGR Plotting
    % 4.1 Time Series (1:100 days)
    % 4.2 Latitude Series 
    % 4.3 EGR Scatter Plots
% 5. Saturation Time
    % 5.1 Saturation Time vs N^2
    % 5.2 Saturation time vs EGR  
% 6. Stratification Plotting
    % 6.1 O'Gorman Plot
    % 6.2 Single lat comparison
%% 1.1 Data Loading - adding paths

clear
cd '/home/users/mborrus/Matlab_HPC'
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('/home/users/mborrus/Matlab_HPC/plots')
addpath('/home/users/mborrus/Matlab_HPC/scripts/cbrewer')
load('./data/EGR/lambda.mat')
load('AM4_Data.mat')

%% 1.2 Data Loading - Struct Creation

% Create the data paths for AM4 data
    OAKPATH = '/oak/stanford/schools/ees/aditis2/mborrus/';
    AM4_Data.path{1} = strcat(OAKPATH, 'AM4/Base/');
    AM4_Data.path{2} = strcat(OAKPATH, 'AM4/co2/');
    AM4_Data.path{3} = strcat(OAKPATH, 'AM4/quarter/');
    AM4_Data.path{4} = strcat(OAKPATH, 'AM4/m3/');
    AM4_Data.path{5} = strcat(OAKPATH, 'AM4/m2/');
    AM4_Data.path{6} = strcat(OAKPATH, 'AM4/m1/');
    AM4_Data.path{7} = strcat(OAKPATH, 'AM4/p1/');
    AM4_Data.path{8} = strcat(OAKPATH, 'AM4/p2/');
    AM4_Data.path{9} = strcat(OAKPATH, 'AM4/p3/');
    AM4_Data.path{10} = strcat(OAKPATH, 'AM4/SST_p4/');
    AM4_Data.path{11} = strcat(OAKPATH, 'AM4/SST_m4/');
    AM4_Data.path{12} = strcat(OAKPATH, 'AM4/SST_p8/');
    AM4_Data.path{13} = strcat(OAKPATH, 'AM4/SST_p12/');
    AM4_Data.path{14} = strcat(OAKPATH, 'AM4/SST_p16/');
    AM4_Data.path{15} = strcat(OAKPATH, 'AM4/SST_m8/');
    AM4_Data.path{16} = strcat(OAKPATH, 'AM4/SST_m12/');
    AM4_Data.path{17} = strcat(OAKPATH, 'AM4/SST_m16/');
    AM4_Data.path{18} = strcat('/scratch/users/mborrus/AM4/TOM/Q1/');
    AM4_Data.path{19} = strcat('/scratch/users/mborrus/AM4/TOM/Q2/');
    AM4_Data.path{20} = strcat('/scratch/users/mborrus/AM4/TOM/Q3/');
    AM4_Data.path{21} = strcat('/scratch/users/mborrus/AM4/TOM/Q4/');

% Create run names, # of files, title names, and plot colors    
    AM4_Data.run_type{1} = '+0C Base';      AM4_Data.File_Numbers{1} = [1:20]; AM4_Data.codename{1} = 'base'; AM4_Data.Color{1}= [0.6980    0.8745    0.5412; 0.2000    0.6275    0.1725]; AM4_Data.Line_Style{1}='-';
    AM4_Data.run_type{2} = '+4C Very Hot';  AM4_Data.File_Numbers{2} = [1:20]; AM4_Data.codename{2} = 'p4'; AM4_Data.Color{2}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{2}='-.';
    AM4_Data.run_type{3} = '-4C Very Cold'; AM4_Data.File_Numbers{3} = [1:20]; AM4_Data.codename{3} = 'm4'; AM4_Data.Color{3}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]; AM4_Data.Line_Style{3}='-.';
    AM4_Data.run_type{4} = '-3C Cold'; AM4_Data.File_Numbers{4} = [1:2]; AM4_Data.codename{4} = 'm3'; AM4_Data.Color{4}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.9; AM4_Data.Line_Style{4}=':';
    AM4_Data.run_type{5} = '-2C Colder'; AM4_Data.File_Numbers{5} = [1:2]; AM4_Data.codename{5} = 'm2'; AM4_Data.Color{5}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.8; AM4_Data.Line_Style{5}=':'; 
    AM4_Data.run_type{6} = '-1C Coldish'; AM4_Data.File_Numbers{6} = [1:2]; AM4_Data.codename{6} = 'm1'; AM4_Data.Color{6}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.7; AM4_Data.Line_Style{6}=':';
    AM4_Data.run_type{7} = '+1C Warmish'; AMta.File_Numbers{7} = [1:2]; AM4_Data.codename{7} = 'p1'; AM4_Data.Color{7}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.9; AM4_Data.Line_Style{7}=':';
    AM4_Data.run_type{8} = '+2C Warmer'; AM4_Data.File_Numbers{8} = [1:2]; AM4_Data.codename{8} = 'p2'; AM4_Data.Color{8}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.8; AM4_Data.Line_Style{8}=':';
    AM4_Data.run_type{9} = '+3C Warm'; AM4_Data.File_Numbers{9} = [1:2]; AM4_Data.codename{9} = 'p3'; AM4_Data.Color{9}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.7; AM4_Data.Line_Style{9}=':';
    AM4_Data.run_type{10} = '+4C QOB (hot)'; AM4_Data.File_Numbers{10} = [1:3]; AM4_Data.codename{10} = 'p4 QOB'; AM4_Data.Color{10}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{10}='-';
    AM4_Data.run_type{11} = '-4C QOB (cold)'; AM4_Data.File_Numbers{11} = [1:3]; AM4_Data.codename{11} = 'm4 QOB'; AM4_Data.Color{11}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]; AM4_Data.Line_Style{11}='-';
    AM4_Data.run_type{12} = '+8C QOB'; AM4_Data.File_Numbers{12} = [1:3]; AM4_Data.codename{12} = 'p8 QOB'; AM4_Data.Color{12}= [1.0000    1.0000    0.6980; 0.9922    0.5529    0.2353]; AM4_Data.Line_Style{12}='-';
    AM4_Data.run_type{13} = '+12C QOB'; AM4_Data.File_Numbers{13} = [1:3]; AM4_Data.codename{13} = 'p12 QOB'; AM4_Data.Color{13}= [0.9961    0.8510    0.4627; 0.9412    0.2314    0.1255]; AM4_Data.Line_Style{13}='-';
    AM4_Data.run_type{14} = '+16C QOB'; AM4_Data.File_Numbers{14} = [1:3]; AM4_Data.codename{14} = 'p16 QOB'; AM4_Data.Color{14}= [0.9961    0.6980    0.2980; 0.7412         0    0.1490]; AM4_Data.Line_Style{14}='-';
    AM4_Data.run_type{15} = '-8C QOB'; AM4_Data.File_Numbers{15} = [1:3]; AM4_Data.codename{15} = 'm8 QOB'; AM4_Data.Color{15}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.7; AM4_Data.Line_Style{15}='-';
    AM4_Data.run_type{16} = '-12C QOB'; AM4_Data.File_Numbers{16} = [1:3]; AM4_Data.codename{16} = 'm12 QOB'; AM4_Data.Color{16}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.8; AM4_Data.Line_Style{16}='-';
    AM4_Data.run_type{17} = '-16C QOB'; AM4_Data.File_Numbers{17} = [1:3]; AM4_Data.codename{17} = 'm16 QOB'; AM4_Data.Color{17}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.9; AM4_Data.Line_Style{17}='-';
    AM4_Data.run_type{18} = 'Q1'; AM4_Data.File_Numbers{18} = [1:4]; AM4_Data.codename{18} = 'Q1'; AM4_Data.Color{18}= [0.5961    0.3059    0.6392; 0.5961    0.3059    0.6392]*.9; AM4_Data.Line_Style{18}='-';
    AM4_Data.run_type{19} = 'Q2'; AM4_Data.File_Numbers{19} = [1:4]; AM4_Data.codename{19} = 'Q2'; AM4_Data.Color{19}= [0.5961    0.3059    0.6392; 0.5961    0.3059    0.6392]*.9; AM4_Data.Line_Style{19}='-';
    AM4_Data.run_type{20} = 'Q3'; AM4_Data.File_Numbers{20} = [1:4]; AM4_Data.codename{20} = 'Q3'; AM4_Data.Color{20}= [0.5961    0.3059    0.6392; 0.5961    0.3059    0.6392]*.9; AM4_Data.Line_Style{20}='-';
    AM4_Data.run_type{21} = 'Q4'; AM4_Data.File_Numbers{21} = [1:4]; AM4_Data.codename{21} = 'Q4'; AM4_Data.Color{21}= [0.5961    0.3059    0.6392; 0.5961    0.3059    0.6392]*.9; AM4_Data.Line_Style{21}='-';

% Anxillary data: Lat ranges, pressure ranges 
    %Lat Ranges
    AM4_Data.lat_range_names{1} = "North High"; AM4_Data.lat_range{1} = find(AM4_Data.lat >  50);
    AM4_Data.lat_range_names{2} = "South High"; AM4_Data.lat_range{2} = find(AM4_Data.lat <  -50);
    AM4_Data.lat_range_names{3} = "North Mid"; AM4_Data.lat_range{3} = find(AM4_Data.lat > 20 & AM4_Data.lat < 50);
    AM4_Data.lat_range_names{4} = "South Mid"; AM4_Data.lat_range{4} = find(AM4_Data.lat < -20 & AM4_Data.lat > -50);
    AM4_Data.lat_range_names{5} = "Equator"; AM4_Data.lat_range{5} = find(AM4_Data.lat > -20 & AM4_Data.lat < 20);
    AM4_Data.lat_range_names{6} = "45 deg N"; AM4_Data.lat_range{6} = find(AM4_Data.lat > 42.5 & AM4_Data.lat < 47.5);
    AM4_Data.lat_range_names{7} = "45 deg S"; AM4_Data.lat_range{7} = find(AM4_Data.lat < -42.5 & AM4_Data.lat > -47.5);
    AM4_Data.lat_range_names{8} = "40 deg N"; AM4_Data.lat_range{8} = find(AM4_Data.lat > 37.5 & AM4_Data.lat < 42.5);
    AM4_Data.lat_range_names{9} = "50 deg N"; AM4_Data.lat_range{9} = find(AM4_Data.lat > 47.5 & AM4_Data.lat < 52.5);
    AM4_Data.lat_range_names{10} = "20 deg N"; AM4_Data.lat_range{10} = find(AM4_Data.lat > 17.5 & AM4_Data.lat < 22.5);
    AM4_Data.lat_range_names{11} = "60-65 deg N"; AM4_Data.lat_range{11} = find(AM4_Data.lat > 57.5 & AM4_Data.lat < 67.5);
    AM4_Data.lat_range_names{12} = "1-10 deg N"; AM4_Data.lat_range{14} = find(AM4_Data.lat > 0 & AM4_Data.lat < 10);    
    AM4_Data.lat_range_names{13} = "20 deg S"; AM4_Data.lat_range{12} = find(AM4_Data.lat > -22.5 & AM4_Data.lat < -17.5);
    AM4_Data.lat_range_names{14} = "60-65 deg S"; AM4_Data.lat_range{13} = find(AM4_Data.lat < -57.5 & AM4_Data.lat > -62.5);
    AM4_Data.lat_range_names{15} = "1-10 deg S"; AM4_Data.lat_range{15} = find(AM4_Data.lat < 0 & AM4_Data.lat > -10);
    
    
    %Pressure Ranges
    AM4_Data.p_range_names{1} = "Lower Tropo"; AM4_Data.p_range{1} = find(AM4_Data.p<850 & AM4_Data.p>500); AM4_Data.p_color{1} = "#a1dab4";
    AM4_Data.p_range_names{2} = "Upper Tropo"; AM4_Data.p_range{2} = find(AM4_Data.p<500 & AM4_Data.p>100); AM4_Data.p_color{2} = "#41b6c4";
    AM4_Data.p_range_names{3} = "Stratosphere"; AM4_Data.p_range{3} = find(AM4_Data.p<100 & AM4_Data.p>1); AM4_Data.p_color{3} = "#225ea8";
    
    save(['AM4_Data.mat'],'AM4_Data')
    
%% 2. Variable Checking
%%%% 2.1 Mean, STD, and RMS Variable Creation
    %U - wind speeds
    for run_selection = [1];
        run_selection
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels> 1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,24,100);
                var_long = zeros(length(File_Numbers),144,90,24,100);
            temp_var = ncread(temp_path,'ucomp');
            var_1 = temp_var(:,:,1:24,1:100);
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'ucomp');
            var_diff(i-1,:,:,:,:) = temp_var(:,:,1:24,1:100) - var_1;
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100); 
            end 
            i
        end
        clear var_1 temp_var
        %%
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff); clear var_RMS var_Mean var_STD
        for lats = 1:15
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp);
%                 var_All(pres,lats,:,:) = squeeze(mean(var_long(:,:,range_temp,p_range{pres},:),[2,3,4]));
%                 var_All_abserr(pres,lats,:,:) = squeeze(mean(var_diff(:,:,range_temp,p_range{pres},:),[2,3,4]));
                var_All_RMS(pres,lats,:,:)   = squeeze(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2));                
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range{pres},:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range{pres},:),1,[1,2,3,4],'omitnan'));
            end
        end
        %%
%             AM4_Data.U_All{run_selection} = var_All;
%             AM4_Data.U_All_abserr{run_selection} = var_All_abserr;
            AM4_Data.U_All_RMS{run_selection} = var_All_RMS;
            AM4_Data.U_RMS{run_selection} = var_RMS;
            AM4_Data.U_Mean{run_selection} = var_Mean;
            AM4_Data.U_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD var_All_RMS var_All var_All_abserr
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% T - Temperatures
    for run_selection = [1 10:21];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,24,100);
                var_long = zeros(length(File_Numbers),144,90,24,100);
            temp_var = ncread(temp_path,'temp');
            var_1 = temp_var(:,:,1:24,1:100);
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'temp');
            var_diff(i-1,:,:,:,:) = temp_var(:,:,1:24,1:100) - var_1;
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff);
        for lats = 1:15
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_All_RMS(pres,lats,:,:) = squeeze(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2));                
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range{pres},:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range{pres},:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.T_All_RMS{run_selection} = var_All_RMS;
            AM4_Data.T_RMS{run_selection} = var_RMS;
            AM4_Data.T_Mean{run_selection} = var_Mean;
            AM4_Data.T_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD var_All_RMS
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% P - Precipitation
    for run_selection = [1 10:21];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,100);
                var_long = zeros(length(File_Numbers),144,90,100);
            temp_var = ncread(temp_path,'precip');
            var_1 = temp_var(:,:,1:100);
            var_long(i,:,:,:) = temp_var(:,:,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'precip');
            var_diff(i-1,:,:,:) = temp_var(:,:,1:100) - var_1;
            var_long(i,:,:,:) = temp_var(:,:,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,:), Nrun,rng*Nlon,Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.P_RMS{run_selection} = var_RMS;
            AM4_Data.P_Mean{run_selection} = var_Mean;
            AM4_Data.P_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% H - Specific Humidity
    for run_selection = [1:17];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTPS.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,100);
                var_long = zeros(length(File_Numbers),144,90,100);
            temp_var = ncread(temp_path,'sphum');
            var_1 = temp_var(:,:,1:100);
            var_long(i,:,:,:) = temp_var(:,:,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTPS.nc');
            temp_var = ncread(temp_path,'sphum');
            var_diff(i-1,:,:,:) = temp_var(:,:,1:100) - var_1;
            var_long(i,:,:,:) = temp_var(:,:,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,:), Nrun,rng*Nlon,Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.P_RMS{run_selection} = var_RMS;
            AM4_Data.P_Mean{run_selection} = var_Mean;
            AM4_Data.P_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% Single Value Means
    for run_selection = [1 10:21];
        for lats = 1:15
            for pres=1:3
                AM4_Data.U_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.U_Mean{run_selection}(pres,lats,:),[1,2,3]));
                AM4_Data.T_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.T_Mean{run_selection}(pres,lats,:),[1,2,3]));
                %AM4_Data.P_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.P_Mean{run_selection}(pres,lats,:),[1,2,3]));
            end
        end
    end
    save(['AM4_Data.mat'],'AM4_Data');
%%%% 2.2 Mean+STD and RMS ploting
%% 2.2 Mean + STD plotting
    %% U - Wind Speed Sat Time
    for run_selection = [1 10:21]; 
        %sat_times = zeros(3,6);
        for lat_band = [1:15];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath)
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
                i;
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 75]), title([AM4_Data.lat_range_names{lat_band},' RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_U_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
            run_selection 
            AM4_Data.U_Sat_Time{run_selection}(:,lat_band)=b_sat(:);
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% T - Temp Sat Time
    for run_selection = [1 10:21]; 
        for lat_band = [1:15];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath);
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.T_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 75]), title([AM4_Data.lat_range_names{lat_band},' T RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_T_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
        end
        run_selection
        AM4_Data.T_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% P - Precip Sat Time
    for run_selection = [15 16 17]; 
        sat_times = zeros(3,6);
        for lat_band = [1:11];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath);
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.P_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.01,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 inf]), title([AM4_Data.lat_range_names{lat_band},' P RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_P_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
        end
        run_selection
        AM4_Data.P_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')
%% 2.3 UTP Scatter Plots
    run_selection = [1:17];  
    for lat_band = [6]
        for p_range = 2;
        UTP_RMS_Sat_Time_Scatter = figure(1); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/Sat_Scatter/'); mkdir(savepath)
        clear x_data y_data
        C = ["#1b9e77";"#d95f02";"#7570b3"];
        for plot_choices = run_selection;
            x_data(1,plot_count,:)=mean(AM4_Data.T_Mean{plot_count}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(1,plot_count,:)=AM4_Data.U_Sat_Time{plot_choices}(:,lat_band);
            x_data(2,plot_count,:)=mean(AM4_Data.T_Mean{plot_count}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(2,plot_count,:)=AM4_Data.T_Sat_Time{plot_choices}(:,lat_band);
            x_data(3,plot_count,:)=mean(AM4_Data.T_Mean{plot_count}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(3,plot_count,:)=AM4_Data.P_Sat_Time{plot_choices}(:,lat_band);
            plot_count=plot_count+1;
        end
        x=squeeze(mean(x_data(:,:,:),[1,3]));
        for i = 1:3
            SPlot(i)=scatter(x(:),y_data(i,:,p_range),'MarkerEdgeColor',AM4_Data.p_color{p_range},'MarkerFaceColor',C(i));
            fitline(i)=plot(x(:),polyval(polyfit(x(:),y_data(i,:,p_range),1),x(:)),'Color',C(i));
            [r,P]=corrcoef(x(:),y_data(i,:,p_range)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),4))];
        end
        ylabel("days"); xlabel("Temp from base"); set(gca,'FontSize',12); axis([-inf inf 10 40]); 
        hleg = legend([SPlot(:);fitline(:)],["Wind";"Temp";"Precip";LegStr(1);LegStr(2);LegStr(3)]);
        title(["U T and P Saturation Times",AM4_Data.p_range_names{p_range}, AM4_Data.lat_range_names{lat_band}])
        saveas(UTP_RMS_Sat_Time_Scatter,strcat(savepath,'UTP_Sat_Time_vs_Temp_',AM4_Data.p_range_names{p_range},AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
        end
    end
%% 3. EGR Calculation
    
    % 3.1 dry and moist values + N^2
    for run_selection = [18:21];
        if run_selection>0
    clear lambda_term T u lat p theta dtheta_dz dtheta_dp du_z DRY_N2 DRY_egr MOIST_N2 MOIST_egr dtheta_dz_eff dtheta_dp_eff
    File_Numbers=AM4_Data.File_Numbers{run_selection};
    AM4_Data_Path = AM4_Data.path{run_selection};

%         if run_selection >= 10;
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
%         else
%             temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTco2.nc');
%         end

    pressure_levels = ncread(temp_path,'pfull');
    lat = ncread(temp_path,'grid_yt');
    days = (1:100);
    Pressure_range = find(pressure_levels<850 & pressure_levels>200);
    p = pressure_levels(Pressure_range); 

    H = 7300; % scale height, m
    z = -H*log(p/1000);
    omega = 2*pi/(24*3600);
    f = 2*omega*sin(2*pi*lat/360);
    g = 9.8;

    lambda_input = interp1(lambda_base(:,1),lambda_base(:,2),lat);

    for run_N = 1:length(File_Numbers)
    %for run_N = 1

%         if run_selection >= 10;
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(run_N)),'/dailyUTP.nc');
%         else
%             temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(run_N)),'/dailyUTco2.nc');
%         end

        %Load U and T
        start_loc = [ 1 1 Pressure_range(1) 1];
        count = [inf inf length(Pressure_range) 100];
        u = ncread(temp_path,'ucomp',start_loc,count);
        T = ncread(temp_path,'temp',start_loc,count);

        %Zonal mean
        T = squeeze(mean(T(:,:,:,:),1)); %average across longitude
        u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
        [la pr d] = size(u);

        if exist('theta') == 0;
            theta = zeros(length(File_Numbers),la, pr, d);
            dtheta_dz = zeros(length(File_Numbers),la, pr, d);
            dtheta_dp = zeros(length(File_Numbers),la, pr, d);
            du_z = zeros(length(File_Numbers),la, pr, d);
            DRY_N2 = zeros(length(File_Numbers),la,pr,d);
            DRY_egr = zeros(length(File_Numbers),la,pr,d);
            MOIST_N2 = zeros(length(File_Numbers),la,pr,d);
            MOIST_egr = zeros(length(File_Numbers),la,pr,d);
            dtheta_dz_eff = zeros(length(File_Numbers),la,pr,d);
            dtheta_dp_eff = zeros(length(File_Numbers),la,pr,d);
            lambda_term = zeros(length(File_Numbers),la,pr,d);

                "vars created"
        else
                "vars already exist"
        end

        %DRY RUNS 
        lambda = 0;
        for lat_N = 1:la
            for day_N = 1:d
                [dtheta_dz(run_N,lat_N,:,day_N),dtheta_dp(run_N,lat_N,:,day_N),lambda_term(run_N,lat_N,:,day_N)] = eff_stat_stab(p', T(lat_N,:,day_N), lambda);
            end
        end

        %MOIST RUNS

        for lat_N = 1:la
        lambda = lambda_input(lat_N);
            for day_N = 1:d
                [dtheta_dz_eff(run_N,lat_N,:,day_N),dtheta_dp_eff(run_N,lat_N,:,day_N),lambda_term(run_N,lat_N,:,day_N)] = eff_stat_stab(p', T(lat_N,:,day_N), lambda);
            end
        end

        pressure_term = (1000./p).^(2/7);

        for j = 1:pr
            theta(run_N,:,j,:) = T(:,j,:).*pressure_term(j);       
        end
                "theta ran"

        for i=2:pr-1
        du_z(run_N,:,i,:) = (u(:,i+1,:)-u(:,i-1,:))/...
            (z(i+1)-z(i-1));       
        end
                "derivatives ran"

        du_z(run_N,:,1,:) = (u(:,2,:)-u(:,1,:))/(z(2)-z(1));
        du_z(run_N,:,pr,:) = (u(:,pr,:)-u(:,pr-1,:))/(z(pr)-z(pr-1));


        for i=1:pr
            DRY_N2(run_N,:,i,:) = (g./theta(run_N,:,i,:)) .* dtheta_dz(run_N,:,i,:);
        end

        for i = 1:la
            DRY_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(run_N,i,:,:))./sqrt(squeeze(DRY_N2(run_N,i,:,:))));
        end

        for i=1:pr
            MOIST_N2(run_N,:,i,:) = (g./theta(run_N,:,i,:)) .* dtheta_dz_eff(run_N,:,i,:);
        end

        for i = 1:la
            MOIST_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(run_N,i,:,:))./sqrt(squeeze(MOIST_N2(run_N,i,:,:))));
        end
        run_N
    end

    %save([strcat('./data/EGR/AM4/AM4_',AM4_Data.codename{run_selection},'.mat')],'lambda_term', 'T','u', 'lat', 'p', 'theta','dtheta_dz','dtheta_dp','du_z','DRY_N2', 'DRY_egr','MOIST_N2','MOIST_egr','dtheta_dz_eff','dtheta_dp_eff');
    AM4_Data.dry_egr{run_selection} = DRY_egr;
    AM4_Data.moist_egr{run_selection} = MOIST_egr;
    AM4_Data.dry_N2{run_selection} = DRY_N2;
    AM4_Data.moist_N2{run_selection} = MOIST_N2;
    "saved"
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')
    
    % 3.2 Mean values
    for run_selection = [18:21];
        if run_selection > 0;
            for i = 1:6
                dry_means(i)=nanmean(AM4_Data.dry_egr{run_selection}(:,AM4_Data.lat_range{i},:,:),[1,2,3,4]);
                moist_means(i)=nanmean(AM4_Data.moist_egr{run_selection}(:,AM4_Data.lat_range{i},:,:),[1,2,3,4]);
            end
            AM4_Data.dry_egr_mean{run_selection} = dry_means;
            AM4_Data.moist_egr_mean{run_selection} = moist_means;
            clear means moist_means dry_means
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')
%% 4. EGR Plotting
% Select which runs you want plotted    
plot_choices = [1 15 16 17];
  %% 4.1 Time Series (1:100 days)
    %Time Series
    for lat_band = 6  
        EGR_TimeSeries = figure(1); clf; hold on; plot_count = 1; legend_names = [];
        for plot_num = plot_choices
            temp_plots(plot_count) = plot(1:100, squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
            temp_plots(plot_count) = plot(1:100, squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(1,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' moist']; plot_count = plot_count +1;
        end
        ylabel("1/day"); xlabel("days"); legend(legend_names,'FontSize',12); legend('boxoff'); set(gca,'FontSize',12)
        title(["EGR Time Series ", AM4_Data.lat_range_names{lat_band}])
        saveas(EGR_TimeSeries,[strcat('./plots/Master/EGR/Time_Series/Time_Series',AM4_Data.lat_range_names{lat_band},'.png')])
    end
    %% Time Series - difference
    for lat_band = 6  
        EGR_TimeSeries_diff = figure(2); clf;subplot(2,1,1); hold on; plot_count = 1; legend_names = [];

        for plot_num = plot_choices
            temp_plots(plot_count) = plot(1:100, (squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) - squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) ) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
        end
        ylabel("1/day"); xlabel("days"); legend(legend_names,'FontSize',12); legend('boxoff'); axis([0 100 -inf inf]); set(gca,'FontSize',12)
        title(["EGR Time Series - Difference ", AM4_Data.lat_range_names{lat_band}])
        
        subplot(2,1,2); hold on; plot_count = 1; legend_names = [];
        for plot_num = plot_choices
            temp_plots(plot_count) = plot(1:100, (squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) - squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) ) ./ squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 3] ) ) *100  ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
        end
        ylabel("% of dry EGR"); xlabel("days"); legend(legend_names,'FontSize',12); legend('boxoff'); axis([0 100 10 30]);set(gca,'FontSize',12)
        title(["EGR Time Series - Percent Difference ", AM4_Data.lat_range_names{lat_band}])
        saveas(EGR_TimeSeries_diff,[strcat('./plots/Master/EGR/Time_Series/Time_Series_Diff_',AM4_Data.lat_range_names{lat_band},'.png')])
    end
   %% 4.2 Latitude Series 
   %Latitude Series
    for lat_band = 6  
        EGR_LatSeries = figure(3); clf; hold on; plot_count = 1; legend_names = [];
        for plot_num = plot_choices
            temp_plots(plot_count) = plot(AM4_Data.lat, squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,:,:,:), [1 3 4] ) ) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
            temp_plots(plot_count) = plot(AM4_Data.lat, squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,:,:,:), [1 3 4] ) ) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(1,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' moist']; plot_count = plot_count +1;
        end
        ylabel("1/day"); xlabel("latitutde"); legend(legend_names,'FontSize',12); legend('boxoff'); set(gca,'FontSize',12)
        title(["EGR Lat Series - days 1:100"])
        saveas(EGR_LatSeries,[strcat('./plots/Master/EGR/Lat_Series/Lat_Series',AM4_Data.lat_range_names{lat_band},'.png')])
    end
   %Latitude Series - difference
    for lat_band = 6  
        EGR_LatSeries_diff = figure(4); clf;subplot(2,1,1); hold on; plot_count = 1; legend_names = [];

        for plot_num = plot_choices
            temp_plots(plot_count) = plot(AM4_Data.lat, (squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,:,:,:), [1 3 4] ) ) - squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,:,:,:), [1 3 4] ) )) *86400 ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
        end
        ylabel("1/day"); xlabel("latitutde"); legend(legend_names,'FontSize',12); legend('boxoff'); axis([-90 90 -inf inf]); set(gca,'FontSize',12)
        title(["EGR Lat Series - Difference ", AM4_Data.lat_range_names{lat_band}])
        
        subplot(2,1,2); hold on; plot_count = 1; legend_names = [];
        for plot_num = plot_choices
            temp_plots(plot_count) = plot(AM4_Data.lat, (squeeze( nanmean( AM4_Data.moist_egr{plot_num}(:,:,:,:), [1 3 4] ) ) - squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,:,:,:), [1 3 4] ) )) ./ squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,:,:,:), [1 3 4] ) )*100  ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
        end
        ylabel("% of dry EGR"); xlabel("latitutde"); legend(legend_names,'FontSize',12); legend('boxoff'); axis([-90 90 0 70]);set(gca,'FontSize',12)
        title(["EGR Lat Series - Percent Difference ", AM4_Data.lat_range_names{lat_band}])
        saveas(EGR_LatSeries_diff,[strcat('./plots/Master/EGR/Lat_Series/Lat_Series_Diff_',AM4_Data.lat_range_names{lat_band},'.png')])
    end
 
   %% 4.3 EGR Scatter Plot 
    run_selection = [1:17];  
    % EGR vs T (moist and dry)
    for lat_band = [6]
        EGR_Scatter = figure(1); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        clear x_data y_data
        C = ["#7fc97f";"#756bb1";"#756bb1"];
        for plot_choices = run_selection;
            x_data(1,plot_count)=mean(AM4_Data.T_Mean{plot_count}(:,lat_band,:)-AM4_Data.T_Mean{1}(:,lat_band,:),[1,3]);y_data(1,plot_count)=AM4_Data.dry_egr_mean{plot_choices}(lat_band)*86400;
            x_data(2,plot_count)=mean(AM4_Data.T_Mean{plot_count}(:,lat_band,:)-AM4_Data.T_Mean{1}(:,lat_band,:),[1,3]);y_data(2,plot_count)=AM4_Data.moist_egr_mean{plot_choices}(lat_band)*86400;
            plot_count=plot_count+1;
        end
        x=squeeze(mean(x_data(:,:,:),[1,3]));
        for i = 1:2
            SPlot(i)=scatter(x(:),y_data(i,:),'MarkerEdgeColor',C(i),'MarkerFaceColor',C(i));
            fitline(i)=plot(x(:),polyval(polyfit(x(:),y_data(i,:),1),x(:)),'Color',C(i));
            [r,P]=corrcoef(x(:),y_data(i,:)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),6))];
        end
        ylabel("EGR 1/Days"); xlabel("Temp from base"); set(gca,'FontSize',12); axis([-inf inf 1 3]); 
        hleg = legend([SPlot(:);fitline(:)],["Dry EGR";"Moist EGR";LegStr(1);LegStr(2)],'Location','NorthWest');
        title(["Dry and Moist EGR", AM4_Data.lat_range_names{lat_band}])
        set(gca,'fontname','times new roman'); set(gcf,'color','w');
        saveas(EGR_Scatter,strcat(savepath,'EGR_Scatter_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
    end
    %EGR difference vs T
    for lat_band = [6]
        EGR_Scatter = figure(1); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        clear x_data y_data
        C = ["#7fc97f";"#756bb1";"#756bb1"];
        for plot_choices = run_selection;
            x_data(1,plot_count)=mean(AM4_Data.T_Mean{plot_count}(:,lat_band,:)-AM4_Data.T_Mean{1}(:,lat_band,:),[1,3]);y_data(1,plot_count)=AM4_Data.dry_egr_mean{plot_choices}(lat_band)*86400;
            x_data(2,plot_count)=mean(AM4_Data.T_Mean{plot_count}(:,lat_band,:)-AM4_Data.T_Mean{1}(:,lat_band,:),[1,3]);y_data(2,plot_count)=AM4_Data.moist_egr_mean{plot_choices}(lat_band)*86400; 
            x_data(3,plot_count)=mean(AM4_Data.T_Mean{plot_count}(:,lat_band,:)-AM4_Data.T_Mean{1}(:,lat_band,:),[1,3]);y_data(3,plot_count)=(y_data(2,plot_count)-y_data(1,plot_count))./y_data(2,plot_count);
            plot_count=plot_count+1;
        end
        x=squeeze(mean(x_data(:,:,:),[1,3]));
        for i = 3
            SPlot(i)=scatter(x(:),y_data(i,:),'MarkerEdgeColor',C(i),'MarkerFaceColor',C(i));
            fitline(i)=plot(x(:),polyval(polyfit(x(:),y_data(i,:),1),x(:)),'Color',C(i));
            [r,P]=corrcoef(x(:),y_data(i,:)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),6))];
        end
        ylabel("% Difference"); xlabel("Temp from base"); set(gca,'FontSize',12); axis([-inf inf .1 .2]); 
        hleg = legend([SPlot(3);fitline(3)],["EGR Difference";LegStr(3)],'Location','NorthWest');
        title(["Dry and Moist EGR", AM4_Data.lat_range_names{lat_band}])
        set(gca,'fontname','times new roman'); set(gcf,'color','w');
        saveas(EGR_Scatter,strcat(savepath,'EGR_Scatter_Difference_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
    end
    
%% 5. Saturation Times
    %% 5.1 Saturation Time vs N^2
    run_selection = [1:17];
    for lat_band = [6]
        N2_Sat_Time_Scatter = figure(1); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/Sat_Scatter/'); mkdir(savepath)
        clear x_data y_data
        for plot_choices = run_selection;
            for p_range = 1:2;
            x_data(p_range,plot_count)=mean(AM4_Data.dry_N2{plot_choices}(:,lat_band,AM4_Data.p_range{p_range},:),[1,2,3,4]);
            y_data(p_range,plot_count)=AM4_Data.U_Sat_Time{plot_count}(p_range,lat_band);
            end
            plot_count=plot_count+1;
        end
        x=x_data;
        for i = 1:2
            SPlot(i)=scatter(x(i,:),y_data(i,:),'MarkerEdgeColor',AM4_Data.p_color{i},'MarkerFaceColor',AM4_Data.p_color{i});
            %fitline(i)=plot(x(i,:),polyval(polyfit(x(i,:),y_data(i,:),1),x(i,:)),'Color',AM4_Data.p_color{i});
            [r,P]=corrcoef(x(i,:),y_data(i,:)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),6))];
        end
        [a,I]=sort([x(1,:) x(2,:)]); b = [y_data(1,:) y_data(2,:)]; b=b(I);     
        ylabel("Saturation Time"); xlabel("N^2 Time Average"); set(gca,'FontSize',12); 
        hleg = legend([SPlot(:);fitline(:)],[AM4_Data.p_range_names{1};AM4_Data.p_range_names{2};LegStr(1);LegStr(2);],'Location','SouthOutside','Numcolumns',2);
        title(["N^2 vs Saturation Time", AM4_Data.lat_range_names{lat_band}])
        saveas(N2_Sat_Time_Scatter,strcat(savepath,'N2_Sat_Time_Scatter_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
    end
    %% 5.2 Saturation time vs EGR  
    run_selection = [1:17];
    day_I = [1 3 5];
    %% 20 day grouping - dry
    for lat_band = [6]
        for p_range = 1:2
        clear x_data y_data; days = [1:20;21:40;41:60;61:80;81:100]; plot_count = 1;
        for plot_choices = run_selection;
            for d = day_I;
            x_data(p_range,plot_count,d)=mean(AM4_Data.dry_egr{plot_choices}(:,AM4_Data.lat_range{lat_band},:,days(d,:)),[1,2,3,4]);
            y_data(p_range,plot_count,d)=AM4_Data.U_Sat_Time{plot_count}(p_range,lat_band);
            end
            plot_count=plot_count+1;
        end
        EGR_Sat_Time_Scatter = figure(p_range); clf; hold on; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        x=x_data*86400;
        C = ["#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6"];
        for i = p_range
            for d = day_I;
            SPlot(i,d)=scatter(x(i,:,d),y_data(i,:,d),'filled','MarkerFaceColor',C(d));
            fitline(i,d)=plot(x(i,:,d),polyval(polyfit(x(i,:,d),y_data(i,:,d),1),x(i,:,d)),'Color',C(d));
            [r,P]=corrcoef(x(i,:,d),y_data(i,:,d)); r2(i,d)=r(2,1)^2;p(i,d)=P(2,1); LegStr(i,d)=["r^2 = " + num2str(round(r2(i,d),2)) + "; p = " + num2str(round(p(i,d),6))];
            end
        end
        ylabel("Saturation Time"); xlabel("EGR Time Average (1/day)"); set(gca,'FontSize',12); 
        hleg = legend([SPlot(p_range,day_I),fitline(p_range,day_I)],["days "+days(day_I,1)+" to "+days(day_I,end);LegStr(p_range,day_I)']','Location','NorthEast','Numcolumns',2);
        title(["Dry EGR vs Saturation Time", AM4_Data.lat_range_names{lat_band}, AM4_Data.p_range_names{p_range}]); axis([1 2.5 14 26])
        saveas(EGR_Sat_Time_Scatter,strcat(savepath,'EGR_Sat_Time_dry_Scatter_',AM4_Data.p_range_names{p_range},'_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
        end
    end
    %% 20 day grouping - moist
    for lat_band = [6]
        for p_range = 1:2
        clear x_data y_data ; days = [1:20;21:40;41:60;61:80;81:100]; plot_count = 1;
        for plot_choices = run_selection;
            for d = day_I;
            x_data(p_range,plot_count,d)=mean(AM4_Data.moist_egr{plot_choices}(:,AM4_Data.lat_range{lat_band},:,days(d,:)),[1,2,3,4]);
            y_data(p_range,plot_count,d)=AM4_Data.U_Sat_Time{plot_count}(p_range,lat_band);
            end
            plot_count=plot_count+1;
        end
        EGR_Sat_Time_Scatter = figure(p_range); clf; hold on; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        x=x_data*86400;
        C = ["#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6"];
        for i = p_range
            for d = day_I;
            SPlot(i,d)=scatter(x(i,:,d),y_data(i,:,d),'filled','MarkerFaceColor',C(d));
            fitline(i,d)=plot(x(i,:,d),polyval(polyfit(x(i,:,d),y_data(i,:,d),1),x(i,:,d)),'Color',C(d));
            [r,P]=corrcoef(x(i,:,d),y_data(i,:,d)); r2(i,d)=r(2,1)^2;p(i,d)=P(2,1); LegStr(i,d)=["r^2 = " + num2str(round(r2(i,d),2)) + "; p = " + num2str(round(p(i,d),6))];
            end
        end
        ylabel("Saturation Time"); xlabel("EGR Time Average (1/day)"); set(gca,'FontSize',12); 
        hleg = legend([SPlot(p_range,day_I),fitline(p_range,day_I)],["days "+days(day_I,1)+" to "+days(day_I,end);LegStr(p_range,day_I)']','Location','NorthEast','Numcolumns',2);
        title(["Moist EGR vs Saturation Time", AM4_Data.lat_range_names{lat_band}, AM4_Data.p_range_names{p_range}]); axis([1.2 3 14 26])
        saveas(EGR_Sat_Time_Scatter,strcat(savepath,'EGR_Sat_Time_Moist_Scatter_',AM4_Data.p_range_names{p_range},'_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
        end
    end
    %% Dry EGR plotting - no grouping
    for lat_band = [6]
        for p_range = 1:3;
        EGR_Sat_Time_Scatter = figure(p_range); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        clear x_data y_data
        for plot_choices = run_selection;         
            x_data(p_range,plot_count)=AM4_Data.dry_egr_mean{plot_choices}(lat_band);
            y_data(p_range,plot_count)=AM4_Data.U_Sat_Time{plot_count}(p_range,lat_band);
            plot_count=plot_count+1;
        end
        x=x_data*86400;
        for i = p_range
            SPlot(i)=scatter(x(i,:),y_data(i,:),'MarkerEdgeColor',AM4_Data.p_color{p_range},'MarkerFaceColor',AM4_Data.p_color{p_range});
            fitline(i)=plot(x(i,:),polyval(polyfit(x(i,:),y_data(i,:),1),x(i,:)),'Color',AM4_Data.p_color{p_range});
            [r,P]=corrcoef(x(i,:),y_data(i,:)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),6))];
        end
        ylabel("Saturation Time"); xlabel("EGR Time Average (1/day)"); set(gca,'FontSize',12); 
        hleg = legend([SPlot(p_range);fitline(p_range)],[AM4_Data.p_range_names{p_range};LegStr(p_range)],'Location','SouthOutside','Numcolumns',2);
        title(["Dry EGR vs Saturation Time", AM4_Data.lat_range_names{lat_band}, AM4_Data.p_range_names{p_range}]); axis([1 3 10 40])
        saveas(EGR_Sat_Time_Scatter,strcat(savepath,'EGR_dry_Sat_Time_Scatter_',AM4_Data.p_range_names{p_range},'_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
        end
    end
    %% Moist EGR plotting - no grouping
    for lat_band = [6]
        for p_range = 1:3;
        EGR_Sat_Time_Scatter = figure(p_range); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('./plots/Master/EGR_Scatter/'); mkdir(savepath)
        clear x_data y_data
        for plot_choices = run_selection;         
            x_data(p_range,plot_count)=AM4_Data.moist_egr_mean{plot_choices}(lat_band);
            y_data(p_range,plot_count)=AM4_Data.U_Sat_Time{plot_count}(p_range,lat_band);
            plot_count=plot_count+1;
        end
        x=x_data*86400;
        for i = p_range
            SPlot(i)=scatter(x(i,:),y_data(i,:),'MarkerEdgeColor',AM4_Data.p_color{p_range},'MarkerFaceColor',AM4_Data.p_color{p_range});
            fitline(i)=plot(x(i,:),polyval(polyfit(x(i,:),y_data(i,:),1),x(i,:)),'Color',AM4_Data.p_color{p_range});
            [r,P]=corrcoef(x(i,:),y_data(i,:)); r2(i)=r(2,1)^2;p(i)=P(2,1); LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),6))];
        end
        ylabel("Saturation Time"); xlabel("EGR Time Average (1/day)"); set(gca,'FontSize',12); 
        hleg = legend([SPlot(p_range);fitline(p_range)],[AM4_Data.p_range_names{p_range};LegStr(p_range)],'Location','SouthOutside','Numcolumns',2);
        title(["Moist EGR vs Saturation Time", AM4_Data.lat_range_names{lat_band}, AM4_Data.p_range_names{p_range}]); axis([1 3 10 40])
        saveas(EGR_Sat_Time_Scatter,strcat(savepath,'EGR_moist_Sat_Time_Scatter_',AM4_Data.p_range_names{p_range},'_',AM4_Data.lat_range_names{lat_band},'.png'))     
        clear x x_data y_data LegStr C P SPlot UTP_RMS_Sat_Time_Scatter fitline hleg p_range r2 r p plot_count plot_choices legend_names savepath i
        end
    end
%% 6. Stratification Plotting
   %% 6.1 O'Gorman Plot
    for plot_num = 4:14
    Ogorman_Tri_Plot = figure(5); clf; subplot(3,1,1); hold on; plot_count = 1; legend_names = [];
    dryN2_temp = squeeze(nanmean(AM4_Data.dry_N2{plot_num}(:,:,:,:),[1,4]))'.*10^5;
    moistN2_temp = squeeze(nanmean(AM4_Data.moist_N2{plot_num}(:,:,1:12,:),[1,4]))'.*10^5;
    limits = [-80 80 50 925]; lines = [5 10 15 20 25 30]; minor_lines = [.2 .3 .4 .5 .6 .7 .8 .9];
    
    subplot(3,1,1)
    [M,c] = contour(AM4_Data.lat,AM4_Data.p, dryN2_temp,lines,'ShowText','on');
    ylabel('pressure'), title(strcat("N^2 dry * 10^5 - ", AM4_Data.codename{plot_num})), c.LineWidth = 2;
    axis(limits), set(gca, 'YDir','reverse')

    subplot(3,1,2)
    [M,c] = contour(AM4_Data.lat,AM4_Data.p, moistN2_temp,lines,'ShowText','on');
    ylabel('pressure'), title(strcat("N^2 dry * 10^5 - ", AM4_Data.codename{plot_num})), c.LineWidth = 2;
    axis(limits), set(gca, 'YDir','reverse')

    subplot(3,1,3)
    [M,c] = contour(AM4_Data.lat,AM4_Data.p, moistN2_temp./dryN2_temp,minor_lines,'ShowText','on');
    xlabel('latitude'), ylabel('pressure'), title('moist/dry'), c.LineWidth = 2; 
    axis(limits); set(gca, 'YDir','reverse')
    
    saveas(Ogorman_Tri_Plot,[strcat('./plots/Master/Strat/Ogorman_Fig_4/Triplot_',AM4_Data.codename{plot_num},'.png')])
    end
   
   %% 6.2 Single lat comparison
   for lat_band = 6
    N2_Comparison = figure(110); clf; subplot(3,1,1), hold on; plot_count = 1; legend_names = [];
    for plot_num = plot_choices
        temp_plots(plot_count) = plot(squeeze( nanmean( AM4_Data.dry_N2{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 4] ) ) .*10^5, AM4_Data.p ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' dry']; plot_count = plot_count +1;
    end
        ylabel('pressure'); legend(legend_names,'FontSize',12); legend('boxoff'); set(gca,'FontSize',12)
        title(["Dry N^2  ", AM4_Data.lat_range_names{lat_band}]), set(gca, 'YDir','reverse'), axis([0 45 200 900])
    subplot(3,1,2); hold on; plot_count = 1; legend_names = [];
    for plot_num = plot_choices
        temp_plots(plot_count) = plot(squeeze( nanmean( AM4_Data.moist_N2{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 4] ) ) .*10^5, AM4_Data.p ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num} ' moist']; plot_count = plot_count +1;
    end 
        title("Moist N^2"); ylabel('pressure'); legend(legend_names,'FontSize',12); legend('boxoff'); set(gca,'FontSize',12), set(gca, 'YDir','reverse'), axis([0 45 200 900])
    subplot(3,1,3); hold on; plot_count = 1; legend_names = [];
    for plot_num = plot_choices
        temp_plots(plot_count) = plot((squeeze( nanmean( AM4_Data.moist_N2{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 4] ) ) .*10^5)./(squeeze( nanmean( AM4_Data.dry_N2{plot_num}(:,AM4_Data.lat_range{lat_band},:,:), [1 2 4] ) ) .*10^5), AM4_Data.p ,'LineWidth',2, 'Color',AM4_Data.Color{plot_num}(2,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num}]; plot_count = plot_count +1;
    end 
        ylabel('pressure'); xlabel("10^5 * N^2"); legend(legend_names,'FontSize',12); legend('boxoff'); set(gca,'FontSize',12),set(gca, 'YDir','reverse'), axis([-inf 1 200 900]), title("Moist/Dry N^2")
    saveas(Ogorman_Tri_Plot,[strcat('./plots/Master/Strat/N2_compare/compare_',AM4_Data.codename{plot_num},'.png')])
   end
    %% 7. EGR, N, and du/dz mean reporting
clear dryegrmean moistegrmean dryNmean moistNmean dryegrrecalc dudzmean
lat = AM4_Data.lat;
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat(AM4_Data.lat_range{6})/360);

 i=1;for run_selection = [17 16 15 11 1 10 12 13 14];
        dryegrmean(:,i) = mean(AM4_Data.dry_egr{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 3 4])*86400;
        moistegrmean(:,i) = mean(AM4_Data.moist_egr{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 3 4])*86400;
        dryNmean(:,i) = mean(AM4_Data.dry_N2{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 3 4]).^(1/2)*86400;
        moistNmean(:,i) = mean(AM4_Data.moist_N2{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 3 4]).^(1/2)*86400;
        dudzmean(:,i) = mean(AM4_Data.du_dz{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 3 4])*86400;
        for latnum = 1:3
        dryegrrecalc(latnum,i) = mean(abs(f(latnum).*(AM4_Data.du_dz{run_selection}(:,AM4_Data.lat_range{6}(latnum),:,1:20))./(AM4_Data.dry_N2{run_selection}(:,AM4_Data.lat_range{6}(latnum),:,1:20).^(1/2))),[1 2 3 4])'*86400;
        end
        i = i+1
 end

mean(dryegrmean)'
mean(dryNmean)'
mean(dudzmean)'

 
backcalcegr = f.*(dudzmean./dryNmean);
mean(backcalcegr,1)'*86400 - mean(dryegrmean)'
mean(dryegrrecalc,1)' - mean(dryegrmean)'