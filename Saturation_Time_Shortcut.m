%% SatTimeSimple
  for run_selection = [1 10:21]; 
        for lat_band = [1:15];
            range_temp = 1:90; rng = length(range_temp);
            clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                i;
            end 
            sat_times(:,lat_band)=b_sat(:) ;
            run_selection 
            AM4_Data.U_Sat_Time{run_selection}(:,lat_band)=b_sat(:);
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% T - Temp Sat Time
    for run_selection = [1 10:21]; 
        for lat_band = [1:15];
            range_temp = 1:90; rng = length(range_temp);
            clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.T_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
            end
            sat_times(:,lat_band)=b_sat(:) ;
        end
        run_selection
        AM4_Data.T_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')

%%
% Scratch Plots
run_selection = [17 16 15 11 1 10 12 13 14];
[sorted, index] = sort(run_selection);
lat_selection = [12 15 10 13 6 7 11 14];

Lower_Tropo = []
Upper_Tropo = []
Strat       = []

j = 1
for i = run_selection;
    Lower_Tropo(j,:) = AM4_Data.U_Sat_Time{i}(1,lat_selection);
    Upper_Tropo(j,:) = AM4_Data.U_Sat_Time{i}(2,lat_selection);
    Strat(j,:)       = AM4_Data.U_Sat_Time{i}(3,lat_selection);
    j=j+1;
end
%%

% for run_selection = [1 10:21];
%     run_selection
%     File_Numbers=AM4_Data.File_Numbers{run_selection};
%     AM4_Data_Path = AM4_Data.path{run_selection};
%     temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
%     
%     % Get the pressure and lats values
%     pressure_levels = ncread(temp_path,'pfull');
%     p_range{1} = find(pressure_levels<850 & pressure_levels>500);
%     p_range{2} = find(pressure_levels<500 & pressure_levels>100);
%     p_range{3} = find(pressure_levels<100 & pressure_levels>1);
%     lat_ranges = ncread(temp_path,'grid_yt');
%     for i = 1:length(File_Numbers)
%         if i == 1
%             var_diff = zeros(length(File_Numbers)-1,144,90,24,100);
%             var_long = zeros(length(File_Numbers),144,90,24,100);
%             temp_var = ncread(temp_path,'ucomp');
%             var_1 = temp_var(:,:,1:24,1:100);
%             var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
%         else
%             temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
%             temp_var = ncread(temp_path,'ucomp');
%             var_diff(i-1,:,:,:,:) = temp_var(:,:,1:24,1:100) - var_1;
%             var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
%         end   
%     end
% end
%     save(['AM4_Data.mat'],'AM4_Data');


clear
cd '/home/users/mborrus/Matlab_HPC'
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('/home/users/mborrus/Matlab_HPC/plots')
addpath('/home/users/mborrus/Matlab_HPC/scripts/cbrewer')
load('./data/EGR/lambda.mat')
load('AM4_Data.mat')
Tmaster=zeros([21, 5, 144,    90,    33-9,   100]);
for run_selection = [19];
    File_Numbers=AM4_Data.File_Numbers{run_selection};
    if length(File_Numbers)>5
        q = File_Numbers(1:5);
    else
        q = File_Numbers;
    end
    q
    for i=q
        AM4_Data_Path = AM4_Data.path{run_selection};
        startLoc = [1 1 1 1]; % Start location along each coordinate
        count  = [Inf Inf 24 100]; % Read until the end of each dimension
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
        Tmaster(run_selection,i,:,:,:,:)=ncread(temp_path,'temp',startLoc,count);
        i
    end
    run_selection
end
save(['Tcomplete.mat'],'Tmaster', '-v7.3');        