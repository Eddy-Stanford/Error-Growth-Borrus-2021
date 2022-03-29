%% RMS_plots V2
% Script to be run in sherlock to make plots
% Marshall Borrus
% This is where the data is held

Dycore_Data_Path = '/scratch/users/mborrus/h0/';
cd /home/users/mborrus/Matlab_HPC
mkdir plots/RMS     %Where these specific plots should go
%load('AM4_Data.mat')
load('dy_data.mat')
%%

% Get the pressure and lats values
% size(P) = 40; size(lat) = 64; size(lon) = 128;

load('axis_stuff_64.mat','lat','lon','P');
days = (1:400);

% Set pressure and lat ranges
whole_P_range = find(P<850 & P>1);
pressure_levels = P(whole_P_range);  clear P;
Nlon = length(lon);
Nlat = length(lat);
Np   = length(pressure_levels);
Ntime= length(days);
Nruns = 21; Nrun = Nruns-1;  %%%%%%%%% CHANGE THIS IF YOU NEED TO %%%%%%%%%%%
lat_ranges = lat;     clear lat;

dy_data.path{1} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/h0/ary'; dy_data.codename{1} = "h0";
dy_data.path{2} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/10/'; dy_data.codename{2} = "10";
dy_data.path{3} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/20/'; dy_data.codename{3} = "20";
dy_data.path{4} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/30/'; dy_data.codename{4} = "30";
dy_data.path{5} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/40/'; dy_data.codename{5} = "40";
dy_data.path{6} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/50/'; dy_data.codename{6} = "50";
dy_data.path{7} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/60/'; dy_data.codename{7} = "60";
dy_data.path{8} = '/oak/stanford/schools/ees/aditis2/mborrus/dycore/dycore/70/'; dy_data.codename{8} = "70";
dy_data.path{9} = '/scratch/users/mborrus/dycore/spin_up/10/'; dy_data.codename{9} = "10s";
dy_data.path{10} = '/scratch/users/mborrus/dycore/spin_up/20/'; dy_data.codename{10} = "20s";
dy_data.path{11} = '/scratch/users/mborrus/dycore/spin_up/30/'; dy_data.codename{11} = "30s";
dy_data.path{12} = '/scratch/users/mborrus/dycore/spin_up/40/'; dy_data.codename{12} = "40s";
dy_data.path{13} = '/scratch/users/mborrus/dycore/spin_up/50/'; dy_data.codename{13} = "50s";
dy_data.path{14} = '/scratch/users/mborrus/dycore/spin_up/60/'; dy_data.codename{14} = "60s";
dy_data.path{15} = '/scratch/users/mborrus/dycore/spin_up/70/'; dy_data.codename{15} = "70s";
dy_data.path{16} = '/scratch/users/mborrus/dycore/spin_up/80/'; dy_data.codename{16} = "80s";
dy_data.path{17} = '/scratch/users/mborrus/dycore/spin_up/90/'; dy_data.codename{17} = "90s";

dy_data.run_type{1} = 'Base';      dy_data.File_Numbers{1} = [0:5];  dy_data.Color{1}= ["#5e3c99"]; dy_data.Line_Style{1}='-';
dy_data.run_type{2} = '10';      dy_data.File_Numbers{2} = [0:2];  dy_data.Color{2}= ["#4575b4"]; dy_data.Line_Style{2}='-';
dy_data.run_type{3} = '20';      dy_data.File_Numbers{3} = [0:2];  dy_data.Color{3}= ["#74add1"]; dy_data.Line_Style{3}='-';
dy_data.run_type{4} = '30';      dy_data.File_Numbers{4} = [0:2];  dy_data.Color{4}= ["#abd9e9"]; dy_data.Line_Style{4}='-';
dy_data.run_type{5} = '40';      dy_data.File_Numbers{5} = [0:2];  dy_data.Color{5}= ["#fee090"]; dy_data.Line_Style{5}='-';
dy_data.run_type{6} = '50';      dy_data.File_Numbers{6} = [0:2];  dy_data.Color{6}= ["#fdae61"]; dy_data.Line_Style{6}='-';
dy_data.run_type{7} = '60';      dy_data.File_Numbers{7} = [0:2];  dy_data.Color{7}= ["#f46d43"]; dy_data.Line_Style{7}='-';
dy_data.run_type{8} = '70';      dy_data.File_Numbers{8} = [0:2];  dy_data.Color{8}= ["#d73027"]; dy_data.Line_Style{8}='-';
dy_data.run_type{9} = '10s';      dy_data.File_Numbers{9} = [0:9];  dy_data.Color{9}= ["#d73027"]; dy_data.Line_Style{9}='-';
dy_data.run_type{10} = '20s';      dy_data.File_Numbers{10} = [0:9];  dy_data.Color{10}= ["#d73027"]; dy_data.Line_Style{10}='-';
dy_data.run_type{11} = '30s';      dy_data.File_Numbers{11} = [0:9];  dy_data.Color{11}= ["#d73027"]; dy_data.Line_Style{11}='-';
dy_data.run_type{12} = '40s';      dy_data.File_Numbers{12} = [0:9];  dy_data.Color{12}= ["#d73027"]; dy_data.Line_Style{12}='-';
dy_data.run_type{13} = '50s';      dy_data.File_Numbers{13} = [0:9];  dy_data.Color{13}= ["#d73027"]; dy_data.Line_Style{13}='-';
dy_data.run_type{14} = '60s';      dy_data.File_Numbers{14} = [0:9];  dy_data.Color{14}= ["#d73027"]; dy_data.Line_Style{14}='-';
dy_data.run_type{15} = '70s';      dy_data.File_Numbers{15} = [0:9];  dy_data.Color{15}= ["#d73027"]; dy_data.Line_Style{15}='-';
dy_data.run_type{16} = '80s';      dy_data.File_Numbers{16} = [1 2 3 4 5 6 8 9];  dy_data.Color{16}= ["#d73027"]; dy_data.Line_Style{16}='-';
dy_data.run_type{17} = '90s';      dy_data.File_Numbers{17} = [0:9];  dy_data.Color{17}= ["#d73027"]; dy_data.Line_Style{17}='-';

%Lat Ranges
dy_data.lat_range_names{1} = "North High";
dy_data.lat_range_names{2} = "South High";
dy_data.lat_range_names{3} = "North Mid";
dy_data.lat_range_names{4} = "South Mid";
dy_data.lat_range_names{5} = "Equator";
dy_data.lat_range_names{6} = "45 deg N";
%Pressure Ranges
dy_data.p_range_names{1} = "Lower Tropo";  dy_data.p_color{1} = "#a1dab4";
dy_data.p_range_names{2} = "Upper Tropo"; dy_data.p_color{2} = "#41b6c4";
dy_data.p_range_names{3} = "Stratosphere";  dy_data.p_color{3} = "#225ea8";


% N_Hi = find(lat_ranges >  50); dy_data.lat_range{1}=N_Hi;
% S_Hi = find(lat_ranges < -50); dy_data.lat_range{2}=S_Hi;
% N_Mid = find(lat_ranges > 20 & lat_ranges < 50); dy_data.lat_range{3}=N_Mid;
% S_Mid = find(lat_ranges < -20 & lat_ranges > -50); dy_data.lat_range{4}=S_Mid;
% Equator = find(lat_ranges > -20 & lat_ranges < 20); dy_data.lat_range{5}=Equator;
% N_45 = find(lat_ranges > 42 & lat_ranges < 48); dy_data.lat_range{6}=N_45;
%%
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('~/Matlab_HPC/scripts/EGR/')
for run_selection = 1:8;
    clear lambda_term T u lat p theta dtheta_dz dtheta_dp du_z DRY_N2 DRY_egr MOIST_N2 MOIST_egr dtheta_dz_eff dtheta_dp_eff
    File_Numbers=dy_data.File_Numbers{run_selection};
    dy_data_Path = dy_data.path{run_selection};
    temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(1)),'/u_interp_01.mat');
    pressure_levels = dy_data.p;
    lat = dy_data.lat;
    days = (1:100);
    Pressure_range = find(pressure_levels<850 & pressure_levels>200);
    p = pressure_levels(Pressure_range);
    
    H = 7300; % scale height, m
    z = -H*log(p/1000);
    omega = 2*pi/(24*3600);
    f = 2*omega*sin(2*pi*lat/360);
    g = 9.8;
    
    %lambda_input = interp1(lambda_base(:,1),lambda_base(:,2),lat);
    
    for run_N = 1:length(File_Numbers)
        %for run_N = 1
        temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(run_N)),'/u_interp_01.mat');
        
        T = load(strcat(dy_data.path{run_selection},num2str(File_Numbers(run_N)),'/T_interp_01.mat'),'T_interp_01');
        u = load(strcat(dy_data.path{run_selection},num2str(File_Numbers(run_N)),'/u_interp_01.mat'),'u_interp_01');
        %%
        %cut off the top 6 layers of atmosphere
        T  = T.T_interp_01(:,:,Pressure_range,:);
        u  = u.u_interp_01(:,:,Pressure_range,:);
        %% Average across pressure and longitude
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
            lambda = 0;
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
    
%     dy_data.dry_egr{run_selection} = DRY_egr;
%     dy_data.moist_egr{run_selection} = MOIST_egr;
%     dy_data.dry_N2{run_selection} = DRY_N2;
%     dy_data.moist_N2{run_selection} = MOIST_N2;
    dy_data.duz{run_selection} = du_z;
    "saved"
    
end
save(['dy_data.mat'],'dy_data')

for run_selection = [1:17];
    if run_selection > 0;
        for i = 1:6
            dry_means(i)=nanmean(dy_data.dry_egr{run_selection}(:,dy_data.lat_range{i},:,:),[1,2,3,4]);
        end
        dy_data.dry_egr_mean{run_selection} = dry_means;
        clear means moist_means dry_means
    end
end

 
%% Grab the U values from each folder
for run_selection=[10:17]
    File_Numbers=dy_data.File_Numbers{run_selection};
    temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(1)),'/u_interp_01.mat');
    temp_u = load(temp_path,'u_interp_01');
    [NNlon, NNlat,NNp,Ndays] = size(temp_u.u_interp_01);
    U_0 = temp_u.u_interp_01(:,:,whole_P_range,Ndays);
    % Get the pressure and lats values
%     for i = 0:length(File_Numbers)-1
    for i = 0:3
        if i == 0
            var_diff = zeros(length(File_Numbers)-1,Nlon,Nlat,Np,Ndays);
            var_long = zeros(length(File_Numbers),Nlon,Nlat,Np,Ndays);
            temp_var = load(temp_path,'u_interp_01');
            temp_var = temp_var.u_interp_01(:,:,whole_P_range,1:Ndays);
            var_1 = temp_var(:,:,:,1:Ndays);
            var_long(i+1,:,:,:,:) = temp_var(:,:,:,1:Ndays);
            "vars made"
        else
            temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(i)),'/u_interp_01.mat');
            temp_var = load(temp_path,'u_interp_01');
            temp_var = temp_var.u_interp_01(:,:,whole_P_range,1:Ndays);
            var_diff(i,:,:,:,:) = temp_var(:,:,:,1:Ndays) - var_1;
            var_long(i+1,:,:,:,:) = temp_var(:,:,:,1:Ndays); i
            "vars filled, not made"
        end
        temp_path
    end
    clear var_1 temp_var
    
    [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff); clear var_RMS var_Mean var_STD
    for lats = 6
        for pres=1:3
            range_temp = dy_data.lat_range{lats}; rng = length(range_temp); p_range_temp = dy_data.p_range{pres};
            var_All_RMS(pres,lats,:,:) = squeeze(rms(reshape(var_diff(:,:,range_temp,p_range_temp,:), Nrun,rng*Nlon*length(p_range_temp),Ntime),2));
            var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range_temp,:), Nrun,rng*Nlon*length(p_range_temp),Ntime),2),1));
            var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range_temp,:),[1,2,3,4]));
            var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range_temp,:),1,[1,2,3,4],'omitnan'));
        end
    end
    dy_data.U_All_RMS{run_selection} = var_All_RMS;
    dy_data.U_RMS{run_selection} = var_RMS;
    dy_data.U_Mean{run_selection} = var_Mean;
    dy_data.U_STD{run_selection} = var_STD;
    save(['dy_data.mat'],'dy_data');
    clear var_diff var_long var_RMS var_Mean var_STD var_All_RMS
    %%%
    temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(1)),'/T_interp_01.mat');
    % Get the pressure and lats values
    p_range{1} = find(pressure_levels<850 & pressure_levels>500);
    p_range{2} = find(pressure_levels<500 & pressure_levels>100);
    p_range{3} = find(pressure_levels<100 & pressure_levels>1);
    for i = 0:length(File_Numbers)-1
        if i == 0
            var_diff = zeros(length(File_Numbers)-1,Nlon,Nlat,Np,Ndays);
            var_long = zeros(length(File_Numbers),Nlon,Nlat,Np,Ndays);
            temp_var = load(temp_path,'T_interp_01');
            temp_var = temp_var.T_interp_01(:,:,whole_P_range,1:Ndays);
            var_1 = temp_var(:,:,:,1:Ndays);
            var_long(i+1,:,:,:,:) = temp_var(:,:,:,1:Ndays);
        else
            temp_path = strcat(dy_data.path{run_selection},num2str(File_Numbers(i)),'/T_interp_01.mat');
            temp_var = load(temp_path,'T_interp_01');
            temp_var = temp_var.T_interp_01(:,:,whole_P_range,1:Ndays);
            var_diff(i,:,:,:,:) = temp_var(:,:,:,1:Ndays) - var_1;
            var_long(i+1,:,:,:,:) = temp_var(:,:,:,1:Ndays); i
        end
    end
    clear var_1 temp_var
    
    [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff); clear var_RMS var_Mean var_STD
     for lats = 1:6
        for pres=1:3
            range_temp = dy_data.lat_range{lats}; rng = length(range_temp); p_range_temp = dy_data.p_range{pres};
            var_All_RMS(pres,lats,:,:) = squeeze(rms(reshape(var_diff(:,:,range_temp,p_range_temp,:), Nrun,rng*Nlon*length(p_range_temp),Ntime),2));
            var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range_temp,:), Nrun,rng*Nlon*length(p_range_temp),Ntime),2),1));
            var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range_temp,:),[1,2,3,4]));
            var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range_temp,:),1,[1,2,3,4],'omitnan'));
        end
    end
    dy_data.T_All_RMS{run_selection} = var_All_RMS;
    dy_data.T_RMS{run_selection} = var_RMS;
    dy_data.T_Mean{run_selection} = var_Mean;
    dy_data.T_STD{run_selection} = var_STD;
    save(['dy_data.mat'],'dy_data');
    clear var_diff var_long var_RMS var_Mean var_STD
    
    %%
    for lat_band = [6]
        for i = 1:3
            a =squeeze(dy_data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
            b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
        end
    end
    sat_times(:,lat_band)=b_sat(:) ;
    %%

    %%
    %for run_selection = 1
    run_selection =1
    sat_times = zeros(3,6);
    for lat_band = [6];
        savepath = strcat('./plots/Master/RMS/dycore/',dy_data.codename{run_selection},'/'); mkdir(savepath)
        run_name = "dycore h0";
        U_RMS_Plot = figure(2); clf, hold on, clear a b plot_rms b_sat
        for i = 1:2;
            a =squeeze(dy_data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10); c=mean(a);d=find(a>c/2,1);
            plot_rms(i) = plot(1:400, a, 'Color', dy_data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
            plot_rms(i) = plot(1:400, b, 'Color', dy_data.p_color{i},'LineWidth',2);
            e=smooth((diff(b(d:end))./b(d:end-1)),10);
            f=find(e<.01,2)+(d);
            %plot((d+1:400),e*100)
            b_sat(i) = f(1);
            plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            xline(b_sat(i))
            i;
            for j = 1:2
                plot(1:300,squeeze(dy_data.U_All_RMS{run_selection}(2,lat_band,j,1:300)))
            end
        end
        hleg = legend(plot_rms(1:i),dy_data.p_range_names{1:i});
        axis([30 300 -inf inf]), title([dy_data.lat_range_names{lat_band},' RMS : U',run_name]), xlabel('Days'), ylabel('RMSE')
        saveas(U_RMS_Plot,strcat(savepath,'_U_',dy_data.lat_range_names{lat_band},'.png'))
        sat_times(1:i,lat_band)=b_sat(1:i) ;
        run_selection
    end
    dy_data.U_Sat_Time{run_selection}=sat_times;
    sat_times
    %end
    
    save(['dy_data.mat'],'dy_data')
    %%
    sat_times = zeros(3,6);
    for lat_band = [6];
        savepath = strcat('./plots/Master/RMS/dycore/',dy_data.codename{run_selection},'/'); mkdir(savepath)
        run_name = "dycore h0";
        T_RMS_Plot = figure(3); clf, hold on, clear a b plot_rms b_sat
        for i = 1:2;
            a =squeeze(dy_data.T_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10); c=mean(a);d=find(a>c/2,1);
            plot_rms(i) = plot(1:400, a, 'Color', dy_data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
            plot_rms(i) = plot(1:400, b, 'Color', dy_data.p_color{i},'LineWidth',2);
            e=smooth((diff(b(d:end))./b(d:end-1)),10);
            f=find(e<0.001,2)+(d);
            %plot((d+1:400),e*100)
            b_sat(i) = f(2);
            plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            xline(b_sat(i))
            i;
            for j = 1:5
                plot(1:100,squeeze(dy_data.T_All_RMS{1}(2,lat_band,j,1:100)))
            end
        end
        hleg = legend(plot_rms(:),dy_data.p_range_names{:});
        axis([1 100 -5 20]), title([dy_data.lat_range_names{lat_band},' RMS : T',run_name]), xlabel('Days'), ylabel('RMSE')
%         saveas(T_RMS_Plot,strcat(savepath,'_T_',dy_data.lat_range_names{lat_band},'.png'))
        sat_times(1:i,lat_band)=b_sat(1:i) ;
        run_selection
    end
    dy_data.T_Sat_Time{run_selection}=sat_times;
    
    save(['dy_data.mat'],'dy_data')
end

%%
for run_selection = 1:8
    lat_band = [6]
    for i = 1:3;
        a =squeeze(dy_data.T_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
    end
end






%%
run_selection = 1
for lat_band = [1:6];
    savepath = strcat('./plots/Master/RMS/AM4+dycore/'); mkdir(savepath)
    run_name = "dycore h0";
    U_duo_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat
    for i = 1:3;
        a =squeeze(dy_data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
        plot_rms(i) = plot(1:100, a, 'Color', dy_data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
        plot_rms(i) = plot(1:100, b, 'Color', dy_data.p_color{i},'LineWidth',2);
        b_sat(i) = find(diff(b(50:end))./b(51:end)<.03,1)+50;
        plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
        i; legname(i) = "dycore " + dy_data.p_range_names{i} ;
    end
    for i = 4:6;
        a =squeeze(dy_data.U_RMS{run_selection}(i-3,lat_band,:)); b = smooth(a,10);
        plot_rms(i) = plot(1:100, a,'--', 'Color', dy_data.p_color{i-3},'LineWidth',2); plot_rms(i).Color(4)=.3;
        plot_rms(i) = plot(1:100, b,'--', 'Color', dy_data.p_color{i-3},'LineWidth',2);
        b_sat(i) = find(diff(b(1:end))./b(2:end)<.03,1)+1;
        plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
        i; legname(i) = "AM4 " + dy_data.p_range_names{i-3} ;
    end
    hleg = legend(plot_rms(:),legname(:));
    axis([1 100 0 75]), title([dy_data.lat_range_names{lat_band},' RMS : U','AM4 and Dycore']), xlabel('Days'), ylabel('RMSE')
    saveas(U_duo_RMS_Plot,strcat(savepath,'_U_',dy_data.lat_range_names{lat_band},'.png'))
    run_selection
end
























%% OLD STUFF
% ucomp = Size:       128x64x40x400
% Dimensions: lon, lat, P, days

% Get initial values so we don't need to save more than we need
temp_path = strcat(Dycore_Data_Path,'ary',num2str(0),'/u_interp_01.mat');
temp_u = load(temp_path,'u_interp_01');
U_0 = temp_u.u_interp_01(:,:,whole_P_range,1:100);

clear temp_path temp_u
U_diff = zeros(20,Nlon,Nlat,Np,Ntime);
for i = 1:Nrun
    
    % The files are spread out in 20 different folders, this gets you to each of them
    temp_path = strcat(Dycore_Data_Path,'ary',num2str(i),'/u_interp_01.mat');
    temp_u = load(temp_path,'u_interp_01');
    temp_u = temp_u.u_interp_01(:,:,whole_P_range,1:100);
    
    % subtract temp_u
    U_diff(i,:,:,:,:) = temp_u - U_0;
end
clear U_1 temp_u

N_low_T = length(Lower_Tropo);
N_Up_T = length(Upper_Tropo);
N_Strat = length(Stratosphere);

%%
c = 'Color';
c1 = ["#a1dab4","#41b6c4","#225ea8"];

% GLOBAL

range_temp = 1:Nlat; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));

MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('Global RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/Global_h0.png'])

% N_Hi
clear U_LT U_UT U_S
range_temp = N_Hi; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('North Hi-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/North_Hi_h0.png'])

% N_Mid
clear U_LT U_UT U_S
range_temp = N_Mid; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('North Mid-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/North_Mid_h0.png'])

% Low_lat
clear U_LT U_UT U_S
range_temp = Equator; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('Low-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/Low_h0.png'])

% Mid_S_lat
clear U_LT U_UT U_S
range_temp = S_Mid; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('South Mid-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/South_Mid_h0.png'])

% High_S_lat
clear U_LT U_UT U_S
range_temp = S_Hi; rng = length(range_temp);
U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
MainPlot = figure(1)
clf, hold on
plot_lt = plot(days, U_LT, c, c1(1));
plot_ut = plot(days, U_UT, c, c1(2));
plot_s = plot(days, U_S , c, c1(3));
plots = [plot_lt(1), plot_ut(1), plot_s(1)];
hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

axis([1 100 0 75]), title('South Hi-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

saveas(MainPlot,['./plots/RMS/South_Hi_h0.png'])


%%
clear dydryegrmean dymoistegrmean dydryNmean dymoistNmean dydryegrrecalc dydudzmean
lat = dy_data.lat;
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat(dy_data.lat_range{6})/360);
daterange = (1:85);
 i=1;
 for run_selection = [2 3 4 5 6 7 8];
     dydryegrmean(:,i) = mean(dy_data.dry_egr{run_selection}(:,dy_data.lat_range{6},:,daterange),[1 3 4])*86400;
     dydryNmean(:,i) = abs(mean(dy_data.dry_N2{run_selection}(:,dy_data.lat_range{6},:,daterange).^(1/2),[1 3 4])*86400);
     dydudzmean(:,i) = mean(dy_data.duz{run_selection}(:,dy_data.lat_range{6},:,daterange),[1 3 4])*86400;
     dydryegrrecalc(:,i) = mean(abs(f.*(dy_data.duz{run_selection}(:,dy_data.lat_range{6},:,daterange))./(dy_data.dry_N2{run_selection}(:,dy_data.lat_range{6},:,daterange)).^(1/2)),[1 3 4])*86400;
     i = i+1
 end
mean(dydryegrmean)'
mean(dydryNmean)'
mean(dydudzmean)'

dybackcalcegr = f.*(dydudzmean./dydryNmean)';
mean(dybackcalcegr,2)*86400 - mean(dydryegrmean)'
mean(dydryegrrecalc,1)' - mean(dydryegrmean)'

mean(abs(f.*(dy_data.duz{run_selection}(:,dy_data.lat_range{6},:,daterange))./(dy_data.dry_N2{run_selection}(:,dy_data.lat_range{6},:,daterange)).^(1/2)),[1 2 3 4])*86400
mean(f.*(dy_data.duz{run_selection}(:,dy_data.lat_range{6},:,daterange))./(abs(dy_data.dry_N2{run_selection}(:,dy_data.lat_range{6},:,daterange)).^(1/2)),[1 2 3 4])*86400


egrtest = ...
    abs(f .* ...
    ( dy_data.duz{run_selection}( : , dy_data.lat_range{6} , : , daterange ) ) ...
    ./ ...
    abs( dy_data.dry_N2{run_selection}( : , dy_data.lat_range{6} , : , daterange ).^( 1/2 ) ));...

dy_data.p(dy_data.p_range{1})
dy_data.p(dy_data.p_range{2})

for plevels = 1:11;
    for days = 1:400
        negs(plevels,days) = length(find(dy_data.dry_N2{3}(:,dy_data.lat_range{6},plevels,days)<0));
    end
end
negs(6,:)
