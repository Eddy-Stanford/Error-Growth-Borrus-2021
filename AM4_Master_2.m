% AM4_Master_2.0
% AM4 documentation for the second iteration of the AM4 runs

clear
cd '/home/users/mborrus/Matlab_HPC'
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('/home/users/mborrus/Matlab_HPC/plots')
addpath('/home/users/mborrus/Matlab_HPC/scripts/cbrewer')
load('./data/EGR/lambda.mat')
load('AM4_Data_v2.mat')

%% 1.2 Data Loading - Struct Creation

% Create the data paths for AM4 data
    SCRATCHPATH = '/scratch/users/mborrus/';
    AM4_Data.path{1} = strcat(SCRATCHPATH, 'AM4/SST_m16/');
    AM4_Data.path{2} = strcat(SCRATCHPATH, 'AM4/SST_m12/');
    AM4_Data.path{3} = strcat(SCRATCHPATH, 'AM4/SST_m8/');
    AM4_Data.path{4} = strcat(SCRATCHPATH, 'AM4/SST_m4/');
    AM4_Data.path{5} = strcat(SCRATCHPATH, 'AM4/SST_m3/');
    AM4_Data.path{6} = strcat(SCRATCHPATH, 'AM4/SST_m2/');
    AM4_Data.path{7} = strcat(SCRATCHPATH, 'AM4/SST_m1/');
    
    AM4_Data.path{8} = strcat(SCRATCHPATH, 'AM4/BASE/');
    
    AM4_Data.path{9} = strcat(SCRATCHPATH, 'AM4/SST_p1/');
    AM4_Data.path{10} = strcat(SCRATCHPATH, 'AM4/SST_p2/');
    AM4_Data.path{11} = strcat(SCRATCHPATH, 'AM4/SST_p3/');
    AM4_Data.path{12} = strcat(SCRATCHPATH, 'AM4/SST_p4/');
    AM4_Data.path{13} = strcat(SCRATCHPATH, 'AM4/SST_p8/');
    AM4_Data.path{14} = strcat(SCRATCHPATH, 'AM4/SST_p12/');
    AM4_Data.path{15} = strcat(SCRATCHPATH, 'AM4/SST_p16/');

% Create run names, # of files, title names, and plot colors    
    idx = 1; AM4_Data.run_type{idx} = '-16C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-16C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 2; AM4_Data.run_type{idx} = '-12C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-12C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 3; AM4_Data.run_type{idx} = '-8C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-8C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 4; AM4_Data.run_type{idx} = '-4C';      AM4_Data.File_Numbers{idx} = [1:9]; AM4_Data.codename{idx} = '-4C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 5; AM4_Data.run_type{idx} = '-3C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-3C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 6; AM4_Data.run_type{idx} = '-2C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-2C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 7; AM4_Data.run_type{idx} = '-1C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '-1C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    
    idx = 8; AM4_Data.run_type{idx} = '+0C Base';      AM4_Data.File_Numbers{idx} = [1:13]; AM4_Data.codename{idx} = 'base'; AM4_Data.Color{idx}= [0.6980    0.8745    0.5412; 0.2000    0.6275    0.1725]; AM4_Data.Line_Style{idx}='-';
  
    idx = 15; AM4_Data.run_type{idx} = '+16C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '+16C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 14; AM4_Data.run_type{idx} = '+12C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '+12C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 13; AM4_Data.run_type{idx} = '+8C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '+8C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 12; AM4_Data.run_type{idx} = '+4C';      AM4_Data.File_Numbers{idx} = [1:9]; AM4_Data.codename{idx} = '+4C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 11; AM4_Data.run_type{idx} = '+3C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '+3C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 10; AM4_Data.run_type{idx} = '+2C';      AM4_Data.File_Numbers{idx} = [1:6]; AM4_Data.codename{idx} = '+2C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    idx = 9; AM4_Data.run_type{idx} = '+1C';      AM4_Data.File_Numbers{idx} = [1:9]; AM4_Data.codename{idx} = '+1C'; AM4_Data.Color{idx}= [0    0    0; 0    0    0]; AM4_Data.Line_Style{idx}='-';
    
    temp_path = strcat(AM4_Data.path{1},num2str(1),'/dailyUTPS.nc');
    AM4_Data.p = ncread(temp_path,'pfull');
    AM4_Data.lat = ncread(temp_path,'grid_yt');
% Anxillary data: Lat ranges, pressure ranges 
    %Lat Ranges
    AM4_Data.lat_range_names{6} = "45 deg N"; AM4_Data.lat_range{6} = find(AM4_Data.lat > 42.5 & AM4_Data.lat < 47.5);
    
    %Pressure Ranges
    AM4_Data.p_range_names{1} = "Lower Tropo"; AM4_Data.p_range{1} = find(AM4_Data.p<850 & AM4_Data.p>500); AM4_Data.p_color{1} = "#a1dab4";
    AM4_Data.p_range_names{2} = "Upper Tropo"; AM4_Data.p_range{2} = find(AM4_Data.p<500 & AM4_Data.p>100); AM4_Data.p_color{2} = "#41b6c4";
    AM4_Data.p_range_names{3} = "Stratosphere"; AM4_Data.p_range{3} = find(AM4_Data.p<100 & AM4_Data.p>1); AM4_Data.p_color{3} = "#225ea8";
    
    save(['/scratch/users/mborrus/AM4/AM4_Data_v2.mat'],'AM4_Data')
    
%%%%%%%%%%%%%%%%%%%%%
%%
% ask for a lot of memory before you do this --mem=32G
% I had a max file size of 256G - sheesh!
a=zeros([15, 20, 144,    90,    33-9,   50]);
for run_selection = [1:2];
    File_Numbers=AM4_Data.File_Numbers{run_selection};
    for i=File_Numbers
        AM4_Data_Path = AM4_Data.path{run_selection};
        startLoc = [1 1 1 1]; % Start location along each coordinate
        count  = [Inf Inf 24 50]; % Read until the end of each dimension
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTPS.nc');
        a(run_selection,i,:,:,:,:)=ncread(temp_path,'ucomp',startLoc,count);
        i
    end
    run_selection
    AM4_Data.codename{run_selection}
end

save(['/scratch/users/mborrus/AM4/UcompleteALL.mat'],'a', '-v7.3');
clear umasterBASE

  %%
  tic
load  /scratch/users/mborrus/AM4/UcompleteALL.mat
  toc
  %%
  
cd '/home/users/mborrus/Matlab_HPC'
load('/scratch/users/mborrus/AM4/AM4_Data_v2.mat'); tic;
load('/scratch/users/mborrus/AM4/UcompleteALL.mat'); toc;

% size = [21     5   144    90    24   100]
% size = [Run   Int  Lon    Lat   p    day]

%% calculate the RMS for each run by cycling through the base run (doesn't calc sat time)
[Ntemps,Nrun,Nlon,Nlat,Np,Ntime]=size(a);

for run_selection=[4 7 12]
    Nrun=length(AM4_Data.File_Numbers{run_selection});
    var_diff=zeros((Nrun-1)*Nrun/2,Nlon,Nlat,Np,Ntime);
    k=1
    q=1
    disp("variable creation")
    disp(AM4_Data.codename{run_selection})
    while k<length(AM4_Data.File_Numbers{run_selection})
            
        for run = k:length(AM4_Data.File_Numbers{run_selection})-1
            var_diff(q,:,:,:,:) = a(run_selection,k,:,:,:,:)-a(run_selection,run+1,:,:,:,:);
            [q, k,run+1]
            q=q+1;           
        end
        k=k+1;
    end
    %%
    %clear a
    %%
    %90 lats = 88 bins
    [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff);
    disp("RMS calculations")
    %% Get RMS
    var_RMS=NaN(Nrun,88,3,Ntime);
    for lat_iteration = 1:88
        for pres=1:3
            range_temp = lat_iteration:2+lat_iteration; rng = length(range_temp);
            
            var_RMS(:,lat_iteration,pres,:) = squeeze(...
                rms(...
                reshape(...
                var_diff(:,:,range_temp,AM4_Data.p_range{pres},:),...
                 Nrun,   rng * Nlon * length(AM4_Data.p_range{pres}),   Ntime)...
                ,2)...
                );
        end
    end
    RMS_storage.U{run_selection}=var_RMS;
    
end

%% Method 2:  calculate the RMS for each run without cycling through the base run (doesn't calc sat time)
 
cd '/home/users/mborrus/Matlab_HPC'
load('/scratch/users/mborrus/AM4/AM4_Data_v2.mat'); tic;
load('/scratch/users/mborrus/AM4/UcompleteALL.mat'); toc;
load('/scratch/users/mborrus/AM4/RMS_storage_v2.mat');

%%
[Ntemps,Nrun,Nlon,Nlat,Np,Ntime]=size(a);

for run_selection=[4 7 12]
    Nrun=length(AM4_Data.File_Numbers{run_selection});
    var_diff=zeros(Nrun,Nrun-1,Nlon,Nlat,Np,Ntime);
    
    baseline_run=1
    difference_run=1
    indexes = AM4_Data.File_Numbers{run_selection};

    disp("variable creation")
    disp(AM4_Data.codename{run_selection})
    
    for run = indexes
        run_list = circshift(indexes,-run)
        var_diff(run,:,:,:,:,:) = a(run_selection,run_list(1:end-1),:,:,:,:)-a(run_selection,run,:,:,:,:);
    end
    
    %90 lats = 88 bins
    [Nrun,Ndiffs,Nlon,Nlat,Np,Ntime] = size(var_diff);
    disp("RMS calculations")
    %% Get RMS
    var_RMS=NaN(Nrun,Ndiffs,88,3,Ntime);
    for lat_iteration = 1:88
        for pres=1:3
            range_temp = lat_iteration:2+lat_iteration; rng = length(range_temp);
            
            var_RMS(:,:,lat_iteration,pres,:) = squeeze(...
                rms(...
                reshape(...
                var_diff(:,:,:,range_temp,AM4_Data.p_range{pres},:),...
                 Nrun, Ndiffs,   rng * Nlon * length(AM4_Data.p_range{pres}),   Ntime)...
                ,3)...
                );
        end
    end
    RMS_storage.U{run_selection}=var_RMS;
    
end
%%
save(['/scratch/users/mborrus/AM4/RMS_storage_v2.mat'],'RMS_storage', '-v7.3');
%% 
load('/scratch/users/mborrus/AM4/RMS_storage_v2.mat');
%% Method 1: This calculates saturation time based on the derivative of a single run,
% whereas in the paper we do it over an average
Nlat=90;

for run_selection=[1:15]
    clear sat_time a c
    a = smoothdata(RMS_storage.U{run_selection},4,'movmean',7);
    c = diff(a,1,4)./a(:,:,:,2:end); 
    for i=1:size(c,1)
        for j=1:Nlat-2
            for k=1:3
                sat_time(i,j,k)=find(c(i,j,k,:)<.03,1)+1;
            end
        end
    end
    s(run_selection)=std(sat_time(:,67,2));
    m(run_selection)=mean(sat_time(:,67,2));
    I = sat_time(:,67,2)<m(run_selection)+s(run_selection) & ...
        sat_time(:,67,2)>m(run_selection)-s(run_selection);
    s_filt(run_selection)=std(sat_time(I,67,2));
    m_filt(run_selection)=mean(sat_time(I,67,2));
end

%%
m(8)
m_filt(8)
(m-m(8))'
(m_filt-m_filt(8))'

s(8)
s_filt(8)
(s-s(8))'
(s_filt-s_filt(8))'

%% Method 2: This calculates the mean based on the mean of all the runs 
% difference is in the third line, where it gets averaged 

for run_selection=[1:15]
    clear sat_time
    a = smoothdata(nanmean(RMS_storage.U{run_selection},1),4,'movmean',2);
    c = diff(a,1,4)./a(:,:,:,2:end); 
    for i=1:size(c,1)
        for j=1:Nlat-2
            for k=1:3
                sat_time(i,j,k)=find(c(i,j,k,:)<.02,1)+1;
            end
        end
    end
    m2(run_selection)=mean(sat_time(:,67,2),1);
end

m2(8)
(m2-m2(8))'

%%
start_in(1) = [1];
end_in(1) = [19];
for i = 2:19
start_in(i) = end_in(i-1)+1;
end_in(i) = start_in(i)+(19-i);
end
%% Method 3: This allowed me to test various smoothing windows 

clear m_compare m_window std_window
for window = 1:20
    for indexes = 1:13
        a = smoothdata(nanmean(RMS_storage.U{8}(start_in(indexes):end_in(indexes),:,:,:),1),4,'movmean',window);
        c = diff(a,1,4)./a(:,:,:,2:end); 
        m_compare(indexes) = find(c(1,67,2,:)<.03,1)+1;
        a_stored(indexes,:,:,:,:) = a;
    end
        m_window(window) = mean(m_compare);
        std_window(window) = std(m_compare);
end
clf; plot([1:window],m_window); hold on
plot([1:window],[m_window+std_window]);
plot([1:window],[m_window-std_window]);
m_window(12)
std_window(12)
legend('mean','std')

%% Method 4: This calculates the mean each baseline seperately (all for 1, then all for 2, etc...)
% will provide as many means as we have
    %clear sat_time;
    sat_time=nan([15    13    88     3]);
for run_selection=[1:15]
    a = smoothdata(squeeze(nanmean(RMS_storage.U{run_selection},2)),4,'movmean',7);    
    c = diff(a,1,4)./a(:,:,:,2:end); 
    for i=1:size(c,1)
        for j=1:90-2
            for k=1:3
                sat_time(run_selection,i,j,k)=find(c(i,j,k,:)<.03,1)+1;
            end
        end
    end
end
m3=nanmean(sat_time(:,:,67,2),2)
%s3=std(sat_time(:,:,67,2),2);

m3(8)
(m3-m3(8))

%% This plots a scatter plot with prediction intervals - averaging
%load T_means.mat
baseline_temp = nanmean(T_mean(8,:,67,2),2);
temps_basic = squeeze(nanmean(T_mean(:,:,67,2),2))-baseline_temp;
sat_times_basic = m3-m3(8);
x=temps_basic;
y=sat_times_basic';

figure(1),clf,hold on
[p,S] = polyfit(x,y,1); 
[y_fit,delta] = polyval(p,x,S);
plot(x,y,'bo')
hold on
plot(x,y_fit,'r-')
plot(x,y_fit+1.15*delta,'m--',x,y_fit-1.15*delta,'m--')
title('Linear Fit of Data with 75% Prediction Interval')
legend('Data','Linear Fit','75% Prediction Interval')
%errorbar(x,y,s3, 'LineStyle','none')
%% This plots the same as above but doesn't average by temp experiement
latband = 67;
baseline_temp = nanmean(T_mean(8,:,latband,2),2);
xx = T_mean(:,1:13,latband,2)-baseline_temp;
x = reshape(xx,[1,size(xx,1)*size(xx,2)]);

baseline_sat = nanmean(sat_time(8,:,latband,2),2);
yy = sat_time(:,:,latband,2)-baseline_sat;
y = reshape(yy,[1,size(yy,1)*size(yy,2)]);

mask = isnan(y);

x=x(~mask)
y=y(~mask)

figure(2),clf,hold on
plot(x,y,'o')
[p,S] = polyfit(x,y,1); 
[y_fit,delta] = polyval(p,x,S);
plot(x,y,'bo')
hold on
plot(x,y_fit,'r-')
plot(x,y_fit+1.15*delta,'m--',x,y_fit-1.15*delta,'m--')
%axis([-7.5 17.5 -5 5])
%% Create a contour plot of the saturation times 
for run_selection=[1:15]
    a = smoothdata(squeeze(nanmean(RMS_storage.U{run_selection},2)),4,'movmean',7);    
    c = diff(a,1,4)./a(:,:,:,2:end); 
    for i=1:size(c,1)
        for j=1:90-2
            for k=1:3
                sat_time(run_selection,i,j,k)=find(c(i,j,k,:)<.03,1)+1;
            end
        end
    end
end

pressure_range_selection=2;
runs = [1:15];
baseline_temp = squeeze(nanmean(T_mean(8,:,:,2),2));   %take the average T for the base run
baseline_sat  = squeeze(nanmean(sat_time(8,:,:,2),2)); %take the average sat time for the base run
baseline_sat_std  = squeeze(std(sat_time(8,:,:,2),0,[2],'omitnan')); %take the std for the base run

for i = 1:90-2 %for each lat band
    lat_mid(i) = mean(AM4_Data.lat(i:i+2)); %Latitude (deg lat) - 88;
    YYY(:,:,i) = T_mean(:,1:13,i,2)-baseline_temp(i); %Temperature (K) - 15 T ranges x 13 runs x 88 lats;
    ZZZ(:,:,i) = sat_time(:,1:13,i,2)-baseline_sat(i); %Saturation Time (days) - 15 T ranges x 13 runs x 88 lats;
end 

X = lat_mid; %Latitude (deg lat) - 88;
YY = reshape(YYY,[size(YYY,1)*size(YYY,2),88]); % flattened T - (15 x 13 = 195) x 88;
ZZ = reshape(ZZZ,size(ZZZ,1)*size(ZZZ,2),88); % flattened sat - (15 x 13 = 195) x 88;
mask = ~isnan(YY(:,45));  % find the non-nan values
clear Z Y
Y = YY(mask,:); % flattened T, no nans
Z = ZZ(mask,:); % flattened sat, no nans

minT=round(min(min(min(Y))),1) %-10.1
maxT=round(max(max(max(Y))),1) %33   

TempAxis=minT:.5:maxT+.25;
count=zeros(88,length(TempAxis));
SatIntGrid=NaN(length(lat_mid),length(TempAxis));
CountGrid=zeros(length(lat_mid),length(TempAxis));
averaged = 0;
single = 0;
for run = 1:size(Z,1)
    for lat_bin = 1:88
        a = round(Y(run,lat_bin),4);
        index = find(min(abs(a-TempAxis))==abs(a-TempAxis));
        count(lat_bin,index)=count(lat_bin,index)+1;
        if count(lat_bin,index) == 1
            SatIntGrid(lat_bin,index)=Z(run,lat_bin);
            single = single+1;
        else count(lat_bin,index) > 1;
            SatIntGrid(lat_bin,index)=(Z(run,lat_bin)+SatIntGrid(lat_bin,index)*count(lat_bin,index))/count(lat_bin,index);
            averaged = averaged+1;
        end
    end
end
sum(count,'all')

for lat_bin=1:88
    x=SatIntGrid(lat_bin,:);
    nanx = isnan(x);
    t    = 1:numel(x);
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    SatIntGrid(lat_bin,:)=x;
end

f = figure(4);clf,hold on
surface(X,TempAxis,SatIntGrid','FaceColor','interp','EdgeColor','interp')
axis([-80 80 -7.5 7.5])
colorbar
c=colormap('turbo');
colormap(f, flipud(colormap(f)))
caxis([-4 4])
xline([0 -45 45],'--')
yline([0],'--')
xlabel('latitude'); ylabel('T from base');title('T and lat vs Sat Time U'); set(gca,'FontSize',14);

f = figure(5);clf,hold on
contourf(X,TempAxis,SatIntGrid',[-4 -2 0 2 4])
colorbar
colormap('turbo');
colormap(f, flipud(colormap(f)))
axis([-80 80 -7.5 7.5])
caxis([-4 4])
xline([0 -45 45],'--')
yline([0],'--')
xlabel('latitude'); ylabel('T from base');title('T and lat vs Sat Time U'); set(gca,'FontSize',14);
%%
for i = 1:90-2
    lat_mid(i) = mean(AM4_Data.lat(i:i+2));
end
%%

pressure_range_selection=2;
runs = [1:15];
baseline_temp = squeeze(nanmean(T_mean(8,:,:,2),2));   %take the average T for the base run
baseline_sat  = squeeze(nanmean(sat_time(8,:,:,2),2)); %take the average sat time for the base run
baseline_sat_std  = squeeze(std(sat_time(8,:,:,2),0,[2],'omitnan')); %take the std for the base run

figure(100), clf, hold on
plot(lat_mid,baseline_sat)
%plot(lat_mid,baseline_sat+baseline_sat_std,'r')
%plot(lat_mid,baseline_sat-baseline_sat_std,'r')

baseline_temp = squeeze(nanmean(T_mean(13,:,:,2),2));   %take the average T for the base run
baseline_sat  = squeeze(nanmean(sat_time(13,:,:,2),2)); %take the average sat time for the base run
baseline_sat_std  = squeeze(std(sat_time(13,:,:,2),0,[2],'omitnan')); %take the std for the base run

plot(lat_mid,baseline_sat)
%plot(lat_mid,baseline_sat+baseline_sat_std,'r')
%plot(lat_mid,baseline_sat-baseline_sat_std,'r')

%% 

pressure_range_selection=2;
runs = [1:15];
baseline_temp = squeeze(nanmean(T_mean(:,:,:,2),2));   %take the average T for the base run
baseline_sat  = squeeze(nanmean(sat_time(:,:,:,2),2)); %take the average sat time for the base run
baseline_sat_std  = squeeze(std(sat_time(:,:,:,2),0,[2],'omitnan')); %take the std for the base run

figure(101), clf, hold on
scatter(lat_mid,baseline_sat,'linewidth',2)
set(0,'DefaultAxesColorOrder',hot(25))
legend(AM4_Data.codename{:})
xlabel('Latitude')
ylabel('Saturation Time')
%plot(lat_mid,baseline_sat+baseline_sat_std,'r')
%plot(lat_mid,baseline_sat-baseline_sat_std,'r')

%%

figure(30), clf, hold on
subplot(2,1,1)
for i = 1:88
x=baseline_temp(:,i);
y=baseline_sat(:,i);

[p,S] = polyfit(x,y,1); 
pp(i)=p(1);
[y_fit,delta] = polyval(p,x,S);
plot(x,y,'bo')
hold on
plot(x,y_fit,'r-')
end
[p,S] = polyfit(baseline_temp(:,:),baseline_sat(:,:),1)
[y_fit,delta] = polyval(p,[205:260],S);
plot([205:260],y_fit,'k-','linewidth',2)

xlabel('Temperature'); ylabel('Saturation Time'); set(gca,'fontsize',14)
subplot(2,1,2)
plot(lat_mid,pp)
yline(0); yline(pp(67),'--');yline(p(1),'-.');
legend('slope', 'no trend', 'trend at 45N')
xlabel('Temperature');ylabel('Slope of Saturation time'); set(gca,'fontsize',14)
%% Store all the temperature profiles, similar to the u values above

b=zeros([15, 20, 144,    90,    33-9,   50]);
for run_selection = [1:15];
    File_Numbers=AM4_Data.File_Numbers{run_selection};
    for i=File_Numbers
        AM4_Data_Path = AM4_Data.path{run_selection};
        startLoc = [1 1 1 1]; % Start location along each coordinate
        count  = [Inf Inf 24 50]; % Read until the end of each dimension
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTPS.nc');
        b(run_selection,i,:,:,:,:)=ncread(temp_path,'temp',startLoc,count);
        i
    end
    run_selection
    AM4_Data.codename{run_selection}
end

save(['/scratch/users/mborrus/AM4/TcompleteALL.mat'],'b', '-v7.3');

%%

load('/scratch/users/mborrus/AM4/TcompleteALL.mat')
% size = [21     5   144    90    24   100]
% size = [Run   Int  Lon    Lat   p    day]
T_mean=nan(15,13,88,3);
for run_selection = [1:15]
    for lat_iteration = 1:88
        for pres=1:3
            range_temp = lat_iteration:2+lat_iteration; 
            p_range_temp = AM4_Data.p_range{pres};
            %
            T_mean(run_selection,AM4_Data.File_Numbers{run_selection},lat_iteration,pres) = squeeze(nanmean(b(run_selection,AM4_Data.File_Numbers{run_selection},:,range_temp,AM4_Data.p_range{pres},:),[3 4 5 6]));

        end
    end
    disp("run_selection")
    %RMS_storage.T_mean{run_selection}=
end

Tmatrix=T_mean(:,1:9,45,2)-nanmean(T_mean(8,1:13,45,2),'all');
Tmatrix(Tmatrix<-200)=NaN
round(Tmatrix-nanmean(Tmatrix,2),1)
save('/scratch/users/mborrus/AM4/T_means.mat','T_mean')

%%
t1=squeeze(nanmean(b(1,1,:,range_temp,AM4_Data.p_range{pres},:),'all'))
t4=squeeze(nanmean(b(4,1:3,:,range_temp,AM4_Data.p_range{pres},:),'all'))
tbase=squeeze(nanmean(b(8,1:3,:,range_temp,AM4_Data.p_range{pres},:),'all'))
t1-t4
toak=ncread('dailyUTP.nc','temp',startLoc,count);
tttoak=squeeze(nanmean(toak(:,range_temp,AM4_Data.p_range{pres},:),'all'))