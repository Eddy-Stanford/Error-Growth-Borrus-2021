%% Surface Plot Work
cd '/home/users/mborrus/Matlab_HPC'
load('AM4_Data.mat')
load('/scratch/users/mborrus/AM4/Ucomplete.mat')

% size = [21     5   144    90    24   100]
% size = [Run   Int  Lon    Lat   p    day]

%%
q=1
for run = 2:3
var_diff(:,q,:,:,:,:) = umaster(:,1,:,:,:,:)-umaster(:,run,:,:,:,:);
q=q+1
end
clear umaster
%%
%90 lats = 88 bins
[Ntemps,Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff);

%% Get RMS
var_RMS=NaN(Ntemps,88,3,Ntime);
for lat_iteration = 1:88
    for pres=1:3
        range_temp = lat_iteration:2+lat_iteration; rng = length(range_temp);
        
        var_RMS([1 10:21],lat_iteration,pres,:) = squeeze(nanmean(...
            rms(...
            reshape(...
            var_diff([1 10:21],:,:,range_temp,AM4_Data.p_range{pres},:),...
            length([1 10:21]), Nrun,   rng * Nlon * length(AM4_Data.p_range{pres}),   Ntime)...
            ,3)...
            ,2)...
            );
        
        var_RMS([4:9],lat_iteration,pres,:) = squeeze(nanmean(...
            rms(...
            reshape(...
            var_diff([4:9],1,:,range_temp,AM4_Data.p_range{pres},:),...
            length([4:9]), Nrun-1,   rng * Nlon * length(AM4_Data.p_range{pres}),   Ntime)...
            ,3)...
            ,2)...
            );
    end
end
clear var_diff

%% Get Saturation Times
sat_time=NaN(Ntemps,Nlat-2,3);
a = smoothdata(var_RMS,4,'movmean',10);
b = diff(a,1,4)./a(:,:,:,2:end);    
for i=1:Ntemps
    for j=1:Nlat-2
        for k=1:3
            if find(b(i,j,k,:)<.03,1)
                sat_time(i,j,k)=find(b(i,j,k,:)<.03,1)+1;
            end
        end
    end
end
save('Sat_Times_U.mat','sat_time')
%% Get Temperature means
load('/scratch/users/mborrus/AM4/Tcomplete.mat')
% size = [21     5   144    90    24   100]
% size = [Run   Int  Lon    Lat   p    day]
for lat_iteration = 1:88
    for pres=1:3
        range_temp = lat_iteration:2+lat_iteration; 
        p_range_temp = AM4_Data.p_range{pres};
        %
        Temperature_mean([1 10:21],lat_iteration,pres) = squeeze(nanmean(Tmaster([1 10:21],1:3,:,range_temp,p_range_temp,:),[2 3 4 5 6]));
        Temperature_mean([4:9],lat_iteration,pres) = squeeze(nanmean(Tmaster([4:9],1:2,:,range_temp,p_range_temp,:),[2 3 4 5 6]));

    end
end
save('T_means.mat','Temperature_mean')

%% For use on Local Computer
%% Get Lat Mid-points
for i = 1:90-2
    lat_mid(i) = mean(AM4_Data.lat(i:i+2));
end
%%
load('T_means.mat')
load('Sat_Times_U.mat')
pressure_range_selection=2;
runs = [1 4:18];

Tmean=Temperature_mean(runs,:,pressure_range_selection);
Sats=sat_time(runs,:,pressure_range_selection);
X = lat_mid; %Latitude (deg lat) - 88;
Y = Tmean-Tmean(1,:); %Temperature (K) - 9x88;
Z = Sats; %Saturation Time (days) - 9x88;

%figure(3),clf, scatter(Y,Sats)
% As Temp goes up, Saturation time goes down
% AND
% As abs(Latitude) goes up, Saturation time goes down


minT=round(min(min(Y)),1); %-14.1
maxT=round(max(max(Y)),1); %34.5    

TempAxis=minT:.5:maxT+.25;
count=zeros(88,length(TempAxis));
SatIntGrid=NaN(length(lat_mid),length(TempAxis));

for run = 1:size(Z,1)
    for lat_bin = 1:88
        a = round(Y(run,lat_bin),1);
        index = find(min(abs(a-TempAxis))==abs(a-TempAxis));
        count(lat_bin,index)=count(index)+1;
        SatIntGrid(lat_bin,index)=Z(run,lat_bin);
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

figure(1),clf,hold on
surface(X,TempAxis,SatIntGrid','FaceColor','interp','EdgeColor','interp')
axis([-80 80 -15 35])
colorbar
caxis([17 25])
xline([0 -45 45],'--')
yline([0],'--')
xlabel('latitude'); ylabel('T from base');title('T and lat vs Sat Time U'); set(gca,'FontSize',14);

figure(2),clf,hold on
contourf(X,TempAxis,SatIntGrid',[17 19 21 23 25 27])
colorbar
axis([-80 80 -15 35])
caxis([17 25])
xline([0 -45 45],'--')
yline([0],'--')
xlabel('latitude'); ylabel('T from base');title('T and lat vs Sat Time U'); set(gca,'FontSize',14);
