%% SST PLAY TIME
%/scratch/users/mborrus/AM4/restart_tester/INPUT/
lon = ncread('hadisst_sst.data.nc','lon');
lat = ncread('hadisst_sst.data.nc','lat');
time = ncread('hadisst_sst.data.nc','time');

sst = ncread('hadisst_sst.data.nc','sst');
mon = ncread('hadisst_sst.data.nc','mon');
day = ncread('hadisst_sst.data.nc','day');
year = ncread('hadisst_sst.data.nc','year');

%1331 is the year and month I want and month I want
%% What is our current SST looking like
figure(1), clf
plot(lat,mean(sst(:,:,1331),1)-273.15,'LineWidth',2)
title("What Our SST's used to look like")
xlabel('Latitude');ylabel('SST (C)');
axis([-90 90 -3 30])
set(gca,'FontSize',16)
%% What does Aquaplanet look like
Control = 27*(1 - sin(1.5.*deg2rad(lat)).^2); Control(abs(lat)>60) = 0;
Flat = 27*(1 - sin(1.5.*deg2rad(lat)).^4); Flat(abs(lat)>60) = 0;
Qobs = (Control + Flat)./2;
figure(2), clf, hold on
plot(lat,Control,'LineWidth',2)
plot(lat,Flat,'LineWidth',2)
plot(lat,Qobs,'LineWidth',2)
xlabel('Latitude');ylabel('SST (C)');
title('What SSTs look like from Aqua Planet')
legend("Control","Flat","Qobs")
axis([-90 90 -3 30])
set(gca,'FontSize',16)
%% What would a hotter SST look like? +4c
Control_new = 31*(1 - sin(1.5.*deg2rad(lat)).^2); Control_new(abs(lat)>60) = 0;
Flat_new = 31*(1 - sin(1.5.*deg2rad(lat)).^4); Flat_new(abs(lat)>60) = 0;
Qobs_new = (Control_new + Flat_new)./2;
figure(3), clf, hold on
plot(lat,mean(sst(:,:,1331),1)-273.15,'LineWidth',2)
plot(lat,Qobs,'LineWidth',2)
plot(lat,Qobs_new,'LineWidth',2)
title("What Our SST's would look like using Qobs")
xlabel('Latitude');ylabel('SST (C)');
legend("Original","Base Qobs", "+4C Qobs")
axis([-90 90 -3 35])
set(gca,'FontSize',16)
%% Take original and add a 4deg Qobs peak to it 
Control_new = 4*(1 - sin(1.5.*deg2rad(lat)).^2); Control_new(abs(lat)>60) = 0;
Flat_new = 4*(1 - sin(1.5.*deg2rad(lat)).^4); Flat_new(abs(lat)>60) = 0;
Qobs_new = (Control_new + Flat_new)./2;
SST_hot  = mean(sst(:,:,1331),1)'+Qobs_new-273.15;
SST_cold = mean(sst(:,:,1331),1)'-Qobs_new-273.15;
figure(4), clf, hold on
plot(lat,mean(sst(:,:,1331),1)-273.15,'LineWidth',2)
plot(lat,SST_hot,'LineWidth',2)
plot(lat,SST_cold,'LineWidth',2)
title("Qobs additions")
xlabel('Latitude');ylabel('SST (C)');
legend("Original","Original + 4deg Qobs")
axis([-90 90 -3 35])
set(gca,'FontSize',16)
%% Plotting Temperature gradient
figure(5), clf, hold on
plot(lat(1:end-1),diff(SST_hot))
plot(lat(1:end-1),diff(SST_cold))
title("Qobs Temperature gradients")
xlabel('Latitude');ylabel('K/deg');
legend("hot","cold")
axis([-90 90 -inf inf])
set(gca,'FontSize',16)
%% Plotting Temperature Gradient %
figure(6), clf, hold on
OG_grad = diff(mean(sst(:,:,1331),1))';
plot(lat(1:end-1),100*(diff(SST_hot)-OG_grad)./OG_grad)
plot(lat(1:end-1),100*(diff(SST_cold)-OG_grad)./OG_grad)
title("% change in T gradient")
xlabel('Latitude');ylabel('% change from original');
legend("hot","cold")
axis([-90 90 -100 100])
set(gca,'FontSize',16)
%% Saving data
clear SST_hot SST_cold sst
sst = ncread('hadisst_sst.data.nc','sst');
for i = 1:360
    for j = 1:1776
%SST_hot(i,:,j) = sst(i,:,j)+Qobs_new';
%SST_cold(i,:,j) = sst(i,:,j)-Qobs_new';
SST_m8(i,:,j) = sst(i,:,j)+ (-2.*Qobs_new)';
SST_m12(i,:,j) = sst(i,:,j)+(-3.*Qobs_new)';
SST_m16(i,:,j) = sst(i,:,j)+(-4.*Qobs_new)';
    end
end
%%
%fileattrib('hadisst_sst.data.nc','+w');

copyfile(which('hadisst_sst.data.nc'),'p4_SST.nc');
copyfile(which('hadisst_sst.data.nc'),'m4_SST.nc');

fileattrib('p4_SST.nc','+w');
fileattrib('m4_SST.nc','+w');

ncwrite('p4_SST.nc','sst',SST_hot)
ncwrite('m4_SST.nc','sst',SST_cold)
%%

copyfile(which('hadisst_sst.data.nc'),'m8_SST.nc');
copyfile(which('hadisst_sst.data.nc'),'m12_SST.nc');
copyfile(which('hadisst_sst.data.nc'),'m16_SST.nc');

fileattrib('m8_SST.nc','+w');
fileattrib('m12_SST.nc','+w');
fileattrib('m16_SST.nc','+w');

ncwrite('m8_SST.nc','sst',SST_m8)
ncwrite('m12_SST.nc','sst',SST_m12)
ncwrite('m16_SST.nc','sst',SST_m16)

%%
sst_og = ncread('hadisst_sst.data.nc','sst');
sstm4 = ncread('m4_SST.nc','sst');
sstm8 = ncread('m8_SST.nc','sst');
sstm12 = ncread('m12_SST.nc','sst');
sstm16 = ncread('m16_SST.nc','sst');
sstp4 = ncread('p4_SST.nc','sst');
sstp8 = ncread('p8_SST.nc','sst');
sstp12 = ncread('p12_SST.nc','sst');
sstp16 = ncread('p16_SST.nc','sst');
%%
savepath = strcat('./plots/Master/Finals/scratch/'); mkdir(savepath)
colormap = [0.1451-.14    0.2039-.14    0.5804-.14;  0.1725-.15    0.4980-.15    0.7216-.15; 0.2549-.15    0.7137-.15    0.7686-.15; ...
            0.4980    0.8039    0.7333; 0.4980    0.8039    0.7333;...
            0.7804    0.9137    0.7059;...
            0.9922    0.8000    0.5412; 0.9922    0.8000    0.5412; ...
            0.9882    0.5529    0.3490; 0.8431    0.1882    0.1216; 0.7020         0         0];
        
SST_plot = figure(1); clf, hold on
subplot(1,2,1), hold on
plot(lat,mean(sstp16(:,:,1331),1)-273.15,'Color',colormap(end,:),'LineWidth',2)
plot(lat,mean(sstp12(:,:,1331),1)-273.15,'Color',colormap(end-1,:),'LineWidth',2)
plot(lat,mean(sstp8(:,:,1331),1)-273.15,'Color',colormap(end-2,:),'LineWidth',2)
plot(lat,mean(sstp4(:,:,1331),1)-273.15,'Color',colormap(end-3,:),'LineWidth',2)
plot(lat,mean(sst_og(:,:,1331),1)-273.15,'k-','LineWidth',3)
plot(lat,mean(sstm4(:,:,1331),1)-273.15,'Color',colormap(4,:),'LineWidth',2)
plot(lat,mean(sstm8(:,:,1331),1)-273.15,'Color',colormap(3,:),'LineWidth',2)
plot(lat,mean(sstm12(:,:,1331),1)-273.15,'Color',colormap(2,:),'LineWidth',2)
plot(lat,mean(sstm16(:,:,1331),1)-273.15,'Color',colormap(1,:),'LineWidth',2)

set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([1 1 1])
legend('boxoff'); set(gca,'FontSize',16)
title("Qobs SST Comparison")
xticks([-90 -60 -30 0 30 60 90])
xlabel('Latitude');ylabel(['SST (',char(176),'C)']); leg=legend('+ 16','+ 12','+ 8','+ 4','Default','- 4','- 8','- 12','- 16','Location','westoutside');
axis([-90 90 -10 50])
set(gca,'FontSize',16)
subplot(1,2,2), hold on
plot(lat,mean(sst_og(:,:,1331),1)-273.15+4,'Color',colormap(end-3,:),'LineWidth',3)
plot(lat,mean(sst_og(:,:,1331),1)-273.15,'k-','LineWidth',3)
plot(lat,mean(sst_og(:,:,1331),1)-273.15-4,'Color',colormap(4,:),'LineWidth',3)
xlabel('Latitude'); leg=legend('+ 4','Default','- 4','Location','eastoutside');
pbaspect([1 1 1])
set(gca,'fontname','times new roman')
set(gcf,'color','w');
legend('boxoff'); set(gca,'FontSize',16)
axis([-90 90 -10 50])
title("Uniform SST Comparison")
xticks([-90 -60 -30 0 30 60 90])

saveas(SST_plot,strcat(savepath,'sup_fig_1_long','.png')) 
saveas(SST_plot,strcat(savepath,'sup_fig_2','.pdf'))
