%% Load Data
addpath '/Users/mborrus/Documents/Stanford/EGR/Final Data and Scripts'
load AM4_Data.mat
load dy_data.mat
addpath '/Users/mborrus/Documents/cbrewer'
%% Figure 1 - h0 vs AM4
for hide=1;
lat_band = 6;
savepath = strcat('/Users/mborrus/Documents/Stanford/EGR/Plots/Scratch/Master/remade/');
AM4_Data.p_color{1} = ["#fdb863"]; AM4_Data.p_color{2}=["#e66101"];
dy_data.p_color{1} = ["#b2abd2"]; dy_data.p_color{2}=["#5e3c99"];
U_duo_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
sp(1) = subplot(2,1,1);
hold on
for i=1:2;
    for j = 1:5;
        a =squeeze(dy_data.U_All_RMS{1}(i,lat_band,j,1:100)); 
        plot_rms(i,j) = plot(1:100, a, 'Color', dy_data.p_color{i},'LineWidth',2);%plot_rms(i,j).Color(4)=.7;
    end
    xline(dy_data.U_Sat_Time{1}(i,lat_band),'-','LineWidth',3,'Color',dy_data.p_color{i})

     i; legname(i) = strcat("Dycore Base ", dy_data.p_range_names{i}) ;
end
for i = 3:4; 
    for j = 1:5;
        a =squeeze(AM4_Data.U_All_RMS{1}(i-2,lat_band,j,1:100)); 
        plot_rms(i,j) = plot(1:100, a, 'Color', AM4_Data.p_color{i-2},'LineWidth',2);%plot_rms(i,j).Color(4)=.7;
    end
    xline(AM4_Data.U_Sat_Time{1}(i-2,lat_band),'-','LineWidth',3,'Color',AM4_Data.p_color{i-2})    
    i; legname(i) = strcat("AM4 Base ", AM4_Data.p_range_names{i-2}) ;
end
sp(1).TitleFontSizeMultiplier = 1.5
uistack(plot_rms(1,:),'top')
uistack(plot_rms(2,:),'top')

axis([1 100 0 30]), ylabel('U RMSE'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')

pbaspect([3 1 1])
sp(1).InnerPosition=sp(1).InnerPosition-[0.05    0.05   -0.05   -0.05];

sp(2) = subplot(2,1,2) ; hold on
sp(2).InnerPosition=sp(2).InnerPosition-[0.05    0.05   -0.05   -0.05];
for i = 1:2; 
    for j = 1:5;
        a =squeeze(dy_data.T_All_RMS{1}(i,lat_band,j,1:100)); 
        plot_rms(i,j) = plot(1:100, a, 'Color', dy_data.p_color{i},'LineWidth',2);%plot_rms(i,j).Color(4)=.7;
    end
    xline(dy_data.T_Sat_Time{1}(i,lat_band),'-','LineWidth',3,'Color',dy_data.p_color{i})
%     i; legname(i) = "Dycore h0 " + dy_data.p_range_names{i} ;
end
for i = 3:4; 
    for j = 1:5;
        a =squeeze(AM4_Data.T_All_RMS{1}(i-2,lat_band,j,1:100)); 
        plot_rms(i,j) = plot(1:100, a, 'Color', AM4_Data.p_color{i-2},'LineWidth',2);%plot_rms(i,j).Color(4)=.7;
    end
    xline(AM4_Data.T_Sat_Time{1}(i-2,lat_band),'-','LineWidth',3,'Color',AM4_Data.p_color{i-2})
%     i; legname(i) = "AM4 Base " + AM4_Data.p_range_names{i-2} ;
end
uistack(plot_rms(1,:),'top')
uistack(plot_rms(2,:),'top')
 hleg = legend(plot_rms(:,1),legname(:),'Location','NorthEast','NumColumns',1,'FontSize',13,'Box','off'); 
axis([1 100 0 15]), xlabel('Days'), ylabel('T RMSE'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([3 1 1])
saveas(U_duo_RMS_Plot,strcat(savepath,'fig_1','.png'))
saveas(U_duo_RMS_Plot,strcat(savepath,'fig_1','.pdf'))     

end

%% Figure 2 - Intra-Dycore Comparison
for hide=2;
lat_band = 6;
savepath = strcat('/Users/mborrus/Documents/Stanford/EGR/Plots/Scratch/Master/remade/');
AM4_Data.p_color{1} = ["#fdb863"]; AM4_Data.p_color{2}=["#e66101"];
dy_data.Color{2} = ["#08519c"]; dy_data.Color{3}=["#6baed6"];
dy_data.Color{4} = ["#9ecae1"];
dycore_Plot = figure(20); clf, hold on, clear a b plot_rms b_sat 
sp(1) = subplot(2,1,1);
hold on
p=1; dasher(1,:) = ['- '];dasher(2,1:2) = ['- '];
for run_selection = 2:8;
    for i=2;
        for j = 1:2;
            a =squeeze(dy_data.U_All_RMS{run_selection}(i,lat_band,j,1:400)); 
            plot_rms(p,j) = plot(1:400, a,dasher(i,1:2), 'Color', dy_data.Color{run_selection},'LineWidth',2);plot_rms(p,j).Color(4)=.7;
        end
        xline(dy_data.U_Sat_Time{run_selection}(2,lat_band),dasher(i,1:2),'LineWidth',4,'Color',dy_data.Color{run_selection})
        i; legname(p) = {"\delta" + dy_data.codename{run_selection} + ": "+ num2str(dy_data.U_Sat_Time{run_selection}(i,lat_band))} ;
        p = p+1;
    end
end
% title(["Dycore Temperature Variations - Upper Troposphere"]) 
hleg = legend(plot_rms(:,1),legname(:),'Location','NorthWest','NumColumns',1,'FontSize',16,'Box','off'); 
axis([1 250 0 30]), ylabel('U RMSE'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([3 1 1])
set(gca,'TitleFontSizeMultiplier',1.5);

sp(1).InnerPosition=sp(1).InnerPosition-[0.05    0.05   -0.05   -0.05];

sp(2) = subplot(2,1,2) ; hold on
sp(2).InnerPosition=sp(2).InnerPosition-[0.05    0.02   -0.05   -0.05];

p=1; dasher(1,:) = ['- '];dasher(2,1:2) = ['- '];
for run_selection = 2:8;
    for i=2;
        for j = 1:2;
            a =squeeze(dy_data.T_All_RMS{run_selection}(i,lat_band,j,1:400)); 
            plot_rms(p,j) = plot(1:400, a,dasher(i,1:2), 'Color', dy_data.Color{run_selection},'LineWidth',2);plot_rms(p,j).Color(4)=.7;
        end
        xline(dy_data.T_Sat_Time{run_selection}(2,lat_band),dasher(i,1:2),'LineWidth',4,'Color',dy_data.Color{run_selection})
        i; legname(p) = {"\delta" + dy_data.codename{run_selection} + ": "+ num2str(dy_data.T_Sat_Time{run_selection}(i,lat_band)) };
        p = p+1;
    end
end
% title(["Dycore Temperature Variations - Upper Troposphere"]) 
hleg = legend(plot_rms(:,1),legname(:),'Location','NorthWest','NumColumns',1,'FontSize',16,'Box','off'); 
axis([1 250 0 30]), ylabel('T RMSE');xlabel('Days'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([3 1 1])
set(gca,'TitleFontSizeMultiplier',1.5);
% set(gca,'Position',[0.06 0.06 .9 .9])

% saveas(dycore_Plot,strcat(savepath,'fig_2','.png'))     
% saveas(dycore_Plot,strcat(savepath,'fig_2','.pdf'))     
%
end

%% Figure 3 - AM4 RMS plot
for hide=3;
run_selection = [17 16 15 11 1 1 1 10 12 13 14];
lat_band = 10;
savepath = strcat('/Users/mborrus/Documents/Stanford/EGR/Plots/Scratch/Master/remade/');
addpath('/Users/mborrus/Documents/Stanford/Eady/cbrewer')
%[colormap] = cbrewer("div","RdYlBu",9); colormap = colormap([1:4 6:9],:);
AM4_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
sp(1) = subplot(2,1,1);

colormap = [0.1451    0.2039    0.5804;  0.1725    0.4980    0.7216; 0.2549    0.7137    0.7686; ...
            0.4980    0.8039    0.7333; 0.4980    0.8039    0.7333;...
            0     0     0 ;...
            0.9922    0.8000    0.5412; 0.9922    0.8000    0.5412; ...
            0.9882    0.5529    0.3490; 0.8431    0.1882    0.1216; 0.7020         0         0];

hold on; plot_num = 1
for run_N = run_selection;
    for i = 3; 
        for j = 1:2;
            a =squeeze(AM4_Data.U_All_RMS{run_N}(i,lat_band,j,1:100)); 
            plot_rms(plot_num,j) = plot(1:100, a,AM4_Data.Line_Style{run_N}, 'Color', colormap(plot_num,:),'LineWidth',2);%plot_rms(plot_num,j).Color(4)=1;
        end
        xlines1(plot_num)=xline(AM4_Data.U_Sat_Time{run_N}(i,lat_band),AM4_Data.Line_Style{run_N},'LineWidth',3,'Color',colormap(plot_num,:)); xlines1(plot_num).Color(4)=1;
%         i; legname(plot_num) = num2str(AM4_Data.codename{run_N} + ": " + AM4_Data.U_Sat_Time{run_N}(i,lat_band));
%         i; legname(plot_num) = {AM4_Data.codename{run_N}};

        plot_num=plot_num+1;
    end
end

legname(6)={"Base"};
sp(1).TitleFontSizeMultiplier = 1.5
axis([1 100 0 40]), ylabel('U RMSE'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
% hleg1 = legend(xlines1,legname(:),'Location','EastOutside','NumColumns',1,'FontSize',12,'Box','off'); 
uistack(plot_rms(6,:),'top')
uistack(xlines1(:))
pbaspect([3 1 1])
sp(1).InnerPosition=sp(1).InnerPosition-[0.05    0.05   -0.05   -0.05];

sp(2) = subplot(2,1,2) ; hold on
sp(2).InnerPosition=sp(2).InnerPosition-[0.05    0.02   -0.05   -0.05];
hold on; plot_num = 1
for run_N = run_selection;
    for i = 3; 
        for j = 1:2;
            a =squeeze(AM4_Data.T_All_RMS{run_N}(i,lat_band,j,1:100)); 
            plot_rms(plot_num,j) = plot(1:100, a,AM4_Data.Line_Style{run_N}, 'Color', colormap(plot_num,:),'LineWidth',2);%plot_rms(plot_num,j).Color(4)=.7;
        end
        xlines2(plot_num)=xline(AM4_Data.T_Sat_Time{run_N}(i,lat_band),AM4_Data.Line_Style{run_N},'LineWidth',3,'Color',colormap(plot_num,:));    
%         i; legname(plot_num) = num2str(AM4_Data.T_Sat_Time{run_N}(i,lat_band)) + ": " + AM4_Data.codename{run_N} ;
%         i; legname(plot_num) = AM4_Data.codename{run_N};

    plot_num=plot_num+1;
    end
end
uistack(plot_rms(6,:),'top')
uistack(xlines2(6),'top')

legname(6)={"Base"};
% hleg2 = legend(xlines2,legname(:),'Location','EastOutside','NumColumns',1,'FontSize',12,'Box','off'); 
axis([1 100 0 20]), xlabel('Days'), ylabel('T RMSE'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([3 1 1])
% saveas(AM4_RMS_Plot,strcat(savepath,'fig_3_100','.png'))  
% saveas(AM4_RMS_Plot,strcat(savepath,'fig_3_100','.pdf'))     
% saveas(AM4_RMS_Plot,strcat(savepath,'test','.png'))  

end

%% Figure 4 - Six Panel N^2 plot
for hide=4;
run_selection = [1 10];
lat_band = 6;
savepath = strcat('/Users/mborrus/Documents/Stanford/EGR/Plots/Scratch/Master/remade/');
addpath('/Users/mborrus/Documents/Stanford/Eady/cbrewer')

[colormap] = cbrewer("div","RdYlBu",9); colormap = colormap([1:4 6:9],:);
N2_Plot = figure(4); clf, hold on

dryN2_base = squeeze(nanmean(AM4_Data.dry_N2{11}(:,:,:,:),[1,4]))'.*10^5;
moistN2_base = squeeze(nanmean(AM4_Data.moist_N2{11}(:,:,1:12,:),[1,4]))'.*10^5;
limits = [-90 90 225 850]; lines = [5 10 15 20 25 30]; minor_lines = [1 1.2 1.4 1.6 .8];

sp(1) = subplot(2,3,1);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, dryN2_base,lines,'ShowText','on');
ylabel('Pressure (hPa)'), title(strcat("Dry : -4 ",char(176),"C Qobs (Cold)")), c.LineWidth = 2;
axis(limits), set(gca, 'YDir','reverse','XTickLabel',[])
xticks([-90 -60 -30 0 30 60 90])

sp(2) = subplot(2,3,4);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, moistN2_base,lines,'ShowText','on');
ylabel('Pressure (hPa)'), title(strcat("Moist : -4 ",char(176),"C Qobs (Cold)")), c.LineWidth = 2;
axis(limits), set(gca, 'YDir','reverse')
xticks([ -60 -30 0 30 60 ])

dryN2_hot = squeeze(nanmean(AM4_Data.dry_N2{10}(:,:,:,:),[1,4]))'.*10^5;
moistN2_hot = squeeze(nanmean(AM4_Data.moist_N2{10}(:,:,1:12,:),[1,4]))'.*10^5;

sp(3) = subplot(2,3,2);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, dryN2_hot,lines,'ShowText','on');
title(["AM4 N^2 Comparison ([1/s^2] \times10^5)",strcat("Dry : +4 ",char(176),"C Qobs (Hot)")]), c.LineWidth = 2;
axis(limits), set(gca, 'YDir','reverse','YTickLabel',[],'XTickLabel',[])

sp(4) = subplot(2,3,5);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, moistN2_hot,lines,'ShowText','on');
title(strcat("Moist : +4 ",char(176),"C Qobs (Hot)")), c.LineWidth = 2;
xticks([ -60 -30 0 30 60 ])
axis(limits), set(gca, 'YDir','reverse','YTickLabel',[]), xlabel('Latitude')

dryN2_ratio = dryN2_hot./dryN2_base;
moistN2_ratio = moistN2_hot./moistN2_base;

sp(5) = subplot(2,3,3);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, dryN2_ratio,minor_lines,'ShowText','on');
title(strcat("Dry Ratio : ", "Hot/Cold")), c.LineWidth = 2;
axis(limits), set(gca, 'YDir','reverse','YTickLabel',[],'XTickLabel',[])
caxis([.6 1.5])
sp(6) = subplot(2,3,6);
[M,c] = contour(AM4_Data.lat,AM4_Data.p, moistN2_ratio,minor_lines,'ShowText','on');
title(strcat("Moist Ratio : Hot/Cold")), c.LineWidth = 2;
axis(limits), set(gca, 'YDir','reverse','YTickLabel',[])
xticks([ -60 -30 0 30 60])
caxis([.6 1.5])
for i = 1:3;
p = 2*i-1;
sp(p).InnerPosition=sp(p).InnerPosition-[0.05    0.05   -0.05   -0.045];
sp(p).FontSize=12; sp(p+1).FontSize=12;
sp(p+1).InnerPosition=sp(p+1).InnerPosition-[0.05    0.03   -0.05  -0.07];
end
for pp = 1:6
    sp(pp)
    sp(pp).FontName="times new roman"
    set(gca,'fontname','times new roman')
    set(gcf,'color','w');
end
saveas(N2_Plot,strcat(savepath,'fig_4','.png'))   
saveas(N2_Plot,strcat(savepath,'fig_4','.pdf'))     

end

%% Figure 5 - EGR by Latitude
lat_band = 6;

EGR_LatSeries = figure(5); clf; hold on; plot_count = 1; legend_names = [];
plot_choices = [17 16 15 3 11 1 2 10 12 13 14];
clear colormap
colormap = [0.1451-.14    0.2039-.14    0.5804-.14;  0.1725-.15    0.4980-.15    0.7216-.15; 0.2549-.15    0.7137-.15    0.7686-.15; ...
            0.4980    0.8039    0.7333; 0.4980    0.8039    0.7333;...
            0    0    0;...
            0.9922    0.8000    0.5412; 0.9922    0.8000    0.5412; ...
            0.9882    0.5529    0.3490; 0.8431    0.1882    0.1216; 0.7020         0         0];
cN = 1;
for plot_num = plot_choices
    temp_plots(plot_count) = plot(AM4_Data.lat, squeeze( nanmean( AM4_Data.dry_egr{plot_num}(:,:,:,1:20), [1 3 4] ) ) *86400 ,'LineWidth',2, 'Color',colormap(cN,:), 'LineStyle', AM4_Data.Line_Style{plot_num} ); legend_names{plot_count} = [AM4_Data.codename{plot_num}];
    temp_plots(plot_count).Color(4) = .7;  plot_count = plot_count +1;
    cN=cN+1;
end
temp_plots(6).LineWidth = 4;
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([2 1 1])
legend_names(6)={"Base"};
xticks([-90 -60 -30 0 30 60 90])
axis([-90 90 0 2.5]); ylabel("EGR (1/Days)"); xlabel("Latitude"); legend(legend_names,'FontSize',14,'Location','BestOutside'); legend('boxoff'); set(gca,'FontSize',16)
%title(["EGR Lat Series - Days 1 through 20"])
saveas(EGR_LatSeries,strcat(savepath,'fig_5_long','.png')) 
saveas(EGR_LatSeries,strcat(savepath,'fig_5_long','.pdf'))
%% Figure 6 - UTP Sat Times
addpath('/Users/mborrus/Documents/cbrewer')
run_selection = [1 10:17];  
run_count = length(run_selection);
if rem(run_count,2)==1
    splitval=(run_count-1)/2;
    [colormap_cold] = cbrewer("seq","Blues",splitval);
    [colormap_hot] = cbrewer("seq","Reds",splitval);
    colormap = [flip(colormap_cold); 0 0 0; colormap_hot];
else
    splitval=(run_count)/2;
    [colormap_cold] = cbrewer("seq","Blues",splitval);
    [colormap_hot] = cbrewer("seq","Reds",splitval);
    colormap = [flip(colormap_hot); 0 0 0; colormap_cold];
end
colormap

lat_band = [6]
p_range = 2;

        UTP_RMS_Sat_Time_Scatter = figure(6); clf; hold on; plot_count = 1; legend_names = [];
        savepath = strcat('/Users/mborrus/Documents/Stanford/EGR/Plots/Scratch/Master/remade/');
        clear x_data y_data
        C = ["#359B73";"#e69f00";"#2271B2"];
        for plot_choices = run_selection;
            x_data(1,plot_count,:)=mean(AM4_Data.T_Mean{plot_choices}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(1,plot_count,:)=AM4_Data.U_Sat_Time{plot_choices}(:,lat_band);
            x_data(2,plot_count,:)=mean(AM4_Data.T_Mean{plot_choices}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(2,plot_count,:)=AM4_Data.T_Sat_Time{plot_choices}(:,lat_band);
            x_data(3,plot_count,:)=mean(AM4_Data.T_Mean{plot_choices}(p_range,lat_band,:)-AM4_Data.T_Mean{1}(p_range,lat_band,:),3);y_data(3,plot_count,:)=AM4_Data.P_Sat_Time{plot_choices}(:,lat_band);
            plot_count=plot_count+1;
        end
       
        x=squeeze(mean(x_data(:,:,:),[1,3]));
        [sorted,I]=sort(x);
        markers= ['s','d','o'];
        for i = 1:3
            c = linspace(1,10,length(x));
            SPlot(i)=scatter(x(I),y_data(i,I,p_range),80,colormap,'filled',markers(i),'LineWidth',1.5)
            fitvals=polyfit(x(I),y_data(i,I,p_range),1);slope(i)=fitvals(1);
            fitline(i)=plot(x(I),polyval(polyfit(x(I),y_data(i,I,p_range),1),x(I)),'Color',C(i));
            [r,P]=corrcoef(x(I),y_data(i,I,p_range)); r2(i)=r(2,1)^2;p(i)=P(2,1); %LegStr(i)=["r^2 = " + num2str(round(r2(i),2)) + "; p = " + num2str(round(p(i),4))];
        end
        
        ylabel("Days to Saturation"); xlabel(["Temperature Difference from Base (" + char(176) + "C)"]); axis([-8 16 10 30]); 
        %hleg = legend([SPlot(:);fitline(:)],["Wind";"Temp";"Precip";LegStr(1);LegStr(2);LegStr(3)]);
        hleg = legend([SPlot(:)],["Wind";"Temperature";"Precipitation"],'fontsize',18);
        %title(["U Temp and Precip Saturation Times","Upper Troposphere at 45" + char(176) + "N"])
        
legend('boxoff'); set(gca,'FontSize',16)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
pbaspect([2 1 1])
%saveas(UTP_RMS_Sat_Time_Scatter,strcat(savepath,'fig_6','.png')) 
%saveas(UTP_RMS_Sat_Time_Scatter,strcat(savepath,'fig_6_colorful_S','.pdf'))
%% Figure 7 - EGR N dudz Sat Times
clear x_data y_data
    colormap = [166,206,227; 31,120,180; 178,223,138; 51,160,44; 20,10,20]./255
run_selections = [17 16 15 11 1 10 12 13 14];
lat_band = [6]

EGR_N_dudz_scatter = figure(7); clf; hold on; plot_count = 1; legend_names = [];
for run_selection = run_selections;
    x_data(plot_count)=mean(AM4_Data.T_Mean{run_selection}(1:2,lat_band,:)-AM4_Data.T_Mean{1}(1:2,lat_band,:),[1 3]);
    dryegrmean(plot_count)=mean(AM4_Data.dry_egr{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 2 3 4])*86400;
    moistegrmean(plot_count)=mean(AM4_Data.moist_egr{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 2 3 4])*86400;
    dryNmean(plot_count)=mean(AM4_Data.dry_N2{run_selection}(:,AM4_Data.lat_range{6},:,1:20).^(1/2),[1 2 3 4])*86400;
    moistNmean(plot_count)=mean(AM4_Data.moist_N2{run_selection}(:,AM4_Data.lat_range{6},:,1:20).^(1/2),[1 2 3 4])*86400;
%     dudzmean(plot_count) = mean(AM4_Data.du_dz{run_selection}(:,AM4_Data.lat_range{6},:,1:20),[1 2 3 4])*86400;
    plot_count=plot_count+1;
end
dudzmean = [  166.3520; 166.9295; 162.7633; 179.3595; 197.9348; 225.3838; 260.6553; 285.3852; 275.2714] 
x_data_norm = (((x_data+273)-(x_data(5)+273))/(273))*100
scatter(x_data_norm,dryegrmean/dryegrmean(5)*100,80,colormap(2,:),'filled')
scatter(x_data_norm,moistegrmean./moistegrmean(5)*100,120,colormap(1,:),'filled')
scatter(x_data_norm,dryNmean./dryNmean(5)*100,80,colormap(4,:),'filled','s')
scatter(x_data_norm,moistNmean./moistNmean(5)*100,120,colormap(3,:),'filled','s')
scatter(x_data_norm,dudzmean./dudzmean(5)*100,80,colormap(5,:),'filled','d')

scatter(x_data_norm,dryegrmean/dryegrmean(5)*100,80,colormap(2,:),'filled')
scatter(x_data_norm,dryNmean./dryNmean(5)*100,80,colormap(4,:),'filled','s')

 axis([-4 6 80 150])
legend('boxoff'); set(gca,'FontSize',16)
legend('Dry EGR','Moist EGR','Dry N','Moist N','du/dz','Location','NorthWest','fontsize',18)
set(gca,'fontname','times new roman')
set(gcf,'color','w');
xlabel('% Temperature Change from Base')
ylabel('% Change from Base')
pbaspect([2 1 1])

