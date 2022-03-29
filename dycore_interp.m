for hlevel = 1
    for run = 4:9
        clear u u_interp ps bk sig T T_interp
        filebase = '/oak/stanford/schools/ees/aditis2/dycore_tree/run_IP_tree/workdir_IP_tree_gauss0/ary_';
        fullfile = strcat(filebase,num2str(hlevel),'_',num2str(run),'/atmos_daily.nc')
        savebase = '/scratch/users/mborrus/dycore/spin_up/';
        savefile = strcat(savebase,num2str(hlevel),'0/',num2str(run))
        tic
        % this will take the model output  and interpolate from sigma levels onto
        % pressure levels and saves in a file called u_interp_the number of the
        % file
        %load axis_stuff_64; % load the values for latitude, longitude etc for the T42 hybrid model setup
        ​
        % sigma interpolation values
        b = [0 0.0000189 .0000472 .000102 .000199 .000358 ...
            .000606 .001 .0015 .0023 .0033 .0046 .0064 ...
            .0086 .0115 .0150 .0194 .0248 .0312 .0390 .0483 ...
            .0592 .0720 .0870 .1044 .1244 .1473 .1736 .2035 ...
            .2373 .2755 .3185 .3666 .4205 .4805 .5471 .6209 ...
            .7025 .7925 .8914 1.0];
        reference_press=100000.0000000000;
        a=zeros(1,41);
        a=reference_press*a;
        P_inter1=[0.0212    0.0330    0.0746    0.1505    0.2785    0.4820    0.8030    1.2500    1.9000    2.8000    3.9500    5.5000    7.5000 ...
            10.0500   13.2500   17.2000   22.1000   28.0000   35.1000   43.6500   53.7500   65.6000   79.5000   95.7000  114.4000  135.8500 ...
            160.4500  188.5500  220.4000  256.4000  297.0000  342.5500  393.5500  450.5000  513.8000  584.0000  661.7000  747.5000  841.9500 925.0000];
        P_inter=P_inter1';
        ​
        u=ncread(fullfile,'ucomp');
        'interpolating u'
        [tt jj kk ll]=size(u);
        u_interp=squeeze(zeros(tt,jj,kk,ll));
        ps=ncread(fullfile,'ps');
        bk=ncread(fullfile,'bk');
        sig=diff(bk)/2+bk(1:end-1);
        ​
        tic
        u_interp=zeros(tt,jj,kk,ll);
        parfor t=1:tt
            u_temp = zeros(tt,jj,kk,ll)
            for l=1:ll
                for j = 1:jj
                    ps1=squeeze(ps(t,j,l));
                    Pr=(a+b.*ps1)./100;Prr=diff(Pr)/2+Pr(1:end-1);
                    u_temp(t,j,:,l)=interp1(Prr,squeeze(u(t,j,:,l)),P_inter,'spline');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % this is to Nan data that is "in" the mountain
                    ​
                    if max(Pr)<max(P_inter)
                        mtn=find(max(Pr)<P_inter);
                        u_temp(t,j,mtn,l)=NaN;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            %t
            u_interp(t,:,:,:)=u_temp(t,:,:,:)
        end
        toc
        clear u
        u_interp_01=u_interp;
        save(strcat(savefile,'/u_interp_01.mat'),'u_interp_01');
        clear u_interp_01
        ​clear u_interp
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % interpolate T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        T=ncread(fullfile,'temp');
        [tt jj kk ll]=size(T);
        'interpolating T'
        tic
        T_interp=zeros(tt,jj,kk,ll);
        parfor t=1:tt
            T_temp = zeros(tt,jj,kk,ll)
            for l=1:ll
                for j = 1:jj
                    ps1=squeeze(ps(t,j,l));
                    Pr=(a+b.*ps1)./100;Prr=diff(Pr)/2+Pr(1:end-1);
                    T_temp(t,j,:,l)=interp1(Prr,squeeze(T(t,j,:,l)),P_inter,'spline');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % this is to Nan data that is "in" the mountain
                    ​
                    if max(Pr)<max(P_inter)
                        mtn=find(max(Pr)<P_inter);
                        T_temp(t,j,mtn,l)=NaN;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            %t
            T_interp(t,:,:,:)=T_temp(t,:,:,:)
        end
        toc
        clear T
        T_interp_01=T_interp;
        save(strcat(savefile,'/T_interp_01'),'T_interp_01');
        clear T_interp_01 T_interp
    end
end
