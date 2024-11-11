%%
% SOCIAL-ECONOMY 1
tic
clear;
load('H:\global-PV-wind\Data\GADM_country120_xz2.mat')
load('H:\global-PV-wind\Data\costland_PV120.mat');  % yuan/grid


% ENERGY Production 2
load('H:\global-PV-wind\Data\CP_PV_all.mat'); % 21600x10800
% Capacity power of PV KW
Rearth =6371; % km
for i = 180*120:-1:(0*120+1)
    gridarea120(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea=gridarea120 * 4 *ones(1,360*30);% 1/120*1/30 单位：Km2
clear gridarea120
CP_Area_PV = CP_PV_all*1000./(gridarea*10^6); % W/m2
% max(max(a))
% min(min(a))
CP_Area_PV_lat(:,1) = [90-1/240:-1/120:-90+1/240]'; 
CP_Area_PV_lat(:,2) = CP_Area_PV(:,1); 

load('H:\global-PV-wind\Data\P_CN_12to20_year.mat'); % Ph_12to18_day kWh/yr
Ph_12to18_day = P_CN_12to20_year/10^6/365;  % GWh/grid/day
% Electricity by PV as an average of 2012 and 2018 (GWh/grid/day)
% Electricity = Area * Ratio_Panel_Grid * Solar Radiation * Factor_shield *
% 16.19% * Factor_temperature * 80.56%
% load('..\data\wind_kineticenergy_100m_mean.mat'); % E_100m_mean(4800x1950)
% Wind kinetic energy for 100-m eight after accounting for slope 1979-2018 (kWh/yr)
% Density 4D*14D, D=126 meters
% 59.3% of kinetic energy is convertable to electricity
% t=(Ph_12to18_day./CP_PV_120_all); image(t,'cdatamapping','scaled'); % kWh/year/m2

% GEOGRAPHIC DATA 5
load('H:\global-PV-wind\Data\LC_CN_xz.mat'); %
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 open Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown

% Transmission of electricity 5
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_w_pv_county.dat','-mat');
powerunit = powerunit_w;
clear powerunit_w
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_w_pv_county.dat','-mat');
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_all_w_pv_county.dat','-mat'); % 电厂id，电厂范围内所有面积
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat'); % lines
unitid = unitid_w;
unitid_all = unitid_all_w;
clear unitid_w
clear unitid_all_w

% Grid network 7
% xnet=zeros(21600,10800);
% ynet=zeros(21600,10800);
Rearth    =  6371.3;      % km average radium of the earth
for i = 180*120:-1:(0*120+1)
    gridarea1200(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea = gridarea1200 * 4 *ones(1,360*30);% 1/120*1/30 单位：km2
clear gridarea1200

load('H:\global-PV-wind\Data\E_120cy_ori.mat')  % GWh/day/county各国家的各县分配的power demand
load('H:\global-PV-wind\Data\Gasoline_120cy_ori.mat')  % GWh/day/county各国家的各县分配的power demand
load('H:\global-PV-wind\Data\Diesel_120cy_ori.mat')  % GWh/day/county 各国家的各县分配的power demand
E_120=E_120cy+Gasoline_120cy+Diesel_120cy;
clear E_120cy
clear Gasoline_120cy
clear Diesel_120cy
load('H:\global-PV-wind\Data\E_120cy.dat','-mat');
% E_120(1:2373,35) = E_120cy;
E_120(1:2373,35) = E_120cy(1:2373,1);

% Suitability Factor of Land Use 8
idx0=find(LC_CN==11 | LC_CN==17);
Ph_12to18_day(idx0)=0; % excluding water land
clear idx0
suitablity = 0.15*ones(21600,10800);
idx01=find(LC_CN==15);
suitablity(idx01)=0.1; % excluding water land
clear idx01
% clear LC_CN
% save('H:\world\code\数据处理\PV\final0325\suitablity_PV_sameC.dat','suitablity');

Ph_12to18_day=Ph_12to18_day.* suitablity;
CP_PV_all=CP_PV_all.* suitablity;
load('H:\global-PV-wind\Data\Carbon_density12030_xz.mat');  % g C m-2 yr-1
% Spawn, S.A., Sullivan, C.C., Lark, T.J. et al. Harmonized global maps of above and belowground biomass carbon density in the year 2010. Sci Data 7, 112 (2020). https://doi.org/10.1038/s41597-020-0444-4

load('H:\global-PV-wind\Data\carbonsink12030_xz.mat');  % g C m-2 yr-1
land_sink = carbonsink12030;
clear carbonsink12030
% Hengmao Wang-2022-Global Terrestrial Ecosystem Carbon Flux Inferred from TanSat XCO2 Retrievals

costs=zeros(10000000,12);
% SR2=zeros(100000000,2); % 1列是日照时长，2列是格子数
% costunits=zeros(21600,10800); % id of the cost unit
id=0; % id of the cost unit
% For a power unit, we divide the grids into cost units based on the
% distance to the center of the power unit. In this way, we can find the
% cheapest grids for solar energy
linecost = 2.676; % million RMB cost of transmission line per 1 km
% linecost1000kV = 7.001; % million RMB cost of 1000kV transmission line per 1 km
trascapacity = 300; % capacity of transformer MW
trascost = trascapacity * 0.912; % million RMB for transformer
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
topon=[1 2 2 4 4];
% iidd=0;
i = 0;
for coun = 1:1:34
    [mp,np]=find(powerunit(:,5)==coun);
    if ~isempty (mp)
        powerunit1 = powerunit(mp,:);
        [mm1,nn1]=find(GADM_country120==coun);
        
        country=zeros(21600,10800);
        country(sub2ind(size(country), mm1, nn1))=coun;
        country2 = country(min(mm1):max(mm1),min(nn1):max(nn1));
        
        
        unitid_all2=zeros(21600,10800);
        unitid_all2(sub2ind(size(unitid_all2), mm1, nn1))=unitid_all(sub2ind(size(unitid_all), mm1, nn1));
        unitid_all22 = unitid_all2(min(mm1):max(mm1),min(nn1):max(nn1));
        clear unitid_all2
        
        unitid2=zeros(21600,10800);
        unitid2(sub2ind(size(unitid2), mm1, nn1))=unitid(sub2ind(size(unitid), mm1, nn1));
        unitid22 = unitid2(min(mm1):max(mm1),min(nn1):max(nn1));
        clear unitid2
        
        
        for ii=1:size(powerunit1,1)
            %     i=i+1;
            nummmm=powerunit1(ii,4);
            cy=powerunit1(ii,6);
            [m,n] = find(unitid_all22==mp(ii));
            mpp = min(m)-1+min(mm1)-1;
            npp = min(n)-1+min(nn1)-1;
            a1 = [1:1:max(m)-min(m)+1]';
            xnet = repmat(a1,1,max(n)-min(n)+1);
            a1 = [1:1:max(n)-min(n)+1];
            ynet = repmat(a1,max(m)-min(m)+1,1);
            clear a1
            %     xnet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
            %     ynet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
            %     for ii11=1:max(m)-min(m)+1
            %         for jj11=1:max(n)-min(n)+1
            %             xnet(ii11,jj11)=ii11;
            %             ynet(ii11,jj11)=jj11;
            %         end
            %     end
            idx=find(unitid22(min(m):max(m),min(n):max(n))==mp(ii));
            idx_all=find(unitid_all22(min(m):max(m),min(n):max(n))==mp(ii));
            centerx=powerunit(mp(ii),1)-min(mm1)+1-min(m)+1; % lat 1/120
            centery=powerunit(mp(ii),2)-min(nn1)+1-min(n)+1; % lon 1/30
            powercy=E_120(cy,coun)/ topon(nummmm);
            % powerdemand_plant_pv(mp(ii),1)*1000/365; % electricity demand divided by the number of power units in this county GWh/day/county
            isnottransmit=0;
            z=floor(abs(xnet(idx)-centerx)+abs(ynet(idx)-centery)*4)+1;
            z_all=floor(abs(xnet(idx_all)-centerx)+abs(ynet(idx_all)-centery)*4)+1;
            zmax=max(z,[],1); % the longest distance (1/120 degree)
            unitlength=sqrt(gridarea(powerunit(mp(ii),1),powerunit(mp(ii),2))*4)/4; % length of a 1/120 deg km
            cumulpower=0; % cumulative power
            cumulpower2=0; % cumulative power
            transid=1;
            id1=id+1; % the first cost unit for this power unit
            % Initializing costs
            costs(id1,6)=trascost; % cost of substation million RMB
            for j=1:zmax
                %         iidd = iidd+1;
                id=id+1;
                costs(id,1)=mp(ii);% id of the power unit
                %         costs(id,1)=index_mineral_pv(mp(ii));% id of the power unit
                costs(id,2)=transid; % id of the transformer unit
                costs(id,7)=2*pi*linecost*j*unitlength; % cost of line million RMB
                costs(id,8)=j; % distance
                idx2=find(z==j);
                idx222=find(z_all==j);
                if size(idx2,1)>0
                    for k=1:size(idx222,1)
                        x=xnet(idx_all(idx222(k)))+mpp;
                        y=ynet(idx_all(idx222(k)))+npp;
                        costs(id,12)=costs(id,12)+costland_PV120(x,y)/1e6*suitablity(x,y); % cost of purchasing the land million RMB
                    end
                    for k=1:size(idx2,1)
                        x=xnet(idx(idx2(k)))+mpp;
                        y=ynet(idx(idx2(k)))+npp;
                        cumulpower=cumulpower+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                        if cumulpower>trascapacity
                            id=id+1;
                            transid=transid+1;
                            cumulpower=cumulpower-trascapacity;
                            costs(id,1)=mp(ii); % id of the power unit
                            costs(id,6)=trascost; % cost of substation million RMB
                            %                     SR2(id,1) = SR2(id-1,1);
                            %                     SR2(id,2) = SR2(id-1,2);
                        end
                        % connection to the major line
                        if isnottransmit==0 % already connected to the major line
                            cumulpower2=cumulpower2+Ph_12to18_day(x,y); % GWh/grid/day
                            if cumulpower2>powercy & powercy>=0
                                id=id+1;
                                costs(id,9)=lines(numlines+mp(ii),8)*unitlength*linecost; % cost of transmission to the major line million RMB
                                costs(id,10)=lines(numlines+mp(ii),5); % id of major line (1-47)
                                isnottransmit=1;
                                %                         SR2(id,1) = SR2(id-1,1);
                                %                         SR2(id,2) = SR2(id-1,2);
                            end
                        end
                        costs(id,1)=mp(ii); % id of the power unit
                        %                 costs(id,1)=index_mineral_pv(mp(ii)); % id of the power unit
                        costs(id,2)=transid; % id of the transformer unit
                        costs(id,3)=costs(id,3)+Ph_12to18_day(x,y)*0.365; % electricity (GWh/grid/day) -> TWh / year
                        costs(id,8)=j; % distance
                        costs(id,11)=costs(id,11)+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                        % Carbon for vegetation only
                        % For wind energy, there is no land-use change emission or land carbon sink change
                        if LC_CN(x,y)<=10 || LC_CN(x,y)==12 || LC_CN(x,y)==14
%                             costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea(x,y)*1000/6.67*suitablity(x,y); % land-use change emission t C/ha -> ton C
                            costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea(x,y)*100*suitablity(x,y); % land-use change emission t C/ha -> ton C
                            costs(id,5)=costs(id,5)+land_sink(x,y)*gridarea(x,y)*suitablity(x,y); % land carbon sink gC/m2/yr -> ton C/yr
                        end
                        %                 costunits(x,y)=id; % id of cost units
                    end
                end
            end
            %     i
        end
    end
    coun
end

for coun = 35 % 1:192
    xnet=zeros(4800,1950);
    ynet=zeros(4800,1950);
    gridarea2=zeros(4800,1950);
    Rearth    =  6371.3;      % km average radium of the earth
    for i=1:4800
        for j=1:1950
            xnet(i,j)=i;
            ynet(i,j)=j;
            gridarea2(i,j)=abs(Rearth^2*(sin(((55-(i-1)/120+1/240)+1/120)*pi/180)-sin((55-(i-1)/120+1/240)*pi/180))*1/30*pi/180); %km2
        end
    end
    
    [mp,np]=find(powerunit(:,5)==coun);
    powerunit1 = powerunit(mp,:);
    [mm1,nn1]=find(GADM_country120==coun);
    
    country=zeros(21600,10800);
    country(sub2ind(size(country), mm1, nn1))=coun;
    country2 = country(35*120+1:75*120,(180+73)*30+1:(180+138)*30);
    
    
    unitid_all2=zeros(21600,10800);
    unitid_all2(sub2ind(size(unitid_all2), mm1, nn1))=unitid_all(sub2ind(size(unitid_all), mm1, nn1));
    unitid_all22 = unitid_all2(35*120+1:75*120,(180+73)*30+1:(180+138)*30);
    clear unitid_all2
    
    unitid2=zeros(21600,10800);
    unitid2(sub2ind(size(unitid2), mm1, nn1))=unitid(sub2ind(size(unitid), mm1, nn1));
    unitid22 = unitid2(35*120+1:75*120,(180+73)*30+1:(180+138)*30);
    clear unitid2
    for i=1:size(powerunit1,1)
        display(i);
        idx=find(unitid22==mp(i));
        idx_all=find(unitid_all22==mp(i));
        centerx=powerunit(mp(i),1)-35*120; % lat 1/120
        centery=powerunit(mp(i),2)-(180+73)*30; % lon 1/30
        cy=powerunit(mp(i),6);
        powercy=E_120cy(cy,1) / topon(powerunit1(i,4)); % electricity demand divided by the number of power units in this county GWh/day/county
        isnottransmit=0;
        z=floor(abs(xnet(idx)-centerx)+abs(ynet(idx)-centery)*4)+1;
        z_all=floor(abs(xnet(idx_all)-centerx)+abs(ynet(idx_all)-centery)*4)+1;
        zmax=max(z,[],1); % the longest distance (1/120 degree)
        unitlength=sqrt(gridarea2(centerx,centery)*4)/4; % length of a 1/120 deg km
        cumulpower=0; % cumulative power
        cumulpower2=0; % cumulative power
        transid=1;
        id1=id+1; % the first cost unit for this power unit
        % Initializing costs
        costs(id1,6)=trascost; % cost of substation million RMB
        for j=1:zmax
            id=id+1;
            costs(id,1)=mp(i); % id of the power unit
            costs(id,2)=transid; % id of the transformer unit
            costs(id,7)=2*pi*linecost*j*unitlength; % cost of line million RMB
            costs(id,8)=j; % distance
            idx2=find(z==j);
            idx222=find(z_all==j);
            if size(idx2,1)>0
                for k=1:size(idx222,1)
                    x=xnet(idx_all(idx222(k)))+35*120;
                    y=ynet(idx_all(idx222(k)))+(180+73)*30;
                    costs(id,12)=costs(id,12)+costland_PV120(x,y)/1e6*suitablity(x,y); % cost of purchasing the land million RMB
                end
                for k=1:size(idx2,1)
                    x=xnet(idx(idx2(k)))+35*120;
                    y=ynet(idx(idx2(k)))+(180+73)*30;
                    cumulpower=cumulpower+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                    if cumulpower>trascapacity
                        id=id+1;
                        transid=transid+1;
                        cumulpower=cumulpower-trascapacity;
                        %                     costs(id,1)=mp(i);
                        costs(id,6)=trascost; % cost of substation million RMB
                    end
                    % connection to the major line
                    if isnottransmit==0 % already connected to the major line
                        cumulpower2=cumulpower2+Ph_12to18_day(x,y); % GWh/grid/day
                        if cumulpower2>powercy
                            id=id+1;
                            costs(id,9)=lines(numlines+mp(i),8)*unitlength*linecost; % cost of transmission to the major line million RMB
                            costs(id,10)=lines(numlines+mp(i),5); % id of major line (1-47)
                            isnottransmit=1;
                        end
                    end
                    costs(id,1)=mp(i); % id of the power unit
                    costs(id,2)=transid; % id of the transformer unit
                    costs(id,3)=costs(id,3)+Ph_12to18_day(x,y)*0.365; % electricity (GWh/grid/day) -> TWh / year
                    costs(id,8)=j; % distance
                    costs(id,11)=costs(id,11)+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                    % Carbon for vegetation only
                    % For wind energy, there is no land-use change emission or land carbon sink change
                    if LC_CN(x,y)<=10 || LC_CN(x,y)==12 || LC_CN(x,y)==14
                        costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea2(x-35*120,y-(180+73)*30)*100*suitablity(x,y); % land-use change emission t C/ha -> ton C
                        costs(id,5)=costs(id,5)+land_sink(x,y)*gridarea2(x-35*120,y-(180+73)*30)*suitablity(x,y); % land carbon sink gC/m2/yr -> ton C/yr
                    end
                end
            end
        end
    end
end

for coun = 36:1:192
    [mp,np]=find(powerunit(:,5)==coun);
    if ~isempty (mp)
        powerunit1 = powerunit(mp,:);
        [mm1,nn1]=find(GADM_country120==coun);
        
        country=zeros(21600,10800);
        country(sub2ind(size(country), mm1, nn1))=coun;
        country2 = country(min(mm1):max(mm1),min(nn1):max(nn1));
        
        
        unitid_all2=zeros(21600,10800);
        unitid_all2(sub2ind(size(unitid_all2), mm1, nn1))=unitid_all(sub2ind(size(unitid_all), mm1, nn1));
        unitid_all22 = unitid_all2(min(mm1):max(mm1),min(nn1):max(nn1));
        clear unitid_all2
        
        unitid2=zeros(21600,10800);
        unitid2(sub2ind(size(unitid2), mm1, nn1))=unitid(sub2ind(size(unitid), mm1, nn1));
        unitid22 = unitid2(min(mm1):max(mm1),min(nn1):max(nn1));
        clear unitid2
        
        
        for ii=1:size(powerunit1,1)
            %     i=i+1;
            nummmm=powerunit1(ii,4);
            cy=powerunit1(ii,6);
            [m,n] = find(unitid_all22==mp(ii));
            mpp = min(m)-1+min(mm1)-1;
            npp = min(n)-1+min(nn1)-1;
            a1 = [1:1:max(m)-min(m)+1]';
            xnet = repmat(a1,1,max(n)-min(n)+1);
            a1 = [1:1:max(n)-min(n)+1];
            ynet = repmat(a1,max(m)-min(m)+1,1);
            clear a1
            %     xnet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
            %     ynet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
            %     for ii11=1:max(m)-min(m)+1
            %         for jj11=1:max(n)-min(n)+1
            %             xnet(ii11,jj11)=ii11;
            %             ynet(ii11,jj11)=jj11;
            %         end
            %     end
            idx=find(unitid22(min(m):max(m),min(n):max(n))==mp(ii));
            idx_all=find(unitid_all22(min(m):max(m),min(n):max(n))==mp(ii));
            centerx=powerunit(mp(ii),1)-min(mm1)+1-min(m)+1; % lat 1/120
            centery=powerunit(mp(ii),2)-min(nn1)+1-min(n)+1; % lon 1/30
            powercy=E_120(cy,coun)/ topon(nummmm);
            % powerdemand_plant_pv(mp(ii),1)*1000/365; % electricity demand divided by the number of power units in this county GWh/day/county
            isnottransmit=0;
            z=floor(abs(xnet(idx)-centerx)+abs(ynet(idx)-centery)*4)+1;
            z_all=floor(abs(xnet(idx_all)-centerx)+abs(ynet(idx_all)-centery)*4)+1;
            zmax=max(z,[],1); % the longest distance (1/120 degree)
            unitlength=sqrt(gridarea(powerunit(mp(ii),1),powerunit(mp(ii),2))*4)/4; % length of a 1/120 deg km
            cumulpower=0; % cumulative power
            cumulpower2=0; % cumulative power
            transid=1;
            id1=id+1; % the first cost unit for this power unit
            % Initializing costs
            costs(id1,6)=trascost; % cost of substation million RMB
            for j=1:zmax
                %         iidd = iidd+1;
                id=id+1;
                costs(id,1)=mp(ii);% id of the power unit
                %         costs(id,1)=index_mineral_pv(mp(ii));% id of the power unit
                costs(id,2)=transid; % id of the transformer unit
                costs(id,7)=2*pi*linecost*j*unitlength; % cost of line million RMB
                costs(id,8)=j; % distance
                idx2=find(z==j);
                idx222=find(z_all==j);
                if size(idx2,1)>0
                    for k=1:size(idx222,1)
                        x=xnet(idx_all(idx222(k)))+mpp;
                        y=ynet(idx_all(idx222(k)))+npp;
                        costs(id,12)=costs(id,12)+costland_PV120(x,y)/1e6*suitablity(x,y); % cost of purchasing the land million RMB
                    end
                    for k=1:size(idx2,1)
                        x=xnet(idx(idx2(k)))+mpp;
                        y=ynet(idx(idx2(k)))+npp;
                        cumulpower=cumulpower+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                        if cumulpower>trascapacity
                            id=id+1;
                            transid=transid+1;
                            cumulpower=cumulpower-trascapacity;
                            costs(id,1)=mp(ii); % id of the power unit
                            costs(id,6)=trascost; % cost of substation million RMB
                            %                     SR2(id,1) = SR2(id-1,1);
                            %                     SR2(id,2) = SR2(id-1,2);
                        end
                        % connection to the major line
                        if isnottransmit==0 % already connected to the major line
                            cumulpower2=cumulpower2+Ph_12to18_day(x,y); % GWh/grid/day
                            if cumulpower2>powercy & powercy>=0
                                id=id+1;
                                costs(id,9)=lines(numlines+mp(ii),8)*unitlength*linecost; % cost of transmission to the major line million RMB
                                costs(id,10)=lines(numlines+mp(ii),5); % id of major line (1-47)
                                isnottransmit=1;
                                %                         SR2(id,1) = SR2(id-1,1);
                                %                         SR2(id,2) = SR2(id-1,2);
                            end
                        end
                        costs(id,1)=mp(ii); % id of the power unit
                        %                 costs(id,1)=index_mineral_pv(mp(ii)); % id of the power unit
                        costs(id,2)=transid; % id of the transformer unit
                        costs(id,3)=costs(id,3)+Ph_12to18_day(x,y)*0.365; % electricity (GWh/grid/day) -> TWh / year
                        costs(id,8)=j; % distance
                        costs(id,11)=costs(id,11)+CP_PV_all(x,y)/1e3; % capacity (kW) -> MW
                        % Carbon for vegetation only
                        % For wind energy, there is no land-use change emission or land carbon sink change
                        if LC_CN(x,y)<=10 || LC_CN(x,y)==12 || LC_CN(x,y)==14
%                             costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea(x,y)*1000/6.67*suitablity(x,y); % land-use change emission t C/ha -> ton C
                            costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea(x,y)*100*suitablity(x,y); % land-use change emission t C/ha -> ton C
                            costs(id,5)=costs(id,5)+land_sink(x,y)*gridarea(x,y)*suitablity(x,y); % land carbon sink gC/m2/yr -> ton C/yr
                        end
                        %                 costunits(x,y)=id; % id of cost units
                    end
                end
            end
            %     i
        end
    end
    coun
end

t=costs;
costs=t(1:id,:); % delete unused data
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\costunits01_2_all2.dat','costs','-v7.3');
% save('F:\wyj\SR_12to18_day_power.dat','SR2','-v7.3');
% save('F:\wyj\2costunits01.dat','costs','-v7.3');
% save('F:\wyj\2SR_12to18_day_power.dat','SR2','-v7.3');

% clear
% toc



