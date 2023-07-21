tic
clear;
rmb2us=1/6.8967; % RMB to USD2019
lifetime_power=25;
load('H:\Global PV and wind\Data\discount_country.mat')
discount = discount/100;
discount(35) = 0.05;
load('H:\Global PV and wind\Data\fossilfuel_emissionfactor.mat')  % kg CO2/kWh
discount1yr=zeros(192,1);
for i = 1:192
    for t=1:lifetime_power
        discount1yr(i,1)=discount1yr(i,1)+1/(1+discount(i))^(t-1);
    end
end
% PV:
OMratio_majorline_PV=0.01;
OMratio_substation_PV=0.01;
% wind
OMratio_majorline_wind=0.03;
OMratio_substation_wind=0.03;

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat');  % lines_IX
% Column 15: Distance between two stations
% assumed ±800kV direct  current power station
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_IX_100GW_3_2_all2_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % 电厂编号
clear powerunit_IX_PV

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\optpowerunit_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\tranmission_lines_IX_100GW_3_2_all_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_num_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind;

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\optpowerunit_offshorewind_100GW_county_5%.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_IX_offshorewind_100GW_county_5%.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_IX_100GW_county_5%.mat');  % lines_IX
lines_IX(size(optpowerunit_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_num_IX_offshorewind_100GW_county_5%.mat');  %
powerunit_num_IX_offshorewind(:,4) = 1;
clear unitid_lcoe
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % 
clear powerunit_IX_offshorewind

%%
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\REG_plant_pv_all.mat')  % country，REG的ID(1-4)，UHV Station的ID
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\REG_plant_ons_all.mat')  % country，REG的ID，UHV Station的ID
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\REG_plant_off_all.mat')  % country，REG的ID，UHV Station的ID
for i = 1:size(optpowerunit_PV,1)
    optpowerunit_PV(i,41:43) = REG_plant_pv(optpowerunit_PV(i,40),1:3);
end
for i = 1:size(optpowerunit_onshorewind,1)
    optpowerunit_onshorewind(i,41:43) = REG_plant_ons(optpowerunit_onshorewind(i,40),1:3);
end
for i = 1:size(optpowerunit_offshorewind,1)
    optpowerunit_offshorewind(i,41:43) = REG_plant_off(optpowerunit_offshorewind(i,40),1:3);
end

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\Country_ID_withUHV.mat')  % 有UHV Sation所在的国家ID
Country_ID_withUHV2 = zeros(192,1);
Country_ID_withUHV2(floor(Country_ID_withUHV),1)=1;
clear Country_ID_withUHV
load('H:\Global PV and wind\ANS\UHV_Station_country_all.mat')

%% mineral limitation
load('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 
load('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 
load('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 

[is,pos]=ismember(index_mineral_pv_time2,optpowerunit_PV(:,40));
optpowerunit_PV2 = optpowerunit_PV(pos,:);
powerunit_num_IX_PV2 = powerunit_num_IX_PV(pos,:);
lines_IX_PV2 = lines_IX_PV(pos,:);
clear optpowerunit_PV
clear powerunit_num_IX_PV
clear lines_IX_PV
optpowerunit_PV = optpowerunit_PV2;
powerunit_num_IX_PV = powerunit_num_IX_PV2;
lines_IX_PV = lines_IX_PV2;
clear optpowerunit_PV2
clear powerunit_num_IX_PV2
clear index_mineral_pv
clear lines_IX_PV2

[is,pos]=ismember(index_mineral_ons_time2,optpowerunit_onshorewind(:,40));
optpowerunit_onshorewind2 = optpowerunit_onshorewind(pos,:);
powerunit_num_IX_onshorewind2 = powerunit_num_IX_onshorewind(pos,:);
lines_IX_onshorewind2 = lines_IX_onshorewind(pos,:);
clear optpowerunit_onshorewind
clear powerunit_num_IX_onshorewind
clear lines_IX_onshorewind
optpowerunit_onshorewind = optpowerunit_onshorewind2;
powerunit_num_IX_onshorewind = powerunit_num_IX_onshorewind2;
lines_IX_onshorewind = lines_IX_onshorewind2;
clear optpowerunit_onshorewind2
clear powerunit_num_IX_onshorewind2
clear index_mineral_ons
clear lines_IX_onshorewind2

[is,pos]=ismember(index_mineral_off_time2,optpowerunit_offshorewind(:,40));
optpowerunit_offshorewind2 = optpowerunit_offshorewind(pos,:);
off_pro_IX2 = off_pro_IX(pos,:);
powerunit_num_IX_offshorewind2 = powerunit_num_IX_offshorewind(pos,:);
lines_IX_offshorewind2 = lines_IX_offshorewind(pos,:);
clear optpowerunit_offshorewind
clear powerunit_num_IX_offshorewind
clear off_pro_IX
clear lines_IX_offshorewind
optpowerunit_offshorewind = optpowerunit_offshorewind2;
off_pro_IX = off_pro_IX2;
powerunit_num_IX_offshorewind = powerunit_num_IX_offshorewind2;
lines_IX_offshorewind = lines_IX_offshorewind2;
clear optpowerunit_offshorewind2
clear off_pro_IX2
clear powerunit_num_IX_offshorewind2
clear lines_IX_offshorewind2

%%
optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
clear optpowerunit_PV
clear optpowerunit_onshorewind
clear optpowerunit_offshorewind
lines_IX_offshorewind(:,12:18)=0;
lines_IX = [lines_IX_PV;lines_IX_onshorewind;lines_IX_offshorewind];
clear lines_IX_PV
clear lines_IX_onshorewind
clear lines_IX_offshorewind
powerunit_num_IX_offshorewind(:,5:6)=0;
powerunit_num_IX = [powerunit_num_IX_PV;powerunit_num_IX_onshorewind;powerunit_num_IX_offshorewind];
clear powerunit_num_IX_PV
clear powerunit_num_IX_onshorewind
clear powerunit_num_IX_offshorewind
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
powerunit_IX(:,1)=IX;
optpowerunit_IX(:,1:43)=optpowerunit(IX,1:43);
lines_IX_IX(:,1:15)=lines_IX(IX,1:15);
powerunit_num_IX_IX(:,1:4)=powerunit_num_IX(IX,1:4);
clear optpowerunit
clear lines_IX
clear powerunit_num_IX
clear IX
% [B,IX]=sort(optpowerunit_IX(:,20),1);
line_IX_all = [lines(:,1:15);lines_IX_IX];
num_lines = size(lines,1);
clear lines
clear lines_IX_IX

load('H:\Global PV and wind\Data\GADM_pro120_xz.mat')
load('H:\Global PV and wind\Data\pro_CN_reg.mat') 
for j=1:size(line_IX_all,1)
    if line_IX_all(j,11)==35
        line_IX_all(j,12)=GADM_pro120(line_IX_all(j,3),line_IX_all(j,4)); % 终点pro
        [mm1,nn1]=find(pro_CN_reg(:,1)==min(line_IX_all(j,10),3018));
        line_IX_all(j,13)=pro_CN_reg(mm1,2); % 起点region(1-4,1-7)
        [mm1,nn1]=find(pro_CN_reg(:,1)==line_IX_all(j,12));
        line_IX_all(j,14)=pro_CN_reg(mm1,2); % 终点region(1-4,1-7)
    end
end

position_c(:,1) = 90-line_IX_all(:,1)/120+1/240; % conversion of 行 to lat
position_c(:,2) = line_IX_all(:,2)/30-180-1/60; % conversion of 列 to lon
position_c(:,3) = 90-line_IX_all(:,3)/120+1/240; % conversion of 行 to lat
position_c(:,4) = line_IX_all(:,4)/30-180-1/60; % conversion of 列 to lon
for i = 1:77
    if abs(position_c(i,2)-position_c(i,4))<=180
        distance_CN(i,1) = distance(position_c(i,1),position_c(i,2),position_c(i,3),position_c(i,4))/180*pi*6371; % 单位：千米
    end
    if abs(position_c(i,2)-position_c(i,4))>180
        if position_c(i,2)<0
            a=180+position_c(i,2);
            b=180-position_c(i,4);
            distance_CN(i,1) = distance(position_c(i,1),a,position_c(i,3),b)/180*pi*6371; % 单位：千米
        end
        if position_c(i,2)>0
            a=180-position_c(i,2);
            b=180+position_c(i,4);
            distance_CN(i,1) = distance(position_c(i,1),a,position_c(i,3),b)/180*pi*6371; % 单位：千米
        end
    end
    i
end
clear position_c

%%
P2019_others = zeros(192,6);
load('H:\Global PV and wind\Data\powerdemand_monhour2070_ele288_1112_2_max.mat')  % 288*34 % TWh/h 电气化后
powerdemand_monhour2070_ele288 = powerdemand_monhour2070_ele288;
load('H:\Global PV and wind\Data\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
powerdemand_monhour2070_ele288_country = zeros(288,192);
% ind1 = [];
for i = 1:192
    [m,n]=find(ID_pro(:,3)==i);
    ind = unique(ID_pro(m,1)+1);
    %     ind1 = [ind1;ind];
    powerdemand_monhour2070_ele288_country(:,i) = sum(powerdemand_monhour2070_ele288(:,ind),2);
    P2019_others_pro(ind,1) = sum(powerdemand_monhour2070_ele288(:,ind))./sum(sum(powerdemand_monhour2070_ele288(:,ind)))*P2019_others(i,1);
    clear m
    clear n
    i
end
powerdemand_pro2060_ele = sum(powerdemand_monhour2070_ele288)'/288*8760;

save('H:\Global PV and wind\ANS\powerdemand_monhour2070_ele288_country.mat','powerdemand_monhour2070_ele288_country')  %  % TWh/h 

P2019_others_pro(find(isnan(P2019_others_pro)==1))=0;

for i = 1:192
    powerdemand_monhour(:,i) = powerdemand_monhour2070_ele288_country(:,i)-P2019_others(i,1)/8760;
    i
end
powerdemand_monhour(powerdemand_monhour<0)=0;

%
load('H:\Global PV and wind\Data\powerdemand_monhour2070_REG0811_2_1112_others0_max.mat')  % 288*size(distance_UHV_Station,1), TWh/h，去除其他清洁能源需电量（2020年实际发电量）
powerdemand_monhour2070_REG = powerdemand_monhour2070_REG;

%
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
powerdemand_monhour_region = zeros(288,11);
for i = 1:11
    [m,n]=find(region_ID==i & Country_ID_withUHV2~=0);
    powerdemand_monhour_region(:,i) = sum(powerdemand_monhour(:,m),2);
    asfg(:,i) = sum(powerdemand_monhour2070_ele288_country(:,m),2);
    end
sum(powerdemand_monhour_region)'/288*8760
sum(asfg)'/288*8760

[m,n]=find(Country_ID_withUHV2~=0);
powerdemand_monhour_world = sum(powerdemand_monhour(:,m),2); % TWh/h
%% 
load('H:\Global PV and wind\Data\powergenerat_monhour_pv_100GW_2_all2.mat')  % 288*3639 % TWh/h
load('H:\Global PV and wind\Data\powergenerat_monhour_onshorewind_100GW_2_all2.mat')  % 288*3639 % TWh/h
load('H:\Global PV and wind\Data\powergenerat_monhour_offshorewind_county2.mat')  % 288*3639 %  TWh/h

powergenerat_monhour_pv(find(isnan(powergenerat_monhour_pv)==1)) = 0;
powergenerat_monhour_onshorewind(find(isnan(powergenerat_monhour_onshorewind)==1)) = 0;
powergenerat_monhour_offshorewind(find(isnan(powergenerat_monhour_offshorewind)==1)) = 0;

Ph_pv_spring = sum(powergenerat_monhour_pv(24*2+1:24*5,:),2);
Ph_pv_summer = sum(powergenerat_monhour_pv(24*5+1:24*8,:),2);
Ph_pv_atumn = sum(powergenerat_monhour_pv(24*8+1:24*11,:),2);
Ph_pv_winnter = [sum(powergenerat_monhour_pv(24*11+1:24*12,:),2);sum(powergenerat_monhour_pv(24*0+1:24*2,:),2)];

Ph_onshorewind_spring = sum(powergenerat_monhour_onshorewind(24*2+1:24*5,:),2);
Ph_onshorewind_summer = sum(powergenerat_monhour_onshorewind(24*5+1:24*8,:),2);
Ph_onshorewind_atumn = sum(powergenerat_monhour_onshorewind(24*8+1:24*11,:),2);
Ph_onshorewind_winnter = [sum(powergenerat_monhour_onshorewind(24*11+1:24*12,:),2);sum(powergenerat_monhour_onshorewind(24*0+1:24*2,:),2)];

Ph_offshorewind_spring = sum(powergenerat_monhour_offshorewind(24*2+1:24*5,:),2);
Ph_offshorewind_summer = sum(powergenerat_monhour_offshorewind(24*5+1:24*8,:),2);
Ph_offshorewind_atumn = sum(powergenerat_monhour_offshorewind(24*8+1:24*11,:),2);
Ph_offshorewind_winnter = [sum(powergenerat_monhour_offshorewind(24*11+1:24*12,:),2);sum(powergenerat_monhour_offshorewind(24*0+1:24*2,:),2)];

a(2,1)=sum(sum(Ph_pv_spring))/72;
a(2,2)=sum(sum(Ph_pv_summer))/72;
a(2,3)=sum(sum(Ph_pv_atumn))/72;
a(2,4)=sum(sum(Ph_pv_winnter))/72;

a(3,1)=sum(sum(Ph_onshorewind_spring))/72;
a(3,2)=sum(sum(Ph_onshorewind_summer))/72;
a(3,3)=sum(sum(Ph_onshorewind_atumn))/72;
a(3,4)=sum(sum(Ph_onshorewind_winnter))/72;

a(4,1)=sum(sum(Ph_offshorewind_spring))/72;
a(4,2)=sum(sum(Ph_offshorewind_summer))/72;
a(4,3)=sum(sum(Ph_offshorewind_atumn))/72;
a(4,4)=sum(sum(Ph_offshorewind_winnter))/72;

Ph_pv_season(1:24,1) = Ph_pv_spring(1:24,1) + Ph_pv_spring(25:24*2,1) + Ph_pv_spring(49:24*3,1);
Ph_pv_season(25:24*2,1) = Ph_pv_summer(1:24,1) + Ph_pv_summer(25:24*2,1) + Ph_pv_summer(49:24*3,1);
Ph_pv_season(49:24*3,1) = Ph_pv_atumn(1:24,1) + Ph_pv_atumn(25:24*2,1) + Ph_pv_atumn(49:24*3,1);
Ph_pv_season(24*3+1:24*4,1) = Ph_pv_winnter(1:24,1) + Ph_pv_winnter(25:24*2,1) + Ph_pv_winnter(49:24*3,1);
Ph_pv_season = Ph_pv_season/3;

Ph_onshorewind_season(1:24,1) = Ph_onshorewind_spring(1:24,1) + Ph_onshorewind_spring(25:24*2,1) + Ph_onshorewind_spring(49:24*3,1);
Ph_onshorewind_season(25:24*2,1) = Ph_onshorewind_summer(1:24,1) + Ph_onshorewind_summer(25:24*2,1) + Ph_onshorewind_summer(49:24*3,1);
Ph_onshorewind_season(49:24*3,1) = Ph_onshorewind_atumn(1:24,1) + Ph_onshorewind_atumn(25:24*2,1) + Ph_onshorewind_atumn(49:24*3,1);
Ph_onshorewind_season(24*3+1:24*4,1) = Ph_onshorewind_winnter(1:24,1) + Ph_onshorewind_winnter(25:24*2,1) + Ph_onshorewind_winnter(49:24*3,1);
Ph_onshorewind_season = Ph_onshorewind_season/3;

Ph_offshorewind_season(1:24,1) = Ph_offshorewind_spring(1:24,1) + Ph_offshorewind_spring(25:24*2,1) + Ph_offshorewind_spring(49:24*3,1);
Ph_offshorewind_season(25:24*2,1) = Ph_offshorewind_summer(1:24,1) + Ph_offshorewind_summer(25:24*2,1) + Ph_offshorewind_summer(49:24*3,1);
Ph_offshorewind_season(49:24*3,1) = Ph_offshorewind_atumn(1:24,1) + Ph_offshorewind_atumn(25:24*2,1) + Ph_offshorewind_atumn(49:24*3,1);
Ph_offshorewind_season(24*3+1:24*4,1) = Ph_offshorewind_winnter(1:24,1) + Ph_offshorewind_winnter(25:24*2,1) + Ph_offshorewind_winnter(49:24*3,1);
Ph_offshorewind_season = Ph_offshorewind_season/3;

Ph_season = Ph_pv_season + Ph_onshorewind_season + Ph_offshorewind_season; % TWh/h
%%
optpowerunit_IX_IX = optpowerunit_IX;
clear optpowerunit_IX
CP_cumsum_all = cumsum(optpowerunit_IX_IX(:,30))/10^6;
AA(:,1)=CP_cumsum_all;%TW
% loss rate of electricity
ef1 = 0.967; % power loss rate during the period of circle transportation
ef2 = 0.983; % power loss rate during UHV transportation within the country
ef3 = 0.99*0.85; % power loss rate during charge, discharge and storage period
ef4 = 1; % power loss rate during UHV transportation between countries

tic
load('H:\Global PV and wind\Data\region_connect0811.mat')  % 1 start；2 end
region_connect = region_connect0811;
clear region_connect0811
load('H:\Global PV and wind\Data\distance_UHV_Station_all.mat')  % km
reg=11;
linescapacity=zeros(numlines,1); % capacity of electricity transmissions
linescapacity_plant=zeros(numlines,size(distance_UHV_Station,1));
etrans=zeros(reg,reg); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_type=zeros(reg,reg,3); %
etrans_cou=zeros(192,192); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_cou_type=zeros(192,192,3); %
etrans_REG=zeros(size(distance_UHV_Station,1)+1,size(distance_UHV_Station,1)); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_REG_type=zeros(size(distance_UHV_Station,1)+1,size(distance_UHV_Station,1),3); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2

etrans1=zeros(reg,reg); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_type1=zeros(reg,reg,3); %
etrans_cou1=zeros(192,192); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_cou_type1=zeros(192,192,3); %
etrans_REG1=zeros(size(distance_UHV_Station,1)+1,size(distance_UHV_Station,1)); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_REG_type1=zeros(size(distance_UHV_Station,1)+1,size(distance_UHV_Station,1),3); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2

powergenerat_monhour = zeros(288,numpowerunit);

P_REG = zeros(288,size(UHV_Station_country,1)); % 
P_country = zeros(288,192); % 
P_REG1 = zeros(288,numpowerunit);
P_country1 = zeros(288,numpowerunit); %
P_country2 = zeros(288,numpowerunit); %
P_region2 = zeros(288,numpowerunit);
P_region = zeros(288,reg); % 
P_world = zeros(288,numpowerunit); % 

storage_UHV_storage = zeros(288,numpowerunit);
storage_UHV_storage_regional = zeros(288,numpowerunit);
storage_UHV_storage_inter = zeros(288,numpowerunit);

powerall_reg = zeros(reg,1);
e0_cou = zeros(192,1);
id_REG = zeros(size(UHV_Station_country,1),289);
id_cou = zeros(192,289);
nnnn = 0;
for i = 1: 1:numpowerunit
    i2=i;
    type(i,1) = optpowerunit_IX_IX(i2,35);
    cy=line_IX_all(numlines+i2,9); % county of power unit
    country=line_IX_all(numlines+i2,11); % country of power
    regg1 = line_IX_all(numlines+i2,13); % region of power (1-4,1-7)
    regg2 = line_IX_all(numlines+i2,14); % region of the nearest UHV station from power (1-4,1-7)
    subs = line_IX_all(numlines+i2,6); % substation id
    country2=optpowerunit_IX_IX(i2,41); % country of power
    REG = optpowerunit_IX_IX(i2,43); % UHV Station ID
    dom =optpowerunit_IX_IX(i2,42); % 1-4
    region=region_ID(country,1); % region of power
    i3 =  optpowerunit_IX_IX(i2,31);
    e0_cou(country,1)=e0_cou(country,1)+optpowerunit_IX_IX(i2,1); % electricity produced by the region TWh / year
    powerall_reg(region,1)=powerall_reg(region,1)+optpowerunit_IX_IX(i2,1); % electricity produced by the region TWh / year
    if type(i,1)==1
        powergenerat_monhour(:,i) = powergenerat_monhour_pv(:,i3);
    else if type(i,1)==2
            powergenerat_monhour(:,i) = powergenerat_monhour_onshorewind(:,i3);
        else if type(i,1)==3
                powergenerat_monhour(:,i) = powergenerat_monhour_offshorewind(:,i3);
            end
        end
    end
    powergenerat_monhour(:,i) = powergenerat_monhour(:,i).*optpowerunit_IX_IX(i2,1)/(sum(powergenerat_monhour(:,i))/288*8760);
    ASD = powergenerat_monhour(:,i);
    ASD(find(isnan(ASD)==1)) = 0;
    powergenerat_monhour(:,i) = ASD;
    powergenerat_monhour_plant(:,i) = ASD;
    clear ASD
    
    if country2~=-1 && country2~=35
        for zz = 1:288
            if id_REG(REG,zz)==0 && id_REG(REG,289)==0 && id_cou(country,zz)==0  && id_cou(country,289)==0&& sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))% zz时pro省有效电量<zz时pro省需电量
                % 1
                P_pro_1 = P_REG(zz,REG);
                a = min(min([powergenerat_monhour(zz,i)*ef1+P_REG(zz,REG);powerdemand_monhour2070_REG(zz,REG)]));
                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                b = (a-P_pro_1)/ef1;
                b(b<0)=0;
                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % 
                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1;
                if a==powerdemand_monhour2070_REG(zz,REG)
                    id_REG(REG,zz)=1;
                    if unique(id_REG(REG,1:288))==1
                        id_REG(REG,289)=1;
                    end
                    if country~=184
                        [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG,4));
                        if unique(id_REG(UHV_Station_country(m,1),zz))==1
                            id_cou(UHV_Station_country(REG,4),zz)=1;
                            if unique(id_cou(UHV_Station_country(REG,4),1:288))==1
                                id_cou(UHV_Station_country(REG,4),289)=1;
                            end
                        end
                    end
                    if country==184
                        [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                        if unique(id_REG(UHV_Station_country(m,1),zz))==1
                            id_cou(floor(UHV_Station_country(REG,4)),zz)=1;
                            if unique(id_cou(floor(UHV_Station_country(REG,4)),1:288))==1
                                id_cou(floor(UHV_Station_country(REG,4)),289)=1;
                            end
                        end
                    end
                end
                P_REG(zz,REG) = P_REG(zz,REG) + P_ef1;
                P_REG1(zz,i) = P_REG1(zz,i) + P_ef1;
                P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                P_country(zz,country) = P_country(zz,country) + P_ef1;
                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                P_region(zz,region) = P_region(zz,region) + P_ef1;
                P_world(zz,i) = P_world(zz,i) + P_ef1;
                
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
            end
        end
    end
    
    if country2==35
        REG2_2 = [];
        idx_cn = [];
        if regg2~=regg1  
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg2);
            REG2_2(1) = ma;
        end
        if isempty(REG2_2)
            nnna = 0;
        else
            nnna = size(REG2_2,1);
        end
        idx=find(line_IX_all(1:77,5)==line_IX_all(subs,5));
        idx(idx<subs)=[];
        regg3=unique(line_IX_all(idx,14)); % region of endpoint 1-7
        for iia = 1:size(regg3,1)
            if regg3(iia)~=regg1 && regg3(iia)~=regg2
                [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg3(iia));
                REG2_2(nnna+1) = ma;
                [mmma,nnna2]=find(line_IX_all(idx,14)==regg3(iia));
                idx_cn(nnna+1,1:size([subs:idx(mmma)],2)) = [subs:idx(mmma)];
                nnna = nnna+1;
            end
        end
        
        [ma,na]=find(UHV_Station_country(:,4)==35);
        [B2,IX2]=sort(P_REG(zz,ma)'-powerdemand_monhour2070_REG(zz,ma)',1);
        regg4=IX2;
        for iix = 1:size(regg4,1)
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg4(iix));
            if ~isempty(REG2_2)
                [is,pos]=ismember(ma,REG2_2);
                if is==0
                    REG2_2(nnna+1) = ma;
                    for j=1:77
                        if line_IX_all(j,14)==regg4(iix)
                            idx_cn(nnna+1,1) = j;
                            if line_IX_all(j,13)==regg2
                                idx_cn(nnna+1,1) = j;
                            end
                        end
                    end
                    nnna = nnna+1;
                end
            end
            if isempty(REG2_2)
                REG2_2(nnna+1) = ma;
                for j=1:77
                    if line_IX_all(j,14)==regg4(iix)
                        idx_cn(1,1) = j;
                        if line_IX_all(j,13)==regg2
                            idx_cn(1,1) = j;
                        end
                    end
                end
                nnna = nnna+1;
            end
        end
        
        for zz = 1:288
            if id_REG(REG,zz)==0 && id_REG(REG,289)==0 && id_cou(country,zz)==0  && id_cou(country,289)==0&& sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))% zz时pro省有效电量<zz时pro省需电量
                % 1
                P_pro_1 = P_REG(zz,REG);
                a = min(min([powergenerat_monhour(zz,i)*ef1+P_REG(zz,REG);powerdemand_monhour2070_REG(zz,REG)]));
                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                b = (a-P_pro_1)/ef1;
                b(b<0)=0;
                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1;
                if a==powerdemand_monhour2070_REG(zz,REG)
                    id_REG(REG,zz)=1;
                    if unique(id_REG(REG,1:288))==1
                        id_REG(REG,289)=1;
                    end
                    if country~=184
                        [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG,4));
                        if unique(id_REG(UHV_Station_country(m,1),zz))==1
                            id_cou(UHV_Station_country(REG,4),zz)=1;
                            if unique(id_cou(UHV_Station_country(REG,4),1:288))==1
                                id_cou(UHV_Station_country(REG,4),289)=1;
                            end
                        end
                    end
                    if country==184
                        [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                        if unique(id_REG(UHV_Station_country(m,1),zz))==1
                            id_cou(floor(UHV_Station_country(REG,4)),zz)=1;
                            if unique(id_cou(floor(UHV_Station_country(REG,4)),1:288))==1
                                id_cou(floor(UHV_Station_country(REG,4)),289)=1;
                            end
                        end
                    end
                end
                
                
                P_REG(zz,REG) = P_REG(zz,REG) + P_ef1;
                P_REG1(zz,i) = P_REG1(zz,i) + P_ef1;
                P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                P_country(zz,country) = P_country(zz,country) + P_ef1;
                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                P_region(zz,region) = P_region(zz,region) + P_ef1;
                P_world(zz,i) = P_world(zz,i) + P_ef1;
                
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
            end
        end
    end
    
    
    if country2==-1 % 没有UHV station
        [m,n]=find(UHV_Station_country(:,4)==country);
        REG = UHV_Station_country(m,1);
        for zz = 1:288
            %         if P_country(zz,country) < powerdemand_monhour(zz,country)  && sum(sum(P_country(:,country))) < sum(sum(powerdemand_monhour(:,country)))% zz时pro省有效电量<zz时pro省需电量
            if id_cou(country,zz)==0 && id_cou(country,289)==0
                P_pro_1 = P_country(zz,country);
                a = min(min([powergenerat_monhour(zz,i)*ef1+P_country(zz,country);powerdemand_monhour(zz,country)]));
                if a==powerdemand_monhour(zz,country)
                    id_cou(country,zz)=1;
                    if unique(id_cou(country,1:288))==1
                        id_cou(country,289)=1;
                    end
                end
                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                b = (a-P_pro_1)/ef1;
                b(b<0)=0;
                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1;
                
                P_country(zz,country) = P_country(zz,country) + P_ef1;
                P_REG1(zz,i) = P_REG1(zz,i) + P_ef1;
                P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                P_world(zz,i) = P_world(zz,i) + P_ef1;
                
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
            end
            
            %         if P_country(zz,country) >= powerdemand_monhour(zz,country)  && sum(sum(P_country(:,country))) < sum(sum(powerdemand_monhour(:,country)))% zz时pro省有效电量<zz时pro省需电量
            if id_cou(country,zz)==1 && id_cou(country,289)==0
                P_pro_1 = sum(P_country(:,country));
                a = min(min([powergenerat_monhour(zz,i)*ef1*ef3+sum(P_country(:,country));sum(powerdemand_monhour(:,country))]));
                if a==sum(sum(powerdemand_monhour(:,country)))
                    %                 id_cou(country,289)=1;
                    id_cou(country,:)=1;
                end
                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                b = (a-P_pro_1)/ef1/ef3;
                b(b<0)=0;
                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef3;
                
                storage_UHV_storage(zz,i) = storage_UHV_storage(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                storage_UHV_storage_regional(zz,i) = storage_UHV_storage_regional(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                storage_UHV_storage_inter(zz,i) = storage_UHV_storage_inter(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                
                P_country(zz,country) = P_country(zz,country) + P_ef1;
                P_REG1(zz,i) = P_REG1(zz,i) + P_ef1;
                P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                P_world(zz,i) = P_world(zz,i) + P_ef1;
                
                etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef3;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef3;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef3;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef3;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1/ef3;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1/ef3;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
            end
        end
    end
    
    %     i
end
1
% domestic UHV
nnnn = 0;
for i = 1:1:numpowerunit
    if rem(i,1000)==1 && floor(i/1000)~=0
        nnnn = floor(i/1000);
    end
    i2=i;
    type(i,1) = optpowerunit_IX_IX(i2,35);
    cy=line_IX_all(numlines+i2,9); % county of power unit
    country=line_IX_all(numlines+i2,11); % country of power
    regg1 = line_IX_all(numlines+i2,13); % region of power (1-4,1-7)
    regg2 = line_IX_all(numlines+i2,14); % region of the nearest UHV station from power (1-4,1-7)
    subs = line_IX_all(numlines+i2,6); % substation id
    country2=optpowerunit_IX_IX(i2,41); % country of power 
    REG = optpowerunit_IX_IX(i2,43); % UHV Station ID
    dom =optpowerunit_IX_IX(i2,42); % 1-4
    region=region_ID(country,1); % region of power
    i3 =  optpowerunit_IX_IX(i2,31);
    
    if country2~=-1 && country2~=35
        for zz = 1:288
            if (id_REG(REG,zz)+id_REG(REG,289))>0 && id_cou(country,zz)==0  && id_cou(country,289)==0&& sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))
                % 2
                [mmm1,nnn1] = find(UHV_Station_country(:,4)==country2 & UHV_Station_country(:,8)==dom);
                [mmm,nnn] = find(UHV_Station_country(:,4)==country2 & UHV_Station_country(:,8)~=dom);
                if country2==184
                    [mmm,nnn] = find(UHV_Station_country(:,4)==184 & UHV_Station_country(:,8)~=dom);
                    [mmma,nnna] = find(UHV_Station_country(:,4)==184.1);
                    mmm = [mmm;mmma];
                end
                if country2==184.1
                    [mmm,nnn] = find(UHV_Station_country(:,4)==184.1 & UHV_Station_country(:,8)~=dom);
                    [mmma,nnna] = find(UHV_Station_country(:,4)==184);
                    mmm = [mmm;mmma];
                end
                dis = distance_UHV_Station(mmm,mmm1);
                [BBB,IXXX] = sort(dis);
                for ii2 = 1:1:size(IXXX,1)
                    if powergenerat_monhour(zz,i)>0
                        REG2 = mmm(IXXX(ii2));
                        P_pro_1 = P_REG(zz,REG2);
                        dis2 = dis(IXXX(ii2));
                        if dis2<=10000 & (id_REG(REG2,zz)+id_REG(REG2,289))==0
                            powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                            a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_REG(zz,REG2);powerdemand_monhour2070_REG(zz,REG2)]));
                            if a==powerdemand_monhour2070_REG(zz,REG2)
                                id_REG(REG2,zz)=1;
                                if unique(id_REG(REG2,1:288))==1
                                    id_REG(REG2,289)=1;
                                    id_REG(REG2,:)=1;
                                end
                                if country~=184
                                    [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG2,4));
                                    if unique(id_REG(UHV_Station_country(m,1),zz))==1
                                        id_cou(UHV_Station_country(REG2,4),zz)=1;
                                        if unique(id_cou(UHV_Station_country(REG2,4),1:288))==1
                                            id_cou(UHV_Station_country(REG2,4),:)=1;
                                        end
                                    end
                                end
                                if country==184
                                    [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                                    if unique(id_REG(UHV_Station_country(m,1),zz))==1
                                        id_cou(floor(UHV_Station_country(REG2,4)),zz)=1;
                                        if unique(id_cou(floor(UHV_Station_country(REG2,4)),1:288))==1
                                            id_cou(floor(UHV_Station_country(REG2,4)),:)=1;
                                        end
                                    end
                                end
                            end
                            b = (a-P_pro_1)/(ef1*ef2);
                            b(b<0)=0;
                            P_ef1 = b*ef1*ef2;
                            
                            if P_ef1*10^6>=1
                                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2;
                                P_REG(zz,REG2) = P_REG(zz,REG2) + P_ef1;
                                P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                                P_country(zz,country)  = P_country(zz,country)  + P_ef1;
                                P_region(zz,region) = P_region(zz,region) + P_ef1;
                                P_world(zz,i) = P_world(zz,i) + P_ef1;
                                
                                etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2;
                                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2;
                                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2;
                                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2;
                                etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                                etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                                
                                etrans1(region,region) = etrans1(region,region) + P_ef1;
                                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                                etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    if country2==35
        REG2_2 = [];
        idx_cn = [];
        if regg2~=regg1 %  
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg2);
            REG2_2(1) = ma;
        end
        if isempty(REG2_2)
            nnna = 0;
        else
            nnna = size(REG2_2,1);
        end
        idx=find(line_IX_all(1:77,5)==line_IX_all(subs,5)); % 寻找同一major line
        idx(idx<subs)=[];
        regg3=unique(line_IX_all(idx,14)); % region of endpoint 1-7
        for iia = 1:size(regg3,1)
            if regg3(iia)~=regg1 && regg3(iia)~=regg2
                [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg3(iia));
                REG2_2(nnna+1) = ma;
                [mmma,nnna2]=find(line_IX_all(idx,14)==regg3(iia));
                idx_cn(nnna+1,1:size([subs:idx(mmma)],2)) = [subs:idx(mmma)];
                nnna = nnna+1;
            end
        end
        
        [ma,na]=find(UHV_Station_country(:,4)==35);
        [B2,IX2]=sort(P_REG(zz,ma)'-powerdemand_monhour2070_REG(zz,ma)',1);
        regg4=IX2;
        for iix = 1:size(regg4,1)
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg4(iix));
            if ~isempty(REG2_2)
                [is,pos]=ismember(ma,REG2_2);
                % is是与B大小一致的向量，如果在A中为1，不在为0
                % pos是B中元素如果在A中出现，出现的位置。
                if is==0
                    REG2_2(nnna+1) = ma;
                    for j=1:77
                        if line_IX_all(j,14)==regg4(iix)
                            idx_cn(nnna+1,1) = j;
                            if line_IX_all(j,13)==regg2
                                idx_cn(nnna+1,1) = j;
                            end
                        end
                    end
                    nnna = nnna+1;
                end
            end
            if isempty(REG2_2)
                REG2_2(nnna+1) = ma;
                for j=1:77
                    if line_IX_all(j,14)==regg4(iix)
                        idx_cn(1,1) = j;
                        if line_IX_all(j,13)==regg2
                            idx_cn(1,1) = j;
                        end
                    end
                end
                nnna = nnna+1;
            end
        end
        for zz = 1:288
            if (id_REG(REG,zz)+id_REG(REG,289))>0 && id_cou(country,zz)==0  && id_cou(country,289)==0&& sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))
                % 2
                for iib = 1:size(REG2_2,2)
                    REG2 = REG2_2(iib);
                    if size(idx_cn,1)>=iib
                        idxxx = idx_cn(iib,:);
                        idxxx(idxxx==0)=[];
                    end
                    if powergenerat_monhour(zz,i)>0
                        P_pro_1 = P_REG(zz,REG2);
                        powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                        a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_REG(zz,REG2);powerdemand_monhour2070_REG(zz,REG2)]));
                        if a==powerdemand_monhour2070_REG(zz,REG2)
                            id_REG(REG2,zz)=1;
                            if unique(id_REG(REG2,1:288))==1
                                id_REG(REG2,289)=1;
                                id_REG(REG2,:)=1;
                            end
                            if country~=184
                                [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG2,4));
                                if unique(id_REG(UHV_Station_country(m,1),zz))==1
                                    id_cou(UHV_Station_country(REG2,4),zz)=1;
                                    if unique(id_cou(UHV_Station_country(REG2,4),1:288))==1
                                        %                         id_cou(UHV_Station_country(REG2,4),289)=1;
                                        id_cou(UHV_Station_country(REG2,4),:)=1;
                                    end
                                end
                            end
                            if country==184
                                [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                                if unique(id_REG(UHV_Station_country(m,1),zz))==1
                                    id_cou(floor(UHV_Station_country(REG2,4)),zz)=1;
                                    if unique(id_cou(floor(UHV_Station_country(REG2,4)),1:288))==1
                                        %                         id_cou(floor(UHV_Station_country(REG2,4)),289)=1;
                                        id_cou(floor(UHV_Station_country(REG2,4)),:)=1;
                                    end
                                end
                            end
                        end
                        b = (a-P_pro_1)/(ef1*ef2);
                        b(b<0)=0;
                        P_ef1 = b*ef1*ef2;
                        
                        if P_ef1*10^6>=1
                            powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                            P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2;
                            P_REG(zz,REG2) = P_REG(zz,REG2) + P_ef1;
                            P_country1(zz,i) = P_country1(zz,i) + P_ef1;
                            P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                            P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                            P_country(zz,country)  = P_country(zz,country)  + P_ef1;
                            P_region(zz,region) = P_region(zz,region) + P_ef1;
                            P_world(zz,i) = P_world(zz,i) + P_ef1;
                            
                            etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2;
                            etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2;
                            etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2;
                            etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2;
                            etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                            etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                            
                            etrans1(region,region) = etrans1(region,region) + P_ef1;
                            etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                            etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                            etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                            etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                            etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                        end
                    end
                end
            end
        end
    end
    
end
2
% domestic storage
nnnn = 0;
for i = 1: 1:numpowerunit
    if rem(i,1000)==1 && floor(i/1000)~=0
        nnnn = floor(i/1000);
    end
    i2=i;
    type(i,1) = optpowerunit_IX_IX(i2,35);
    cy=line_IX_all(numlines+i2,9); % county of power unit
    country=line_IX_all(numlines+i2,11); % country of power
    regg1 = line_IX_all(numlines+i2,13); % region of power (1-4,1-7)
    regg2 = line_IX_all(numlines+i2,14); % region of the nearest UHV station from power (1-4,1-7)
    subs = line_IX_all(numlines+i2,6); % substation id
    country2=optpowerunit_IX_IX(i2,41); % country of power 
    REG = optpowerunit_IX_IX(i2,43); % UHV Station ID
    dom =optpowerunit_IX_IX(i2,42); % 1-4
    region=region_ID(country,1); % region of power
    i3 =  optpowerunit_IX_IX(i2,31);
    
    if country2~=-1 && country2~=35
        [mmm1,nnn1] = find(UHV_Station_country(:,5)==region & UHV_Station_country(:,4)==country & UHV_Station_country(:,8)==dom);
        [mmm,nnn] = find(UHV_Station_country(:,5)==region & floor(UHV_Station_country(:,4))==country);
        dis = distance_UHV_Station(mmm,mmm1);
        [BBB,IXXX] = sort(dis);
        for zz = 1:288
            if id_cou(country,zz)==1 && id_cou(country,289)==0&& sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))
                %             3
                for ii2 = 1:1:size(IXXX,1)
                    if powergenerat_monhour(zz,i)>0
                        REG2 = mmm(IXXX(ii2));
                        P_pro_1 = sum(P_REG(:,REG2));
                        dis2 = dis(IXXX(ii2));
                        if dis2<=10000 & id_REG(REG2,289)==0
                            powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                            a = min(min([powergenerat_monhour(zz,i)*ef1*ef2*ef3+sum(P_REG(:,REG2));sum(powerdemand_monhour2070_REG(:,REG2))]));
                            if a==sum(powerdemand_monhour2070_REG(:,REG2))
                                %                 id_REG(REG2,289)=1;
                                id_REG(REG2,:)=1;
                                if country~=184
                                    [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG2,4));
                                    if unique(id_REG(UHV_Station_country(m,1),289))==1
                                        %                     id_cou(UHV_Station_country(REG2,4),289)=1;
                                        id_cou(UHV_Station_country(REG2,4),:)=1;
                                    end
                                end
                                if country==184
                                    [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                                    if unique(id_REG(UHV_Station_country(m,1),289))==1
                                        %                     id_cou(floor(UHV_Station_country(REG2,4)),289)=1;
                                        id_cou(floor(UHV_Station_country(REG2,4)),:)=1;
                                    end
                                end
                            end
                            b = (a-P_pro_1)/(ef1*ef2*ef3);
                            b(b<0)=0;
                            P_ef1 = b*ef1*ef2*ef3;
                            
                            if P_ef1*10^6>=1
                                powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2*ef3;
                                
                                storage_UHV_storage(zz,i) = storage_UHV_storage(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                storage_UHV_storage_regional(zz,i) = storage_UHV_storage_regional(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                storage_UHV_storage_inter(zz,i) = storage_UHV_storage_inter(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                P_REG(zz,REG2) = P_REG(zz,REG2) + P_ef1;
                                P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                                P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                                P_country(zz,country) = P_country(zz,country) + P_ef1;
                                P_region(zz,region) = P_region(zz,region) + P_ef1;
                                P_world(zz,i) = P_world(zz,i) + P_ef1;
                                
                                etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2/ef3;
                                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2/ef3;
                                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2/ef3;
                                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2/ef3;
                                etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2/ef3;
                                etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2/ef3;
                                
                                etrans1(region,region) = etrans1(region,region) + P_ef1;
                                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                                etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    if country2==35
        REG2_2 = [];
        idx_cn = [];
        if regg2~=regg1 % 最近的UHV station所在的region
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg2);
            REG2_2(1) = ma;
        end
        if isempty(REG2_2)
            nnna = 0;
        else
            nnna = size(REG2_2,1);
        end
        idx=find(line_IX_all(1:77,5)==line_IX_all(subs,5)); % 寻找同一major line
        idx(idx<subs)=[];
        regg3=unique(line_IX_all(idx,14)); % region of endpoint 1-7
        for iia = 1:size(regg3,1)
            if regg3(iia)~=regg1 && regg3(iia)~=regg2
                [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg3(iia));
                REG2_2(nnna+1) = ma;
                [mmma,nnna2]=find(line_IX_all(idx,14)==regg3(iia));
                idx_cn(nnna+1,1:size([subs:idx(mmma)],2)) = [subs:idx(mmma)];
                nnna = nnna+1;
            end
        end
        
        [ma,na]=find(UHV_Station_country(:,4)==35);
        [B2,IX2]=sort(P_REG(zz,ma)'-powerdemand_monhour2070_REG(zz,ma)',1);
        regg4=IX2;
        for iix = 1:size(regg4,1)
            [ma,na]=find(UHV_Station_country(:,4)==35 & UHV_Station_country(:,8)==regg4(iix));
            if ~isempty(REG2_2)
                [is,pos]=ismember(ma,REG2_2);
                % is是与B大小一致的向量，如果在A中为1，不在为0
                % pos是B中元素如果在A中出现，出现的位置。
                if is==0
                    REG2_2(nnna+1) = ma;
                    for j=1:77
                        if line_IX_all(j,14)==regg4(iix)
                            idx_cn(nnna+1,1) = j;
                            if line_IX_all(j,13)==regg2
                                idx_cn(nnna+1,1) = j;
                            end
                        end
                    end
                    nnna = nnna+1;
                end
            end
            if isempty(REG2_2)
                REG2_2(nnna+1) = ma;
                for j=1:77
                    if line_IX_all(j,14)==regg4(iix)
                        idx_cn(1,1) = j;
                        if line_IX_all(j,13)==regg2
                            idx_cn(1,1) = j;
                        end
                    end
                end
                nnna = nnna+1;
            end
        end
        for zz = 1:288
            %         if roundn(P_country(zz,country),-10)>=roundn(powerdemand_monhour(zz,country),-10) && roundn(sum(sum(P_country(:,country))),-10) < roundn(sum(sum(powerdemand_monhour(:,country))),-10)
            if id_cou(country,zz)==1 && id_cou(country,289)==0 && sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region)))
                for ii2 = 1:1:size(REG2_2,2)
                    if powergenerat_monhour(zz,i)>0
                        REG2 = REG2_2(ii2);
                        if size(idx_cn,1)>=ii2
                            idxxx = idx_cn(ii2,:);
                            idxxx(idxxx==0)=[];
                        end
                        P_pro_1 = sum(P_REG(:,REG2));
                        powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                        a = min(min([powergenerat_monhour(zz,i)*ef1*ef2*ef3+sum(P_REG(:,REG2));sum(powerdemand_monhour2070_REG(:,REG2))]));
                        if a==sum(powerdemand_monhour2070_REG(:,REG2))
                            %                 id_REG(REG2,289)=1;
                            id_REG(REG2,:)=1;
                            if country~=184
                                [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==UHV_Station_country(REG2,4));
                                if unique(id_REG(UHV_Station_country(m,1),289))==1
                                    %                     id_cou(UHV_Station_country(REG2,4),289)=1;
                                    id_cou(UHV_Station_country(REG2,4),:)=1;
                                end
                            end
                            if country==184
                                [m,n]=find(UHV_Station_country(1:size(UHV_Station_country,1),4)==184 | UHV_Station_country(1:size(UHV_Station_country,1),4)==184.1);
                                if unique(id_REG(UHV_Station_country(m,1),289))==1
                                    %                     id_cou(floor(UHV_Station_country(REG2,4)),289)=1;
                                    id_cou(floor(UHV_Station_country(REG2,4)),:)=1;
                                end
                            end
                        end
                        b = (a-P_pro_1)/(ef1*ef2*ef3);
                        b(b<0)=0;
                        P_ef1 = b*ef1*ef2*ef3;
                        
                        if P_ef1*10^6>=1
                            powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                            P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2*ef3;
                            
                            storage_UHV_storage(zz,i) = storage_UHV_storage(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                            storage_UHV_storage_regional(zz,i) = storage_UHV_storage_regional(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                            storage_UHV_storage_inter(zz,i) = storage_UHV_storage_inter(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                            P_REG(zz,REG2) = P_REG(zz,REG2) + P_ef1;
                            P_country2(zz,i) = P_country2(zz,i) + P_ef1;
                            P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                            P_country(zz,country) = P_country(zz,country) + P_ef1;
                            P_region(zz,region) = P_region(zz,region) + P_ef1;
                            P_world(zz,i) = P_world(zz,i) + P_ef1;
                            
                            etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2/ef3;
                            etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2/ef3;
                            etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2/ef3;
                            etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2/ef3;
                            etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2/ef3;
                            etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2/ef3;
                            
                            etrans1(region,region) = etrans1(region,region) + P_ef1;
                            etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                            etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                            etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                            etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                            etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                        end
                    end
                end
            end
        end
        
    end
    %     i
end
3

nnnn = 0;
for i = 1:1:numpowerunit
    if rem(i,1000)==1 && floor(i/1000)~=0
        nnnn = floor(i/1000);
    end
    
    i2=i;
    type(i,1) = optpowerunit_IX_IX(i2,35);
    cy=line_IX_all(numlines+i2,9); % county of power unit
    country=line_IX_all(numlines+i2,11); % country of power
    regg1 = line_IX_all(numlines+i2,13); % region of power (1-4,1-7)
    regg2 = line_IX_all(numlines+i2,14); % region of the nearest UHV station from power (1-4,1-7)
    subs = line_IX_all(numlines+i2,6); % substation id
    country2=optpowerunit_IX_IX(i2,41); % country of power  % 没有UHV的国家是-1
    REG = optpowerunit_IX_IX(i2,43); % UHV Station ID
    dom =optpowerunit_IX_IX(i2,42); % 1-4
    region=region_ID(country,1); % region of power
    i3 =  optpowerunit_IX_IX(i2,31);
    [mkk,nkk]=find(region_ID==region);
    reg_uni = mean(id_cou(mkk,289));

    if country2~=-1 && country2~=35
        for zz = 1:288
            if powergenerat_monhour(zz,i)>0
                if id_cou(country,289)==1 && reg_uni~=1
                    [mmm1,nnn1] = find(UHV_Station_country(:,5)==region & UHV_Station_country(:,4)==country & UHV_Station_country(:,8)==dom);
                    [mmm,nnn] = find(UHV_Station_country(:,5)==region & floor(UHV_Station_country(:,4))~=country);
                    dis = distance_UHV_Station(mmm,mmm1);
                    [BBB,IXXX] = sort(dis);
                    for ii2 = 1:1:size(IXXX,1)
                        if powergenerat_monhour(zz,i)>0
                            REG2 = mmm(IXXX(ii2));
                            cou2 = floor(UHV_Station_country(mmm(IXXX(ii2)),4));
                            P_country_1 = P_country(zz,cou2);
                            P_country_2 = sum(P_country(:,cou2));
                            dis2 = dis(IXXX(ii2));
                            %                     if dis2<=10000000000 && sum(sum(P_country(:,cou2))) < sum(sum(powerdemand_monhour(:,cou2))) && sum(sum(P_country(zz,cou2))) < sum(sum(powerdemand_monhour(zz,cou2)))
                            if dis2<=10000 && id_cou(cou2,289)==0 && id_cou(cou2,zz)==0
                                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                                a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_country(zz,cou2);powerdemand_monhour(zz,cou2)]));
                                if a==powerdemand_monhour(zz,cou2)
                                    id_cou(cou2,zz)=1;
                                    [m,n]=find(floor(UHV_Station_country(1:size(UHV_Station_country,1),4))==cou2);
                                    id_REG(m,zz)=1;
                                    if unique(id_cou(cou2,1:288))==1
                                        id_cou(cou2,289)=1;
                                        id_REG(m,:)=1;
                                    end
                                end
                                b1 = (a-P_country_1)/(ef1*ef2);
                                b1(b1<0)=0;
                                
                                a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_country_2;sum(powerdemand_monhour(:,cou2))]));
                                b2 = (a-P_country_2)/(ef1*ef2);
                                b2(b2<0)=0;
                                
                                b=min([b1;b2]);
                                P_ef1 = b*ef1*ef2;
                                
                                if P_ef1*10^6>=1
                                    powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                    P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2;
                                    
                                    P_country(zz,cou2) = P_country(zz,cou2) + P_ef1;
                                    P_region(zz,region) = P_region(zz,region) + P_ef1;
                                    P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                                    P_world(zz,i) = P_world(zz,i) + P_ef1;
                                    
                                    etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2;
                                    etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2;
                                    etrans_cou(country,cou2) = etrans_cou(country,cou2) + P_ef1/ef1/ef2;
                                    etrans_cou_type(country,cou2,type(i,1)) = etrans_cou_type(country,cou2,type(i,1)) + P_ef1/ef1/ef2;
                                    etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                                    etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                                    
                                    etrans1(region,region) = etrans1(region,region) + P_ef1;
                                    etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                                    etrans_cou1(country,cou2) = etrans_cou1(country,cou2) + P_ef1;
                                    etrans_cou_type1(country,cou2,type(i,1)) = etrans_cou_type1(country,cou2,type(i,1)) + P_ef1;
                                    etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                    etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                                end
                            end
                        end
                    end
                end
                
                
                %         if roundn(sum(sum(P_country(:,country))),-5) >= roundn(sum(sum(powerdemand_monhour(:,country))),-10) && sum(sum(P_region(:,region))) < sum(sum(powerdemand_monhour_region(:,region))) && sum(sum(P_region(zz,region))) >= sum(sum(powerdemand_monhour_region(zz,region)))
                if id_cou(country,289)==1 && reg_uni~=1
                    [mmm1,nnn1] = find(UHV_Station_country(:,5)==region & UHV_Station_country(:,4)==country & UHV_Station_country(:,8)==dom);
                    [mmm,nnn] = find(UHV_Station_country(:,5)==region & floor(UHV_Station_country(:,4))~=country);
                    dis = distance_UHV_Station(mmm,mmm1);
                    [BBB,IXXX] = sort(dis);
                    for ii2 = 1:1:size(IXXX,1)
                        if powergenerat_monhour(zz,i)>0
                            REG2 = mmm(IXXX(ii2));
                            cou2 = floor(UHV_Station_country(mmm(IXXX(ii2)),4));
                            P_country_1 = sum(P_country(:,cou2));
                            dis2 = dis(IXXX(ii2));
                            %                     if dis2<=10000000000 && sum(sum(P_country(:,cou2))) < sum(sum(powerdemand_monhour(:,cou2)))
                            if dis2<=10000 && id_cou(cou2,289)==0 && id_cou(cou2,zz)==1
                                powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                                a = min(min([powergenerat_monhour(zz,i)*ef1*ef2*ef3*ef4+sum(sum(P_country(:,cou2)));sum(sum(powerdemand_monhour(:,cou2)))]));
                                if a==sum(sum(powerdemand_monhour(:,cou2)))
                                    %                 id_cou(cou2,289)=1;
                                    id_cou(cou2,:)=1;
                                    [m,n]=find(floor(UHV_Station_country(1:size(UHV_Station_country,1),4))==cou2);
                                    %                 id_REG(m,289)=1;
                                    id_REG(m,:)=1;
                                end
                                b = (a-P_country_1)/(ef1*ef2*ef3*ef4);
                                b(b<0)=0;
                                P_ef1 = b*ef1*ef2*ef3*ef4;
                                
                                if P_ef1*10^6>=1
                                    powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                    P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2*ef3*ef4;
                                    storage_UHV_storage_regional(zz,i) = storage_UHV_storage_regional(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                    storage_UHV_storage_inter(zz,i) = storage_UHV_storage_inter(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                    
                                    P_country(zz,cou2) = P_country(zz,cou2) + P_ef1;
                                    P_region(zz,region) = P_region(zz,region) + P_ef1;
                                    P_region2(zz,i) = P_region2(zz,i) + P_ef1;
                                    P_world(zz,i) = P_world(zz,i) + P_ef1;
                                    
                                    etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2/ef3/ef4;
                                    etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                    etrans_cou(country,cou2) = etrans_cou(country,cou2) + P_ef1/ef1/ef2/ef3/ef4;
                                    etrans_cou_type(country,cou2,type(i,1)) = etrans_cou_type(country,cou2,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                    etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2/ef3/ef4;
                                    etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                    
                                    etrans1(region,region) = etrans1(region,region) + P_ef1;
                                    etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                                    etrans_cou1(country,cou2) = etrans_cou1(country,cou2) + P_ef1;
                                    etrans_cou_type1(country,cou2,type(i,1)) = etrans_cou_type1(country,cou2,type(i,1)) + P_ef1;
                                    etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                    etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                                end
                            end
                        end
                    end
                end
                
                if reg_uni==1 && sum(sum(P_region)) < sum(sum(powerdemand_monhour_region))
                    [m,n]=find(region_connect(:,1)==region);
                    reggg = unique(region_connect(m,2)); % 终点
                    mmm = [];
%                     if  sum(sum(P_region(:,reggg))) < sum(sum(powerdemand_monhour_region(:,reggg)))
                        for iiii = 1:size(reggg,1)
                            [mmma,nnna] = find(UHV_Station_country(:,5)==reggg(iiii));
                            mmm = [mmm;mmma];
                        end
                        [mmm1,nnn1] = find(UHV_Station_country(:,5)==region & UHV_Station_country(:,4)==country & UHV_Station_country(:,8)==dom);
                        dis = distance_UHV_Station(mmm,mmm1);
                        [BBB,IXXX] = sort(dis);
                        for ii2 = 1:1:size(IXXX,1)
                            if powergenerat_monhour(zz,i)>0
                                REG2 = mmm(IXXX(ii2));
                                cou2 = floor(UHV_Station_country(mmm(IXXX(ii2)),4));
                                reg2 = UHV_Station_country(mmm(IXXX(ii2)),5);
                                P_country_1 = P_country(zz,cou2);
                                P_country_2 = sum(P_country(:,cou2));
                                dis2 = dis(IXXX(ii2));
                                %                     if dis2<=10000000000 && sum(sum(P_country(:,cou2))) < sum(sum(powerdemand_monhour(:,cou2))) && sum(sum(P_region(:,reg2))) < sum(sum(powerdemand_monhour_region(:,reg2)))
                                if dis2<=10000 && id_cou(cou2,289)==0 && id_cou(cou2,zz)==0
                                    powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                                    a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_country(zz,cou2);powerdemand_monhour(zz,cou2)]));
                                    if a==powerdemand_monhour(zz,cou2)
                                        id_cou(cou2,zz)=1;
                                        [m,n]=find(floor(UHV_Station_country(1:size(UHV_Station_country,1),4))==cou2);
                                        id_REG(m,zz)=1;
                                        if unique(id_cou(cou2,1:288))==1
                                            id_cou(cou2,289)=1;
                                            id_REG(m,:)=1;
                                        end
                                    end
                                    b1 = (a-P_country_1)/(ef1*ef2);
                                    b1(b1<0)=0;
                                    
                                    a = min(min([powergenerat_monhour(zz,i)*ef1*ef2+P_country_2;sum(powerdemand_monhour(:,cou2))]));
                                    b2 = (a-P_country_2)/(ef1*ef2);
                                    b2(b2<0)=0;
                                    
                                    b=min([b1;b2]);
                                    P_ef1 = b*ef1*ef2;
                                    
                                    if P_ef1*10^6>=1
                                        powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                        P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2;
                                        P_country(zz,cou2) = P_country(zz,cou2) + P_ef1;
                                        P_region(zz,reg2) = P_region(zz,reg2) + P_ef1;
                                        P_world(zz,i) = P_world(zz,i) + P_ef1;
                                        
                                        etrans(region,reg2) = etrans(region,reg2) + P_ef1/ef1/ef2;
                                        etrans_type(region,reg2,type(i,1)) = etrans_type(region,reg2,type(i,1)) + P_ef1/ef1/ef2;
                                        etrans_cou(country,cou2) = etrans_cou(country,cou2) + P_ef1/ef1/ef2;
                                        etrans_cou_type(country,cou2,type(i,1)) = etrans_cou_type(country,cou2,type(i,1)) + P_ef1/ef1/ef2;
                                        etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                                        etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                                        
                                        etrans1(region,reg2) = etrans1(region,reg2) + P_ef1;
                                        etrans_type1(region,reg2,type(i,1)) = etrans_type1(region,reg2,type(i,1)) + P_ef1;
                                        etrans_cou1(country,cou2) = etrans_cou1(country,cou2) + P_ef1;
                                        etrans_cou_type1(country,cou2,type(i,1)) = etrans_cou_type1(country,cou2,type(i,1)) + P_ef1;
                                        etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                        etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                                    end
                                end
                                
                            end
                        end
%                     end
                end
                
                if reg_uni==1 && sum(sum(P_region)) < sum(sum(powerdemand_monhour_region))
                    [m,n]=find(region_connect(:,1)==region);
                    reggg = unique(region_connect(m,2)); % 终点
                    mmm = [];
%                     if  sum(sum(P_region(:,reggg))) < sum(sum(powerdemand_monhour_region(:,reggg)))
                        for iiii = 1:size(reggg,1)
                            [mmma,nnna] = find(UHV_Station_country(:,5)==reggg(iiii));
                            mmm = [mmm;mmma];
                        end
                        [mmm1,nnn1] = find(UHV_Station_country(:,5)==region & UHV_Station_country(:,4)==country & UHV_Station_country(:,8)==dom);
                        dis = distance_UHV_Station(mmm,mmm1);
                        [BBB,IXXX] = sort(dis);
                        for ii2 = 1:1:size(IXXX,1)
                            if powergenerat_monhour(zz,i)>0
                                REG2 = mmm(IXXX(ii2));
                                cou2 = floor(UHV_Station_country(mmm(IXXX(ii2)),4));
                                reg2 = UHV_Station_country(mmm(IXXX(ii2)),5);
                                P_country_1 = sum(P_country(:,cou2));
                                dis2 = dis(IXXX(ii2));
                                %                     if dis2<=10000000000 && sum(sum(P_country(:,cou2))) < sum(sum(powerdemand_monhour(:,cou2))) && sum(sum(P_region(:,reg2))) < sum(sum(powerdemand_monhour_region(:,reg2)))
                                if dis2<=10000 && id_cou(cou2,289)==0 && id_cou(cou2,zz)==1
                                    powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                                    a = min(min([powergenerat_monhour(zz,i)*ef1*ef2*ef3*ef4+P_country_1;sum(sum(powerdemand_monhour(:,cou2)))]));
                                    if a==sum(sum(powerdemand_monhour(:,cou2)))
                                        %                 id_cou(cou2,289)=1;
                                        %                 [m,n]=find(floor(UHV_Station_country(1:size(UHV_Station_country,1),4))==cou2);
                                        %                 id_REG(m,289)=1;
                                        id_cou(cou2,:)=1;
                                        [m,n]=find(floor(UHV_Station_country(1:size(UHV_Station_country,1),4))==cou2);
                                        id_REG(m,:)=1;
                                    end
                                    powergenerat_monhour_zz1 = powergenerat_monhour(zz,i);
                                    b = (a-P_country_1)/(ef1*ef2*ef3*ef4); %
                                    b(b<0)=0;
                                    P_ef1 = b*ef1*ef2*ef3*ef4;
                                    
                                    if P_ef1*10^6>=1
                                        powergenerat_monhour(zz,i) = powergenerat_monhour(zz,i)-b; % i电厂剩余可用电量
                                        P_ef1 = (powergenerat_monhour_zz1-powergenerat_monhour(zz,i))*ef1*ef2*ef3*ef4;
                                        storage_UHV_storage_inter(zz,i) = storage_UHV_storage_inter(zz,i) + powergenerat_monhour_zz1-powergenerat_monhour(zz,i);
                                        
                                        P_country(zz,cou2) = P_country(zz,cou2) + P_ef1;
                                        P_region(zz,reg2) = P_region(zz,reg2) + P_ef1;
                                        P_world(zz,i) = P_world(zz,i) + P_ef1;
                                        
                                        etrans(region,reg2) = etrans(region,reg2) + P_ef1/ef1/ef2/ef3/ef4;
                                        etrans_type(region,reg2,type(i,1)) = etrans_type(region,reg2,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                        etrans_cou(country,cou2) = etrans_cou(country,cou2) + P_ef1/ef1/ef2/ef3/ef4;
                                        etrans_cou_type(country,cou2,type(i,1)) = etrans_cou_type(country,cou2,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                        etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2/ef3/ef4;
                                        etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2/ef3/ef4;
                                        
                                        etrans1(region,reg2) = etrans1(region,reg2) + P_ef1;
                                        etrans_type1(region,reg2,type(i,1)) = etrans_type1(region,reg2,type(i,1)) + P_ef1;
                                        etrans_cou1(country,cou2) = etrans_cou1(country,cou2) + P_ef1;
                                        etrans_cou_type1(country,cou2,type(i,1)) = etrans_cou_type1(country,cou2,type(i,1)) + P_ef1;
                                        etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                        etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                                    end
                                end
                            end
                        end
%                     end
                end
            end
        end
    end
    
    
    utilize_ratio2060(i,1) = sum(P_REG1(:,i))./sum(powergenerat_monhour_plant(:,i));
    utilize_ratio2060_UHV(i,1) = sum(P_country1(:,i))./sum(powergenerat_monhour_plant(:,i));
    utilize_ratio2060_UHV_storage(i,1) = sum(P_country2(:,i))./sum(powergenerat_monhour_plant(:,i));
    utilize_ratio2060_UHV_storage_reg(i,1) = sum(P_region2(:,i))./sum(powergenerat_monhour_plant(:,i));
    utilize_ratio2060_UHV_storage_inter(i,1) = sum(P_world(:,i))./sum(powergenerat_monhour_plant(:,i));
    i
end
toc


storage_max_plant_UHV_storage = max(storage_UHV_storage)'; % TWh/h
storage_year_plant_UHV_storage = sum(storage_UHV_storage)'/288*8760; % TWh/year

storage_max_plant_UHV_storage_reg = max(storage_UHV_storage_regional)'; % TWh/h
storage_year_plant_UHV_storage_reg = sum(storage_UHV_storage_regional)'/288*8760; % TWh/year

storage_max_plant_UHV_storage_inter = max(storage_UHV_storage_inter)'; % TWh/h
storage_inter_year_plant_UHV_storage = sum(storage_UHV_storage_inter)'/288*8760; % TWh/year

save('H:\Global PV and wind\ANS\utilize_ratio2060_county_all_pro.mat','utilize_ratio2060')
save('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_county_all_pro.mat','utilize_ratio2060_UHV')
save('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_storage_county_all_pro.mat','utilize_ratio2060_UHV_storage')
save('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_storage_reg_county_all_pro.mat','utilize_ratio2060_UHV_storage_reg')
save('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_pro.mat','utilize_ratio2060_UHV_storage_inter')

save('H:\Global PV and wind\ANS\storage_max_plant_UHV_storage_county_all_pro.mat','storage_max_plant_UHV_storage')   % TWh/h
save('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_county_all_pro.mat','storage_year_plant_UHV_storage')  % TWh/year
save('H:\Global PV and wind\ANS\storage_max_plant_UHV_storage_reg_county_all_pro.mat','storage_max_plant_UHV_storage_reg')   % TWh/h
save('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_reg_county_all_pro.mat','storage_year_plant_UHV_storage_reg')  % TWh/year
save('H:\Global PV and wind\ANS\storage_max_plant_UHV_storage_inter_county_all_pro.mat','storage_max_plant_UHV_storage_inter')   % TWh/h
save('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_inter_county_all_pro.mat','storage_inter_year_plant_UHV_storage')  % TWh/year
