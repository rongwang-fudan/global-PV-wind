tic
clear;
rmb2us=1/6.8967; % RMB to USD2019
lifetime_power=25;
load('H:\global-PV-wind\Data\discount_country.mat')
discount = discount/100;
discount(35)=0.05;
load('H:\global-PV-wind\Data\fossilfuel_emissionfactor.mat')  %kg CO2/kWh
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

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat');  % lines
% Column 15: Distance between two stations
% assumed ±800kV direct  current power station
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_IX_100GW_3_2_all2_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % 电厂编号
clear powerunit_IX_PV

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\optpowerunit_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\tranmission_lines_IX_100GW_3_2_all_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_num_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind;

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\optpowerunit_offshorewind_100GW_county_5%.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_IX_offshorewind_100GW_county_5%.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_IX_100GW_county_5%.mat');  % lines_IX
lines_IX(size(optpowerunit_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_num_IX_offshorewind_100GW_county_5%.mat');  %
powerunit_num_IX_offshorewind(:,4) = 1;
% load('H:\global-PV-wind\ANS\数据处理\unitid_lcoe_offshorewind1227_2_xz_100GW.mat');
% unitid_lcoe_offshorewind = unitid_lcoe;
clear unitid_lcoe
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %


% unitid_lcoe_offshorewind = unitid_lcoe;
optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % 电厂编号
clear powerunit_IX_offshorewind

%%
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\REG_plant_pv_all.mat')  % 各电厂所在的国家，REG的ID(1-4)，UHV Station的ID
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\REG_plant_ons_all.mat')  % 各电厂所在的国家，REG的ID，UHV Station的ID
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\REG_plant_off_all.mat')  % 各电厂所在的国家，REG的ID，UHV Station的ID
for i = 1:size(optpowerunit_PV,1)
    optpowerunit_PV(i,41:43) = REG_plant_pv(optpowerunit_PV(i,40),1:3);
end
for i = 1:size(optpowerunit_onshorewind,1)
    optpowerunit_onshorewind(i,41:43) = REG_plant_ons(optpowerunit_onshorewind(i,40),1:3);
end
for i = 1:size(optpowerunit_offshorewind,1)
    optpowerunit_offshorewind(i,41:43) = REG_plant_off(optpowerunit_offshorewind(i,40),1:3);
end
% 各电厂所在的国家，REG的ID(1-4)，如果该国家没有UHV staton，则赋值为-1

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\Country_ID_withUHV.mat')  % 有UHV Sation所在的国家ID
Country_ID_withUHV2 = zeros(192,1);
Country_ID_withUHV2(floor(Country_ID_withUHV),1)=1;
clear Country_ID_withUHV
load('H:\global-PV-wind\ANS\UHV_Station_country_all.mat')
% 1Substation ID; 2行；3列；4国家ID;
% 5region ID; 6pro ID(0-3638); 7该序号所分配的power demand (2060年电气化后全部需电量，TWh/year); 8REG(1-4)

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\global-PV-wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\global-PV-wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\global-PV-wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

[is,pos]=ismember(index_mineral_pv_time2,optpowerunit_PV(:,40));
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
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
% lines_IX = [lines_IX_PV(:,1:15);lines_IX_onshorewind;lines_IX_offshorewind];
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
% for i=1:numpowerunit
%     i2=IX(i);
%     powerunit_IX(i,1)=i2;
%     optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); %  optpowerunit_IX(i2,1); electricity used by the county TWh / year
%     lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
%     powerunit_num_IX_IX(i,1:4)=powerunit_num_IX(i2,1:4); % lat lon
% end
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
clear lines
clear lines_IX_IX


load('H:\global-PV-wind\Data\GADM_pro120_xz.mat')
load('H:\global-PV-wind\Data\pro_CN_reg.mat') % 第一列是pro ID，第二列是对应中国国内的region ID (1-7)
for j=1:size(line_IX_all,1)
    if line_IX_all(j,11)==35
        line_IX_all(j,12)=GADM_pro120(line_IX_all(j,3),line_IX_all(j,4)); % 终点pro
        [mm1,nn1]=find(pro_CN_reg(:,1)==min(3018,line_IX_all(j,10)));
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
CHO_MAC = 5*ones(size(optpowerunit_IX,1),1);

%%
% load('H:\global-PV-wind\ANS\P2070_others.mat'); % TWh/year
% P2019_others = sum(P2070_others,2);
P2019_others = zeros(192,6);
load('H:\global-PV-wind\Data\powerdemand_monhour2070_ele288_1112_2_max.mat')  % 288*34 % TWh/h 电气化后
load('H:\global-PV-wind\Data\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
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

P2019_others_pro(find(isnan(P2019_others_pro)==1))=0;

for i = 1:192
    powerdemand_monhour(:,i) = powerdemand_monhour2070_ele288_country(:,i)-P2019_others(i,1)/8760;
    i
end
powerdemand_monhour(powerdemand_monhour<0)=0;

%
load('H:\global-PV-wind\Data\powerdemand_monhour2070_REG0811_2_1112_others0_max.mat')  % 288*size(distance_UHV_Station,1), TWh/h，去除其他清洁能源需电量（2020年实际发电量）

for i =1 :1:192
    [m,n]=find(floor(UHV_Station_country(1:583,4))==i);
    if ~isempty(m)
        pp(:,i) = sum(powerdemand_monhour2070_REG(:,m),2);
        pp3(:,i) = powerdemand_monhour(:,i) ;
        iii3(:,i) = powerdemand_monhour(:,i) ;
        
    end
end



%
load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
powerdemand_monhour_region = zeros(288,11);
for i = 1:11
    [m,n]=find(region_ID==i & Country_ID_withUHV2~=0);
    powerdemand_monhour_region(:,i) = sum(powerdemand_monhour(:,m),2);
    [m,n]=find(region_ID==i);
    powerdemand_monhour_region2(:,i) = sum(powerdemand_monhour(:,m),2);
    powerdemand_others_region2(1,i) = sum(P2019_others(m,1));
    powerdemand_all_region2(:,i) = sum(powerdemand_monhour2070_ele288_country(:,m),2);
end

sum(powerdemand_monhour_region)'/288*8760
sum(powerdemand_monhour_region2)'/288*8760
sum(powerdemand_all_region2)/288*8760

[m,n]=find(Country_ID_withUHV2~=0);
powerdemand_monhour_world = sum(powerdemand_monhour(:,m),2); % TWh/h
%% 以UTC时间计算
% 单位需要验证
load('H:\global-PV-wind\Data\powergenerat_monhour_pv_100GW_2_all2.mat')  % 288*3639 % TWh/h
load('H:\global-PV-wind\Data\powergenerat_monhour_onshorewind_100GW_2_all2.mat')  % 288*3639 % TWh/h
load('H:\global-PV-wind\Data\powergenerat_monhour_offshorewind_county2.mat')  % 288*3639 %  TWh/h

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
CF_IX = optpowerunit_IX_IX(:,1)./(optpowerunit_IX_IX(:,30)/10^6*8760);
CP_cumsum_all = cumsum(optpowerunit_IX_IX(:,30))/10^6;
AA(:,1)=CP_cumsum_all;%TW
%
% load('H:\global-PV-wind\ANS\ID_pro3.mat') % pro FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
% a = unique(line_IX_all(:,10));
% for i = 1:size(ID_pro,1)
%     [m,n]=find(line_IX_all(:,10)==ID_pro(i,1));
%     if ~isempty(m)
%         line_IX_all(m,16) = ID_pro(i,3);
%     end
% end

load('H:\global-PV-wind\ANS\cho_county_all_pro2_8_2070.mat'); %
% 0是无需储能，1是机械储能，-1是化学储能
cho(cho==-1)=2;
cho(cho==0)=2;
% 1是机械储能，2是化学储能，没有的默认为化学储能

% loss rate of electricity
ef1 = 0.967; % power loss rate during the period of circle transportation
ef2 = 0.983; % power loss rate during UHV transportation within the country
ef3 = [0.7;0.99*0.85]; % power loss rate during charge, discharge and storage period
% ef4 = 0.983; % power loss rate during UHV transportation between countries
ef4 = 1; % power loss rate during UHV transportation between countries

reg = 11;
tic

load('H:\global-PV-wind\Data\region_connect0811.mat')  % 1 起点；2 终点。
region_connect = region_connect0811;
clear region_connect0811
% load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all.mat')
% optpowerunit_IX_IX(:,1)=optpowerunit_IX_IX(:,1).*utilize_ratio2060_UHV_storage_inter;
load('H:\global-PV-wind\ANS\distance_UHV_Station_all.mat')  % 千米

% linescapacity=zeros(numlines,1); % capacity of electricity transmissions
% linescapacity_plant=zeros(numlines,size(distance_UHV_Station,1));
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

etrans_cou_num=zeros(numpowerunit,192,192); %
etrans_cou1_num=zeros(numpowerunit,192,192); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
% etrans_cou1_num2=zeros(numpowerunit,192,192);
% etrans_cou1_num3=zeros(numpowerunit,192,192);

powergenerat_monhour = zeros(288,numpowerunit);

P_REG = zeros(288,size(UHV_Station_country,1)); % 省的有效电量
P_country = zeros(288,192); % 各国的有效电量
P_REG1 = zeros(288,numpowerunit);
P_country1 = zeros(288,numpowerunit); % 仅有UHV情景下各国的有效发电量
P_country2 = zeros(288,numpowerunit); % 有UHV和storage的情景下各国的有效发电量
P_region2 = zeros(288,numpowerunit);
% P_region = zeros(288,numpowerunit); % 有UHV，storage和国际电力运输时世界的有效发电量
P_region = zeros(288,reg); % 有UHV，storage和国际电力运输时世界的有效发电量
P_world = zeros(288,numpowerunit); % 有UHV，storage和国际电力运输时世界的有效发电量

storage_UHV_storage = zeros(288,numpowerunit);
storage_UHV_storage_regional = zeros(288,numpowerunit);
storage_UHV_storage_inter = zeros(288,numpowerunit);

% linescapacity=zeros(size(distance_UHV_Station,1),size(distance_UHV_Station,1)); % capacity of electricity transmissions
% linescapacity_plant=zeros(size(distance_UHV_Station,1),size(distance_UHV_Station,1),10000);
% linescapacity_CN=zeros(77,1); % capacity of electricity transmissions
% linescapacity_plant_CN=zeros(77,numpowerunit);

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
    country2=optpowerunit_IX_IX(i2,41); % country of power  % 没有UHV的国家是-1
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
                
                if CHO_MAC(i,1)~=0
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                etrans_cou_num(i,country,country) = etrans_cou_num(i,country,country)+P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
                etrans_cou1_num(i,country,country) = etrans_cou1_num(i,country,country)+P_ef1;
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
            %         if P_REG(zz,REG)<powerdemand_monhour2060_REG(zz,REG) && P_country(zz,country) < powerdemand_monhour(zz,country)  && sum(sum(P_country(:,country))) < sum(sum(powerdemand_monhour(:,country)))% zz时pro省有效电量<zz时pro省需电量
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
                
                if CHO_MAC(i,1)~=0
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                etrans_cou_num(i,country,country) = etrans_cou_num(i,country,country)+P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
                etrans_cou1_num(i,country,country) = etrans_cou1_num(i,country,country)+P_ef1;
                 end
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
                
                if CHO_MAC(i,1)~=0
                etrans(region,region) = etrans(region,region) + P_ef1/ef1;
                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1;
                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1;
                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1;
                etrans_REG(REG,REG) = etrans_REG(REG,REG) + P_ef1/ef1;
                etrans_REG_type(REG,REG,type(i,1)) = etrans_REG_type(REG,REG,type(i,1)) + P_ef1/ef1;
                etrans_cou_num(i,country,country) = etrans_cou_num(i,country,country)+P_ef1/ef1;
                
                etrans1(region,region) = etrans1(region,region) + P_ef1;
                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                etrans_REG1(REG,REG) = etrans_REG1(REG,REG) + P_ef1;
                etrans_REG_type1(REG,REG,type(i,1)) = etrans_REG_type1(REG,REG,type(i,1)) + P_ef1;
                etrans_cou1_num(i,country,country) = etrans_cou1_num(i,country,country)+P_ef1;
                end
            end
        end
    end
end
1
etrans_cou1_num1 = etrans_cou1_num/288*8760; % 实际有效发电量, TWh/year
etrans_cou_num1 = etrans_cou_num/288*8760; % 全部发电量, TWh/year
save('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_8_2070.mat','etrans_cou1_num1','-v7.3')% 实际有效发电量, TWh/year
save('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_8_2070.mat','etrans_cou_num1','-v7.3')% 全部发电量, TWh/year
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
    country2=optpowerunit_IX_IX(i2,41); % country of power  % 没有UHV的国家是-1
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
                                
                                if CHO_MAC(i,1)>1
                                etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2;
                                etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2;
                                etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2;
                                etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2;
                                etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                                etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                                etrans_cou_num(i,country,country) = etrans_cou_num(i,country,country)+P_ef1/ef1/ef2;
                                
                                etrans1(region,region) = etrans1(region,region) + P_ef1;
                                etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                                etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                                etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                                etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                                etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                                etrans_cou1_num(i,country,country) = etrans_cou1_num(i,country,country)+P_ef1;
                                end
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
                            
                            if CHO_MAC(i,1)>1
                            etrans(region,region) = etrans(region,region) + P_ef1/ef1/ef2;
                            etrans_type(region,region,type(i,1)) = etrans_type(region,region,type(i,1)) + P_ef1/ef1/ef2;
                            etrans_cou(country,country) = etrans_cou(country,country) + P_ef1/ef1/ef2;
                            etrans_cou_type(country,country,type(i,1)) = etrans_cou_type(country,country,type(i,1)) + P_ef1/ef1/ef2;
                            etrans_REG(REG,REG2) = etrans_REG(REG,REG2) + P_ef1/ef1/ef2;
                            etrans_REG_type(REG,REG2,type(i,1)) = etrans_REG_type(REG,REG2,type(i,1)) + P_ef1/ef1/ef2;
                            etrans_cou_num(i,country,country) = etrans_cou_num(i,country,country)+P_ef1/ef1/ef2;
                            
                            etrans1(region,region) = etrans1(region,region) + P_ef1;
                            etrans_type1(region,region,type(i,1)) = etrans_type1(region,region,type(i,1)) + P_ef1;
                            etrans_cou1(country,country) = etrans_cou1(country,country) + P_ef1;
                            etrans_cou_type1(country,country,type(i,1)) = etrans_cou_type1(country,country,type(i,1)) + P_ef1;
                            etrans_REG1(REG,REG2) = etrans_REG1(REG,REG2) + P_ef1;
                            etrans_REG_type1(REG,REG2,type(i,1)) = etrans_REG_type1(REG,REG2,type(i,1)) + P_ef1;
                            etrans_cou1_num(i,country,country) = etrans_cou1_num(i,country,country)+P_ef1;
                            end
                        end
                    end
                end
            end
            
            
        end
        
    end
    
    
    %     i
end
2
etrans_cou1_num1 = etrans_cou1_num/288*8760; % 实际有效发电量, TWh/year
etrans_cou_num1 = etrans_cou_num/288*8760; % 全部发电量, TWh/year
save('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_8_2070.mat','etrans_cou1_num1','-v7.3')% 实际有效发电量, TWh/year
save('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_8_2070.mat','etrans_cou_num1','-v7.3')% 全部发电量, TWh/year

