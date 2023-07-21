tic
clear;
rmb2us=1/6.8967; % RMB to USD2019
lifetime_power=25;
load('H:\Global PV and wind\Data\discount_country.mat')
discount = discount/100;
discount(35) = 0.05;
load('H:\Global PV and wind\Data\fossilfuel_emissionfactor.mat')  %kg CO2/kWh
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
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % 电厂编号
clear powerunit_IX_offshorewind

powerunit_num_IX_PV_ori = powerunit_num_IX_PV(:,5);
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind(:,5);
off_pro_IX_ori = off_pro_IX(:,2);
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
a = zeros(size(powerunit_num_IX_PV,1),1);
b = zeros(size(powerunit_num_IX_onshorewind,1),1);
c = zeros(size(off_pro_IX,1),1);
for country =  1:1:192
    [m,n] = find(powerunit_num_IX_PV(:,5)==country);
    a(m,1) = region_ID(country,1);
    [m,n] = find(powerunit_num_IX_onshorewind(:,5)==country);
    b(m,1) = region_ID(country,1);
    [m,n] = find(off_pro_IX(:,2)==country);
    c(m,1) = region_ID(country,1);
end
powerunit_num_IX_PV(:,5)=a;
powerunit_num_IX_onshorewind(:,5)=b;
powerunit_num_IX_offshorewind(:,5)=c;
clear a
clear b
clear c
clear m

%%
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\REG_plant_pv_all.mat')  % 各电厂所在的国家，REG的ID(1-4)，UHV Station的ID
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\REG_plant_ons_all.mat')  % 各电厂所在的国家，REG的ID，UHV Station的ID
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\REG_plant_off_all.mat')  % 各电厂所在的国家，REG的ID，UHV Station的ID
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

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\Country_ID_withUHV.mat')  % 有UHV Sation所在的国家ID
Country_ID_withUHV2 = zeros(192,1);
Country_ID_withUHV2(floor(Country_ID_withUHV),1)=1;
clear Country_ID_withUHV
load('H:\Global PV and wind\ANS\UHV_Station_country_all.mat')
% 1Substation; 2行；3列；4国家ID; 
% 5region ID; 6pro ID(0-3638); 7该序号所分配的power demand (2060年电气化后全部需电量，TWh/year); 8REG(1-4)

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

[is,pos]=ismember(index_mineral_pv_time2,optpowerunit_PV(:,40));
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
optpowerunit_PV2 = optpowerunit_PV(pos,:);
powerunit_num_IX_PV2 = powerunit_num_IX_PV(pos,:);
powerunit_num_IX_PV_ori2 = powerunit_num_IX_PV_ori(pos,:);
lines_IX_PV2 = lines_IX_PV(pos,:);
clear optpowerunit_PV
clear powerunit_num_IX_PV
clear lines_IX_PV
optpowerunit_PV = optpowerunit_PV2;
powerunit_num_IX_PV = powerunit_num_IX_PV2;
powerunit_num_IX_PV_ori = powerunit_num_IX_PV_ori2;
lines_IX_PV = lines_IX_PV2;
clear optpowerunit_PV2
clear powerunit_num_IX_PV2
clear index_mineral_pv
clear lines_IX_PV2
clear powerunit_num_IX_PV_ori2

[is,pos]=ismember(index_mineral_ons_time2,optpowerunit_onshorewind(:,40));
optpowerunit_onshorewind2 = optpowerunit_onshorewind(pos,:);
powerunit_num_IX_onshorewind2 = powerunit_num_IX_onshorewind(pos,:);
lines_IX_onshorewind2 = lines_IX_onshorewind(pos,:);
powerunit_num_IX_onshorewind_ori2 = powerunit_num_IX_onshorewind_ori(pos,:);
clear optpowerunit_onshorewind
clear powerunit_num_IX_onshorewind
clear powerunit_num_IX_onshorewind_ori
clear lines_IX_onshorewind

optpowerunit_onshorewind = optpowerunit_onshorewind2;
powerunit_num_IX_onshorewind = powerunit_num_IX_onshorewind2;
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind_ori2;
lines_IX_onshorewind = lines_IX_onshorewind2;
clear optpowerunit_onshorewind2
clear powerunit_num_IX_onshorewind2
clear index_mineral_ons
clear lines_IX_onshorewind2
clear powerunit_num_IX_onshorewind_ori2

[is,pos]=ismember(index_mineral_off_time2,optpowerunit_offshorewind(:,40));
optpowerunit_offshorewind2 = optpowerunit_offshorewind(pos,:);
off_pro_IX2 = off_pro_IX(pos,:);
off_pro_IX_ori2 = off_pro_IX_ori(pos,:);
powerunit_num_IX_offshorewind2 = powerunit_num_IX_offshorewind(pos,:);
lines_IX_offshorewind2 = lines_IX_offshorewind(pos,:);
clear optpowerunit_offshorewind
clear powerunit_num_IX_offshorewind
clear off_pro_IX
clear off_pro_IX_ori
clear lines_IX_offshorewind
optpowerunit_offshorewind = optpowerunit_offshorewind2;
off_pro_IX = off_pro_IX2;
off_pro_IX_ori = off_pro_IX_ori2;
powerunit_num_IX_offshorewind = powerunit_num_IX_offshorewind2;
lines_IX_offshorewind = lines_IX_offshorewind2;
clear optpowerunit_offshorewind2
clear off_pro_IX2
clear off_pro_IX_ori2
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
powerunit_num_IX_offshorewind(:,6)=0;
powerunit_num_IX = [powerunit_num_IX_PV;powerunit_num_IX_onshorewind;powerunit_num_IX_offshorewind];
powerunit_country_IX = [powerunit_num_IX_PV_ori;powerunit_num_IX_onshorewind_ori;off_pro_IX_ori];
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
powerunit_num_IX_IX(:,1:5)=powerunit_num_IX(IX,1:5);
countryy_IX = powerunit_country_IX(IX,1);
clear optpowerunit
clear lines_IX
clear powerunit_num_IX
clear IX
% [B,IX]=sort(optpowerunit_IX(:,20),1);
line_IX_all = [lines(:,1:15);lines_IX_IX];
clear lines
clear lines_IX_IX

powerunit_num_IX_IX_region = powerunit_num_IX_IX(:,5);
save('H:\Global PV and wind\ANS\powerunit_num_IX_IX_region_county_all_pro.mat','powerunit_num_IX_IX_region'); % TWh/year
save('H:\Global PV and wind\ANS\countryy_IX_county_all_pro.mat','countryy_IX'); % TWh/year
save('H:\Global PV and wind\ANS\optpowerunit_IX_county_all_pro.mat','optpowerunit_IX'); % TWh/year

type_IX = optpowerunit_IX(:,35);
save('H:\Global PV and wind\ANS\type_IX_county_all_pro.mat','type_IX'); %
