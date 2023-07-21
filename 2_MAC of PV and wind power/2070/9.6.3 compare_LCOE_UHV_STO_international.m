tic
clear
load('H:\Global PV and wind\ANS\countryy_IX_county_all_pro2_8_2070.mat') 
load('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_inter_county_all_pro2_8_2070.mat')  % Average annual electricity storage, TWh/year
storage_year_plant = storage_inter_year_plant_UHV_storage;
load('H:\Global PV and wind\Data\discount_country.mat')
discount = discount/100;
discount(35)=0.05;

lifetime_power = 25; % year
lifetime_inverter=10; % renewed per 10 years
OMratio_PV = 0.01;
OMratio_wind=0.03;

discountinverter=zeros(192,1);
for i = 1:192
for t=1:floor(lifetime_power/lifetime_inverter)
    discountinverter(i,1)=discountinverter(i,1)+1/(1+discount(i,1))^((t-1)*lifetime_inverter);
end
end

discount40yr=zeros(4,192);
for i = 1:192
for t=1:lifetime_power
    discount40yr(4,i)=discount40yr(4,i)+1/(1+discount(i,1))^(t-1);
end
end

discountinverter_1 = zeros(size(storage_year_plant,1),1);
discount40yr_1 = zeros(4,size(storage_year_plant,1));
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discountinverter_1(m,1) = discountinverter(i,1);
    discount40yr_1(4,m) = discount40yr(4,i);
end
clear discountinverter
clear discount40yr
discountinverter=discountinverter_1;
discount40yr=discount40yr_1;

disrate_PV = ones(4,size(storage_year_plant,1))+OMratio_PV.*discount40yr;
disrate_inverter_PV = (discountinverter.*ones(1,4))'+OMratio_PV.*discount40yr;
disrate_wind = ones(4,size(storage_year_plant,1))+OMratio_wind.*discount40yr;

discount_country = zeros(size(storage_year_plant,1),1);
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discount_country(m,1) = discount(i,1);
end

degration_PV = 0.00;
degration_onshorewind = 0.00;
degration_offshorewind = 0.0; 
degrat40yr_PV=zeros(4,size(storage_year_plant,1));
degrat40yr_onshorewind=zeros(4,size(storage_year_plant,1));
degrat40yr_offshorewind=zeros(4,size(storage_year_plant,1));
for t=1:lifetime_power
    degrat40yr_PV(4,:)=degrat40yr_PV(4,:)+((1-degration_PV).^(t-1)./(1+discount_country).^(t-1))';
    degrat40yr_onshorewind(4,:)=degrat40yr_onshorewind(4,:)+((1-degration_onshorewind).^(t-1)./(1+discount_country).^(t-1))';
    degrat40yr_offshorewind(4,:)=degrat40yr_offshorewind(4,:)+((1-degration_offshorewind).^(t-1)./(1+discount_country).^(t-1))';
end

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
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind;

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

[is,pos]=ismember(index_mineral_pv_time2,optpowerunit_PV(:,40));
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
optpowerunit_PV2 = optpowerunit_PV(pos,:);
lines_IX_PV2 = lines_IX_PV(pos,:);
clear optpowerunit_PV
clear lines_IX_PV
optpowerunit_PV = optpowerunit_PV2;
lines_IX_PV = lines_IX_PV2;
clear optpowerunit_PV2
clear index_mineral_pv
clear lines_IX_PV2

[is,pos]=ismember(index_mineral_ons_time2,optpowerunit_onshorewind(:,40));
optpowerunit_onshorewind2 = optpowerunit_onshorewind(pos,:);
lines_IX_onshorewind2 = lines_IX_onshorewind(pos,:);
clear optpowerunit_onshorewind
clear lines_IX_onshorewind
optpowerunit_onshorewind = optpowerunit_onshorewind2;
lines_IX_onshorewind = lines_IX_onshorewind2;
clear optpowerunit_onshorewind2
clear lines_IX_onshorewind2

[is,pos]=ismember(index_mineral_off_time2,optpowerunit_offshorewind(:,40));
optpowerunit_offshorewind2 = optpowerunit_offshorewind(pos,:);
lines_IX_offshorewind2 = lines_IX_offshorewind(pos,:);
clear optpowerunit_offshorewind
clear lines_IX_offshorewind
optpowerunit_offshorewind = optpowerunit_offshorewind2;
lines_IX_offshorewind = lines_IX_offshorewind2;
clear optpowerunit_offshorewind2
clear lines_IX_offshorewind2

%%
optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
clear optpowerunit_PV
clear optpowerunit_onshorewind
clear optpowerunit_offshorewind
lines_IX_offshorewind(:,12:18)=0;
lines_IX = [lines_IX_PV;lines_IX_onshorewind;lines_IX_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
    lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
end

%%
load('H:\Global PV and wind\ANS\CO2_mecha_county_all_pro2_8_2070.mat'); % USD2019/kWh
load('H:\Global PV and wind\ANS\B_mecha_county_all_pro2_8_2070.mat'); % USD2019/kWh
load('H:\Global PV and wind\ANS\LCOEE_mecha_county_all_pro2_8_2070.mat'); % USD2019/kWh
CO2_mecha = CO2all_c_utilize_trans_storage;
B_mecha = B_utilize_trans_storage;
LCOE_mecha = LCOEE_all_utilize_trans_storage;

load('H:\Global PV and wind\ANS\CO2_Battery_county_all_pro2_8_2070.mat'); % USD2019/kWh
load('H:\Global PV and wind\ANS\B_Battery_county_all_pro2_8_2070.mat'); % USD2019/kWh
load('H:\Global PV and wind\ANS\LCOEE_Battery_county_all_pro2_8_2070.mat'); % USD2019/kWh
CO2_Battery = CO2all_c_utilize_trans_storage;
B_Battery = B_utilize_trans_storage;
LCOE_Battery = LCOEE_all_utilize_trans_storage;

CO2=zeros(size(CO2_Battery,1),1);
B=zeros(size(CO2_Battery,1),1);
for i = 1:size(CO2_Battery,1)
    if LCOE_mecha(i,1)<LCOE_Battery(i,1)
        cho(i,1)=1; %机械储能
        B(i,1) = B_mecha(i,1);
        CO2(i,1) = CO2_mecha(i,1);
    end
    if LCOE_mecha(i,1)>LCOE_Battery(i,1)
        cho(i,1)=-1; % 化学储能
        B(i,1) = B_Battery(i,1);
        CO2(i,1) = CO2_Battery(i,1);
    end
    if LCOE_mecha(i,1)==LCOE_Battery(i,1)
        cho(i,1)=0; % 无需储能
        B(i,1) = B_Battery(i,1);
        CO2(i,1) = CO2_Battery(i,1);
    end
end
save('H:\Global PV and wind\ANS\cho_county_all_pro2_8_2070.mat','cho'); % 

CHO_ALL = optpowerunit_IX(:,35);
CHO_ALL(:,2) = optpowerunit_IX(:,40);
CHO_ALL(:,3) = cho;
save('H:\Global PV and wind\ANS\CHO_ALL_pro2_8_2070.mat','CHO_ALL'); % 

B_Battery(B_Battery==-Inf)=0;
CC_battery = -B_Battery.*CO2_Battery;
CC_mecha = -B_mecha.*CO2_mecha;
CC = -B.*CO2;

(sum(CC_mecha)-sum(CC))/sum(CC_mecha)
(sum(CC_battery)-sum(CC))/sum(CC_battery)

load('G:\Code1123\ANS1114\storage_max_plant_xz.mat')   % Maximum charging power, TWh/h  
load('G:\Code1123\ANS1114\storage_year_plant_xz.mat')  % Average annual electricity storage, TWh/year
Ip = storage_max_plant/(0.99*0.85*0.983*0.967)*10^9; % kW  power capacity
