tic
clear;
load('H:\Global PV and wind\Data\r_module_hardwarecost_pv.mat');
% Module占代码中Moduleprice2的比例	
% 硬件成本占代码中Moduleprice2的比例
r_module_hardwarecost_ons = 0.704960836;
r_fixednomodule_hardwarecost_ons = 0.7493-r_module_hardwarecost_ons;

load('H:\Global PV and wind\ANS\cho_county_all_pro2_8_2070.mat'); % 
load('H:\Global PV and wind\ANS\countryy_IX_county_all_pro2_8_2070.mat') 
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];
load('H:\Global PV and wind\ANS\cost_trans_8_2040_xzxz.mat')% million $和optpowerunit_IX顺序一样
load('H:\Global PV and wind\ANS\CP_trans_8_2040_xzxz.mat')% MW

load('H:\Global PV and wind\ANS\utilize_ratio_8_2040_xzxz.mat')  % power use efficiency with global electric transport and storage
utilize_ratio2060 =utilize_ratio2060_UHV_storage_inter;
utilize_ratio2060(find(isnan(utilize_ratio2060)==1))=0;
clear utilize_ratio2060_UHV_storage_inter
numpowerunit = size(utilize_ratio2060,1);

load('H:\Global PV and wind\ANS\storage_max_plant_8_2040_xzxz.mat')   % Maximum charging power, TWh/h  
load('H:\Global PV and wind\ANS\storage_year_plant_8_2040_xzxz.mat')  % Average annual electricity storage, TWh/year
storage_max_plant = storage_max_plant_UHV_storage_inter;
storage_year_plant = storage_inter_year_plant_UHV_storage;
clear storage_max_plant_UHV_storage_inter
clear storage_year_plant_UHV_storage_inter

lifetime_s = 50; % year
lifetime_power = 25; % year
Ip = storage_max_plant/(0.7*0.983*0.967)*10^9; % kW  power capacity
Ie = storage_year_plant/(0.7*0.983*0.967)*10^9*15/6000; % kWh  energy capacity
for i = 1:numpowerunit
    if Ip(i,1)>Ie(i,1)
        Ie(i,1) = Ip(i,1);
    end
end

p_dis = storage_year_plant*10^9/(0.983*0.967); % annual discharge, kWh/year 
p_char = storage_year_plant*10^9/(0.7*0.983*0.967); % annual charge, kWh/year 
Co = 1.5/1000; % $/kWh
Cp=1200; % $/kW
Ce = 100; % 2018$/kWh

load('H:\Global PV and wind\Data\discount_country.mat')
discount = discount/100;
discount(35)=0.05;
discount1yr_1 = zeros(192,1);
for i = 1:192
for t=1:lifetime_power
    discount1yr_1(i,1)=discount1yr_1(i,1)+1/(1+discount(i))^(t-1);
end
end

discount1yr = zeros(numpowerunit,1);
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discount1yr(m,1) = discount1yr_1(i,1);
end
clear discount1yr_1
r_y = lifetime_power/lifetime_s;
%

Cost_m = (Ce.*Ie+Cp.*Ip)/10^6*r_y+Co*(p_dis+p_char)/10^6.*discount1yr; %million $
Cost_m_ini = (Ce.*Ie+Cp.*Ip)/10^6*r_y; %million $
Cost_m_OM = Co*(p_dis+p_char)/10^6.*discount1yr; %million $
clear discount1yr
Cost_mechanical=Cost_m;
clear Cost_m

%
% 0是无需储能，1是机械储能，-1是化学储能
[mmm,nnn]=find(cho==1);
Cost_mechanical1=zeros(numpowerunit,1);
Ip_mechanical=zeros(numpowerunit,1);
Cost_m_ini1 = zeros(numpowerunit,1); %
Cost_m_OM1 = zeros(numpowerunit,1); %
cost_trans=zeros(numpowerunit,1);
CP_trans=zeros(numpowerunit,1);
utilize_ratio2060_2=zeros(numpowerunit,1);
Cost_mechanical1(mmm,1) = Cost_mechanical(mmm,1);
Cost_m_ini1(mmm,1) = Cost_m_ini(mmm,1); %million $
Cost_m_OM1(mmm,1) = Cost_m_OM(mmm,1); %million $
Ip_mechanical(mmm,1) = Ip(mmm,1); %kW
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1); 
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);
% clear Inc
% clear CP_trans_IX
% clear cost_trans_IX
% clear Cost_mechanical
% clear Ip_mechanical
% clear utilize_ratio2060

save('H:\Global PV and wind\ANS\Cost_mechanical1_centralcase_county_all_pro2_8_2040_8_2040_2s_2050_test6xzsuper.mat','Cost_mechanical1') 

%%
% load('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all.mat')  % power use efficiency with global electric transport and storage
% utilize_ratio2060 =utilize_ratio2060_UHV_storage_inter;
% clear utilize_ratio2060_UHV_storage_inter
% 
% load('H:\Global PV and wind\ANS\storage_max_plant_UHV_storage_inter_county_all.mat')   % Maximum charging power, TWh/h  
% load('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_inter_county_all.mat')  % Average annual electricity storage, TWh/year
% storage_max_plant = storage_max_plant_UHV_storage_inter;
% storage_year_plant = storage_inter_year_plant_UHV_storage;
% clear storage_max_plant_UHV_storage_inter
% clear storage_year_plant_UHV_storage_inter
% load('H:\Global PV and wind\ANS\cost_trans_IX_all_battery_mechanical_all.mat')% million $和optpowerunit_IX顺序一样
% load('H:\Global PV and wind\ANS\CP_trans_IX_all_battery_mechanical_all.mat')% MW

lifetime_s = 15; % year
lifetime_power = 25; % year
Ip = storage_max_plant/(0.99*0.85*0.983*0.967)*10^9; % kW  power capacity
Ie = storage_year_plant/(0.99*0.85*0.983*0.967)*10^9*lifetime_s/6000; % kWh  energy capacity
for i = 1:numpowerunit
    if Ip(i,1)>Ie(i,1)
        Ie(i,1) = Ip(i,1);
    end
end

p_dis = storage_year_plant*10^9/(0.983*0.967); % annual discharge, kWh/year 
p_char = storage_year_plant*10^9/(0.99*0.85*0.983*0.967); % annual charge, kWh/year 
Co = 1.5/1000; % $/kWh
Cp=595.73; % $/kW
Ce = 345; % 2018$/kWh

discount1yr_1 = zeros(192,1);
for i = 1:192
for t=1:lifetime_power
    discount1yr_1(i,1)=discount1yr_1(i,1)+1/(1+discount(i))^(t-1);
end
end
discount2yr_1 = zeros(192,1);
for i = 1:192
for t=16
    discount2yr_1(i,1) = discount2yr_1(i,1)+1/(1+discount(i))^(t-1);
end
end

discount1yr = zeros(numpowerunit,1);
discount2yr = zeros(numpowerunit,1);
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discount1yr(m,1) = discount1yr_1(i,1);
    discount2yr(m,1) = discount2yr_1(i,1);
end
clear discount1yr_1
clear discount2yr_1

r_y = (lifetime_power-lifetime_s)/lifetime_s;
%
load('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat'); %% rIp=[595.73, 374.45, 327.22,280.78]./595.73;
unitmin_CHA = unitmin;
unitmin(unitmin==1) = 6;
unitmin(unitmin==2) = 7;
rIp=[595.73,440,374.45,374.45-(374.45-327.22)/2, 327.22,327.22-(327.22-280.78)/2, 280.78, 280.78-(280.78-234.343)/2, 234.343, 234.343-(280.78-234.343)/2,234.343-(280.78-234.343)]./595.73;
rIe=[345,242,198,186, 174, 161, 149, 149-(149-124)/2, 124, 124-(149-124)/2, 124-(149-124)]./345;
for i = 1:10
    [m,n]=find(unitmin==i);
    learnprice_sto_stage_p(m,1) = rIp(i);
    learnprice_sto_stage_e(m,1) = rIe(i);
end

Cost_storage1 = (Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6+(Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6.*discount2yr*r_y+Co*(p_dis+p_char)/10^6.*discount1yr; %million $
Cost_storage1_ini = (Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6+(Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6.*discount2yr*r_y; %million $
Cost_storage1_OM = Co*(p_dis+p_char)/10^6.*discount1yr; %million $
clear discount1yr
Cost_storage=Cost_storage1;
clear Cost_storage1

% 0是无需储能，1是机械储能，-1是化学储能
[mmm,nnn]=find(cho==-1);
Cost_storage1=zeros(numpowerunit,1);
Ip_storage=zeros(numpowerunit,1);
Cost_storage1(mmm,1) = Cost_storage(mmm,1);
Cost_storage1_ini1(mmm,1) = Cost_storage1_ini(mmm,1);
Cost_storage1_OM1(mmm,1) = Cost_storage1_OM(mmm,1);
Ip_storage(mmm,1) = Ip(mmm,1); %kW
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1); 
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);
% clear Inc
clear Cost_storage
clear Ip_storage
Cost_storage=Cost_storage1;
clear Cost_storage1
save('H:\Global PV and wind\ANS\Cost_storage_centralcase_county_all_pro2_8_2040_2s_2050_test6xzsuper.mat','Cost_storage') 


Cost_sto_ini = Cost_m_ini1;
Cost_sto_OM = Cost_m_OM1;
[m,nnn]=find(cho==-1);
Cost_sto_ini(m) = Cost_storage1_ini(m); % million $/y
Cost_sto_OM(m) = Cost_storage1_OM1(m); % million $/y
save('H:\Global PV and wind\ANS\Cost_sto_ini_centralcase_county_all_pro2_8_2040_2s_2050_test6xzsuper.mat','Cost_sto_ini') 
save('H:\Global PV and wind\ANS\Cost_sto_OM_centralcase_county_all_pro2_8_2040_2s_2050_test6xzsuper.mat','Cost_sto_OM') 

%%
[mmm,nnn]=find(cho==0);
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1); 
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);
clear utilize_ratio2060
clear CP_trans_IX
clear cost_trans_IX
% utilize_ratio2060_trans_plant_alone = utilize_ratio2060_2;
% utilize_ratio2060_trans_storage_plant_alone = utilize_ratio2060_2;
utilize_ratio2060 = utilize_ratio2060_2;
clear utilize_ratio2060_2

%%
LR_PV2 = 0.18; %0.18; % 0.37; %  0.37 PV module
CP0_PV = sum(sum(CP_trans))/10; % MW
learnprice_trans(1)=1;
for i = 2:10
    [m,n]=find(unitmin_CHA<=i-1);
    learnprice_trans(i,1) = ((sum(CP_trans(m))+CP0_PV)./CP0_PV).^(log(1- LR_PV2)/log(2));
end
learnprice_trans(11,1) = ((sum(CP_trans)+CP0_PV)./CP0_PV).^(log(1- LR_PV2)/log(2));
learnprice_trans(learnprice_trans>1)=1;

%%
CO2_C=0.2727;
lifetime_power=25;
lifetime_inverter=10; % renewed per 10 years
OMratio_PV = 0.01;
OMratio_wind = 0.03;
discountinverter=zeros(192,1);

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

% % 
% load('H:\Global PV and wind\ANS\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
% for i = 1:size(off_pro_IX,1)
%     a = off_pro_IX(i,1);
%     [m1,n1]=find(ID_pro(:,1)==a);
%     off_pro_IX(i,2) = unique(ID_pro(m1,3));
%     i
% end
% clear ID_pro2

% powerunit_country_IX = [powerunit_num_IX_PV(:,5);powerunit_num_IX_onshorewind(:,5);off_pro_IX(:,2)];
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
off_pro_IX(:,3)=c;
clear a
clear b
clear c
clear m

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

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
clear powerunit_num_IX_PV_ori
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
% powerunit_num_IX_offshorewind2 = powerunit_num_IX_offshorewind(pos,:);
lines_IX_offshorewind2 = lines_IX_offshorewind(pos,:);
clear optpowerunit_offshorewind
% clear powerunit_num_IX_offshorewind
clear off_pro_IX
clear off_pro_IX_ori
clear lines_IX_offshorewind
optpowerunit_offshorewind = optpowerunit_offshorewind2;
off_pro_IX = off_pro_IX2;
off_pro_IX_ori = off_pro_IX_ori2;
% powerunit_num_IX_offshorewind = powerunit_num_IX_offshorewind2;
lines_IX_offshorewind = lines_IX_offshorewind2;
clear optpowerunit_offshorewind2
clear off_pro_IX2
clear off_pro_IX_ori2
% clear powerunit_num_IX_offshorewind2
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
powerunit_num_IX = [powerunit_num_IX_PV(:,5);powerunit_num_IX_onshorewind(:,5);off_pro_IX(:,3)];
powerunit_country_IX = [powerunit_num_IX_PV_ori;powerunit_num_IX_onshorewind_ori;off_pro_IX_ori];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
powerunit_IX(:,1)=IX;
optpowerunit_IX(:,1:40)=optpowerunit(IX,1:40);
lines_IX_IX(:,1:15)=lines_IX(IX,1:15);
powerunit_num_IX_IX(:,:)=powerunit_num_IX(IX,:);
powerunit_country_IX_IX(:,:)=powerunit_country_IX(IX,:);
clear optpowerunit
clear lines_IX
clear powerunit_num_IX
clear powerunit_country_IX
clear IX

% %%
% [m,n]=find(utilize_ratio2060<=0);
% utilize_ratio2060(m,:)=[];
% powerunit_IX(m,:)=[];
% optpowerunit_IX(m,:)=[];
% lines_IX_IX(m,:)=[];
% powerunit_num_IX_IX(m,:)=[];
% powerunit_country_IX_IX(m,:)=[];
% numpowerunit = size(optpowerunit_IX,1);
% Cost_storage(m,:)=[];
% Cost_mechanical1(m,:)=[];
% Cost_mechanical(m,:)=[];
% cost_trans(m,:)=[];
% countryy_IX(m,:)=[];

%%
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
% for t=1:30
%     discount40yr(3,1)=discount40yr(3,1)+1/(1+discount)^(t-1);
% end
% for t=1:20
%     discount40yr(2,1)=discount40yr(2,1)+1/(1+discount)^(t-1);
% end
% for t=1:10
%     discount40yr(1,1)=discount40yr(1,1)+1/(1+discount)^(t-1);
% end
discountinverter_1 = zeros(numpowerunit,1);
discount40yr_1 = zeros(4,numpowerunit);
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discountinverter_1(m,1) = discountinverter(i,1);
    discount40yr_1(4,m) = discount40yr(4,i);
end
clear discountinverter
clear discount40yr
discountinverter=discountinverter_1;
discount40yr=discount40yr_1;


% disrate_PV = ones(4,numpowerunit)+OMratio_PV.*discount40yr;
% disrate_inverter_PV = discountinverter.*ones(4,numpowerunit)+OMratio_PV.*discount40yr;
% disrate_wind = ones(4,numpowerunit)+OMratio_wind.*discount40yr;
disrate_PV = ones(4,numpowerunit)+OMratio_PV.*discount40yr;
disrate_inverter_PV = (discountinverter.*ones(1,4))'+OMratio_PV.*discount40yr;
disrate_wind = ones(4,numpowerunit)+OMratio_wind.*discount40yr;

discount_country = zeros(numpowerunit,1);
for i = 1:192
    [m,n] = find(countryy_IX==i);
    discount_country(m,1) = discount(i,1);
end

degration_PV = 0.00;
degration_onshorewind = 0.00;
degration_offshorewind = 0.0; 
degrat40yr_PV=zeros(4,numpowerunit);
degrat40yr_onshorewind=zeros(4,numpowerunit);
degrat40yr_offshorewind=zeros(4,numpowerunit);
for t=1:lifetime_power
    degrat40yr_PV(4,:)=degrat40yr_PV(4,:)+((1-degration_PV).^(t-1)./(1+discount_country).^(t-1))';
    degrat40yr_onshorewind(4,:)=degrat40yr_onshorewind(4,:)+((1-degration_onshorewind).^(t-1)./(1+discount_country).^(t-1))';
    degrat40yr_offshorewind(4,:)=degrat40yr_offshorewind(4,:)+((1-degration_offshorewind).^(t-1)./(1+discount_country).^(t-1))';
end

%%
CP_cumsum_all = cumsum(optpowerunit_IX(:,30))/10^6;
CF=sum(optpowerunit_IX(:,1))./sum(optpowerunit_IX(:,30))/8760*10^6;
[B,IX]=sort(optpowerunit_IX(:,20),1);

load('H:\Global PV and wind\Data\fossilfuel_emissionfactor.mat')  % kg CO2/kWh
EF_country = zeros(numpowerunit,1);
for i = 1:1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    EF_country(m,1) = fossilfuel_emissionfactor(i,1);
end
clear fossilfuel_emissionfactor
%%
learnprice_region2_2 = ones(192,6,11);
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro22_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,2) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro32_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,3) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro42_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,4) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro52_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,5) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro62_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,6) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro72_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,7) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro82_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,8) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro92_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,9) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro102_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,10) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro112_8_2040xz.mat'); % 10 regions * 6
learnprice_region2_2(:,:,11) = learnprice_region;
clear learnprice_region
learnprice_region_nodiffuse = learnprice_region2_2;
clear learnprice_region2_2

learnprice_region2_2 = ones(192,6,11);
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro22_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,2) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro32_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,3) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro42_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,4) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro52_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,5) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro62_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,6) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro72_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,7) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro82_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,8) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro92_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,9) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro102_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,10) = learnprice_region;
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro112_8_2040xzsuper.mat'); % 10 regions * 6
learnprice_region2_2(:,:,11) = learnprice_region;
clear learnprice_region
learnprice_region = learnprice_region2_2;


line_IX_all = [lines(:,1:15);lines_IX_IX];
% [m,n] = find(line_IX_all(:,9)==0);
% line_IX_all(m,9) = 1362; % offshorewind, Shenzhen

% PV
% LR_PV1 = 0.2; 
% LR_PV2 = 0.18;
CP0_PV = 253000;  % cumulative capacity potential in China in 2020, MW
CP_PV = zeros(numpowerunit,1);
CP_PV(1,1) = CP0_PV;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==1
        CP_PV(i+1,1) =  CP_PV(i,1) + optpowerunit_IX(i,30);
    else
        CP_PV(i+1,1) =  CP_PV(i,1);
    end    
end
CP0_PV_cum = CP_PV(2:end,1)-CP0_PV;
% r_pv = CP_PV./CP0_PV;
% for i = 1:size(r_pv,1)
%         LR_PV(i,1) = LR_PV1;
%         learnprice_PV(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_PV(i,1))/log(2));
%         LR_PV_2(i,1) = LR_PV2;
%         learnprice_PV_2(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_PV_2(i,1))/log(2));        
% end

for i = 1:size(CP_PV,1)-1
    CP_PV22(i,1) = CP_PV(i+1,1)-CP_PV(1,1);
    CP_PV22(i,2) = (CP_PV(i+1,1)-CP_PV(1,1))/(CP_PV(end,1)-CP_PV(1,1));
end
cost_PV = zeros(numpowerunit,19);
cost_PV_alone = zeros(numpowerunit,30);
cost_PV_trans = zeros(numpowerunit,1);
CO2_PV_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==1
    cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_PV_alone(1,1:30) = optpowerunit_IX(1,1:30);
    cost_PV_trans(1,1) = cost_trans(1,1);
    CO2_PV_alone(1,1:3)  = optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==1
        cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_PV_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        cost_PV_trans(i,1) = cost_trans(i,1);
        CO2_PV_alone(i,1:3)  = optpowerunit_IX(i,8:10);
    else
        cost_PV(i,11:19) =  cost_PV(i-1,11:19);
        cost_PV_alone(i,1:30) =  0;
        cost_PV_trans(i,1) = 0;
        CO2_PV_alone(i,1:3)  = 0;
    end    
end

% onshorewind
% LR_onshorewind1= 0.074; 
% LR_onshorewind2= 0.18;
CP0_onshorewind = 272010;% MW
CP_onshorewind = zeros(numpowerunit,1);
CP_onshorewind(1,1) = CP0_onshorewind;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==2
        CP_onshorewind(i+1,1) =  CP_onshorewind(i,1) + optpowerunit_IX(i,30);
    else
        CP_onshorewind(i+1,1) =  CP_onshorewind(i,1);
    end    
end
CP0_onshorewind_cum = CP_onshorewind(2:end,1)-CP0_onshorewind;
% r_onshorewind = CP_onshorewind./CP0_onshorewind;
% for i = 1:size(r_onshorewind,1)
%         LR_onshorewind(i,1) = LR_onshorewind1;
%         learnprice_onshorewind(i,1) = (CP_onshorewind(i,1)./CP0_onshorewind).^(log(1- LR_onshorewind(i,1))/log(2));
%         LR_onshorewind_2(i,1) = LR_onshorewind2;
%         learnprice_onshorewind_2(i,1) = (CP_onshorewind(i,1)./CP0_onshorewind).^(log(1- LR_onshorewind_2(i,1))/log(2));
% end

for i = 1:size(CP_onshorewind,1)-1
    CP_onshorewind22(i,1) = CP_onshorewind(i+1,1)-CP_onshorewind(1,1);
    CP_onshorewind22(i,2) = (CP_onshorewind(i+1,1)-CP_onshorewind(1,1))/(CP_onshorewind(end,1)-CP_onshorewind(1,1));
end
cost_onshorewind = zeros(numpowerunit,19);
cost_onshorewind_alone = zeros(numpowerunit,30);
cost_onshorewind_trans = zeros(numpowerunit,1);
CO2_onshorewind_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==2
    cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_onshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
    cost_onshorewind_trans(1,1) = cost_trans(1,1);
    CO2_onshorewind_alone(1,1:3) =  optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==2
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_onshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        cost_onshorewind_trans(i,1) = cost_trans(i,1);
        CO2_onshorewind_alone(i,1:3) =  optpowerunit_IX(i,8:10);
    else
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
        cost_onshorewind_alone(i,1:30) =  0;
        cost_onshorewind_trans(i,1) = 0;
        CO2_onshorewind_alone(i,1:3) =  0;
    end    
end

% offhorewind
% LR_offshorewind1 = 0.074;
% LR_offshorewind2 = 0.18;
CP0_offshorewind = 8990;
CP_offshorewind = zeros(numpowerunit,1);
CP_offshorewind(1,1) = CP0_offshorewind;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==3
        CP_offshorewind(i+1,1) =  CP_offshorewind(i,1) + optpowerunit_IX(i,30);
    else
        CP_offshorewind(i+1,1) =  CP_offshorewind(i,1);
    end    
end
CP0_offshorewind_cum = CP_offshorewind(2:end,1)-CP0_offshorewind;
for i = 1:size(CP_offshorewind,1)-1 
    CP_offshorewind22(i,1) = CP_offshorewind(i+1,1)-CP_offshorewind(1,1);
    CP_offshorewind22(i,2) = (CP_offshorewind(i+1,1)-CP_offshorewind(1,1))/(CP_offshorewind(end,1)-CP_offshorewind(1,1));
end
% r_offshorewind = CP_offshorewind./CP0_offshorewind;
CP_offshorewind=CP_offshorewind;

% for i = 1:size(r_offshorewind,1)
%         LR_offshorewind(i,1) = LR_offshorewind1;
%         learnprice_offshorewind(i,1) = (CP_offshorewind(i,1)./CP0_offshorewind).^(log(1- LR_offshorewind(i,1))/log(2));
%         LR_offshorewind_2(i,1) = LR_offshorewind2;
%         learnprice_offshorewind_2(i,1) = (CP_offshorewind(i,1)./CP0_offshorewind).^(log(1- LR_offshorewind_2(i,1))/log(2));
% end

cost_offshorewind = zeros(numpowerunit,10);
cost_offshorewind_alone = zeros(numpowerunit,30);
cost_offshorewind_trans = zeros(numpowerunit,1);
if optpowerunit_IX(1,35)==3
    cost_offshorewind(1,6:10) =  optpowerunit_IX(1,6:10);
    cost_offshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
    cost_offshorewind_trans(1,1) = cost_trans(1,1);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==3
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10) + optpowerunit_IX(i,6:10);
        cost_offshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);        
        cost_offshorewind_trans(i,1) = cost_trans(i,1);
    else
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10);
        cost_offshorewind_alone(i,1:30) =  0;
        cost_offshorewind_trans(i,1) = 0;
    end    
end

optpowerunit1= zeros(numpowerunit,30);
optpowerunit1(1,30)=optpowerunit_IX(1,30);
for i=2:numpowerunit
    optpowerunit1(i,30)=optpowerunit1(i-1,30)+optpowerunit_IX(i,30); % cumulative capacity
end
optpowerunit_IX(:,36)=optpowerunit1(:,30)/optpowerunit1(end,30);

%%
lcoee_PV = zeros(numpowerunit,1);
lcoee_onshorewind = zeros(numpowerunit,1);
lcoee_offshorewind = zeros(numpowerunit,1);

CP_PV_alone(:,1) = cost_PV_alone(:,30); % Capacity potential, MW
CP_onshorewind_alone(:,1) = cost_onshorewind_alone(:,30); % Capacity potential, MW
CP_offshorewind_alone(:,1) = cost_offshorewind_alone(:,30); % Capacity potential, MW

CO2_PV_year(:,1)= cost_PV_alone(:,8)./CO2_C+ cost_PV_alone(:,9)./CO2_C; % Mton CO2
CO2_onshorewind_year(:,1)=cost_onshorewind_alone(:,8)./CO2_C+cost_onshorewind_alone(:,9)./CO2_C; % Mton CO2
CO2_offshorewind_year(:,1)=cost_offshorewind_alone(:,5)./CO2_C; % Mton CO2

CO2_PV_year_abated(:,1)= cost_PV_alone(:,8)./CO2_C; % Mton CO2
CO2_onshorewind_year_abated(:,1)=cost_onshorewind_alone(:,8)./CO2_C; % Mton CO2
CO2_offshorewind_year_abated(:,1)=cost_offshorewind_alone(:,5)./CO2_C; % Mton CO2
CO2_year_abated = CO2_PV_year_abated+CO2_onshorewind_year_abated+CO2_offshorewind_year_abated;

CO2_PV_year_ls(:,1)= cost_PV_alone(:,9)./CO2_C; % Mton CO2
CO2_onshorewind_year_ls(:,1)= cost_onshorewind_alone(:,9)./CO2_C; % Mton CO2
CO2_offshorewind_year_ls(:,1)= 0; % Mton CO2

CO2_PV_lcc(:,1)= cost_PV_alone(:,10)./CO2_C; % Mton CO2
CO2_onshorewind_lcc(:,1)=cost_onshorewind_alone(:,10)./CO2_C; % Mton CO2

% 4  land cost
lcoee_PV(:,4) = cost_PV_alone(:,16)./disrate_PV(4,:)' ; % cost million USD 
% 6  abated annual CO2 emission
lcoee_PV(:,6) = cost_PV_alone(:,17); % cost million USD 

% 7  land carbon sink
lcoee_PV(:,7) = cost_PV_alone(:,18); % cost million USD 
% 8  land use change
lcoee_PV(:,8) = cost_PV_alone(:,19); % cost million USD 
% 4  land cost
lcoee_onshorewind(:,4) = cost_onshorewind_alone(:,16)./disrate_wind(4,:)' ; % cost million USD 
% 6  abated annual CO2 emission
lcoee_onshorewind(:,6) = cost_onshorewind_alone(:,17); % cost million USD 

% 7  land carbon sink
lcoee_onshorewind(:,7) = cost_onshorewind_alone(:,18); % cost million USD 
% 8  land use change
lcoee_onshorewind(:,8) = cost_onshorewind_alone(:,19); % cost million USD 
% 6  abated annual CO2 emission
lcoee_offshorewind(:,6) = cost_offshorewind_alone(:,7); % cost million USD 

% power generation
phh_PV(:,1) = cost_PV_alone(:,1).*degrat40yr_PV(4,:)'; 

% sto_PV(:,1) = cost_PV_alone(:,1).*(1-utilize_ratio2060_trans_plant_alone); 
phh_PV_alone(:,1) = cost_PV_alone(:,1); 
% sto_PV_sto2(:,1) = cost_PV_alone(:,1).*(1-utilize_ratio2060); 

% power generation
phh_offshorewind(:,1) = cost_offshorewind_alone(:,1).*degrat40yr_offshorewind(4,:)'; 

% power generation
phh_onshorewind_utilize_trans_storage(:,1) = cost_onshorewind_alone(:,1).*degrat40yr_onshorewind(4,:)'.*utilize_ratio2060; 

load('H:\Global PV and wind\Data\initialcost_ratio_country_0111low_off.mat')  % offshorewind的成本与China的对比
for region = 1:1:192
    module_price_offshore = initialcost_ratio_country(region).*1.008; % USD/W
    for nnna = 1:10
    [mm,nn]=find(powerunit_country_IX_IX==region & unitmin==nnna);
% PV
% 1  fixed system cost
lcoee_PV(mm,1) = cost_PV_alone(mm,14)./disrate_PV(4,mm)'.*r_module_hardwarecost_pv(region,1) .*learnprice_region(region,1,nnna) + cost_PV_alone(mm,14)./disrate_PV(4,mm)'.*(r_module_hardwarecost_pv(region,2)-r_module_hardwarecost_pv(region,1)) .*learnprice_region(region,2,nnna)+ cost_PV_alone(mm,14)./disrate_PV(4,mm)'.*(1-r_module_hardwarecost_pv(region,2)) .*learnprice_region_nodiffuse(region,2,nnna) + cost_PV_alone(mm,15)./disrate_inverter_PV(4,mm)'.*discountinverter(mm,1) .*learnprice_region(region,2,nnna) ; % cost million USD 
% 2  substation
lcoee_PV(mm,2) = cost_PV_alone(mm,12)./disrate_PV(4,mm)'.*learnprice_region(region,2,nnna); % cost million USD 
% 3  Power plant cable
lcoee_PV(mm,3) = cost_PV_alone(mm,13)./disrate_PV(4,mm)'.*learnprice_region(region,2,nnna); % cost million USD 
% 5  connection to national grid
lcoee_PV(mm,5) = cost_PV_alone(mm,11)./disrate_PV(4,mm)'.*learnprice_region(region,2,nnna); % cost million USD 

% onshorewind
% 1  fixed system cost
lcoee_onshorewind(mm,1) = cost_onshorewind_alone(mm,14)./disrate_wind(4,mm)'.*r_module_hardwarecost_ons.*learnprice_region(region,3,nnna)+cost_onshorewind_alone(mm,14)./disrate_wind(4,mm)'.*r_fixednomodule_hardwarecost_ons.*learnprice_region(region,4,nnna)+cost_onshorewind_alone(mm,14)./disrate_wind(4,mm)'.*(1-r_module_hardwarecost_ons-r_fixednomodule_hardwarecost_ons).*learnprice_region_nodiffuse(region,4,nnna); % cost million USD 
% 2  substation
lcoee_onshorewind(mm,2) = cost_onshorewind_alone(mm,12)./disrate_wind(4,mm)'.*learnprice_region(region,4,nnna); % cost million USD 
% 3  Power plant cable
lcoee_onshorewind(mm,3) = cost_onshorewind_alone(mm,13)./disrate_wind(4,mm)'.*learnprice_region(region,4,nnna); % cost million USD 
% 5  connection to national grid
lcoee_onshorewind(mm,5) = cost_onshorewind_alone(mm,11)./disrate_wind(4,mm)'.*learnprice_region(region,4,nnna); % cost million USD 

% offshorewind
% 1  module cost
lcoee_offshorewind(mm,1) = module_price_offshore.*CP_offshorewind_alone(mm,1).*learnprice_region(region,5,nnna); % cost million USD 
% 2  others
lcoee_offshorewind(mm,2) = (cost_offshorewind_alone(mm,6)./disrate_wind(4,mm)'-module_price_offshore.*CP_offshorewind_alone(mm,1)).*learnprice_region(region,6,nnna); % cost million USD 
% UHV transport
lcoee_PV(mm,11) = cost_PV_trans(mm,1).*learnprice_trans(nnna,1);
lcoee_onshorewind(mm,11) = cost_onshorewind_trans(mm,1).*learnprice_trans(nnna,1);
lcoee_offshorewind(mm,11) = cost_offshorewind_trans(mm,1).*learnprice_trans(nnna,1);
    end
end
% 9  O&M cost
lcoee_PV(:,9) = lcoee_PV(:,11).*(disrate_PV(4,:)-1)'+lcoee_PV(:,4).*(disrate_PV(4,:)-1)'+(sum(lcoee_PV(:,1:3),2)+lcoee_PV(:,5)).*(disrate_PV(4,:)-1)'; % cost million USD 
lcoee_onshorewind(:,9) = lcoee_onshorewind(:,11).*(disrate_wind(4,:)-1)'+lcoee_onshorewind(:,4).*(disrate_wind(4,:)-1)'+(sum(lcoee_onshorewind(:,1:3),2)+lcoee_onshorewind(:,5)).*(disrate_wind(4,:)-1)'; % cost million USD 
lcoee_offshorewind(:,9) = lcoee_offshorewind(:,11).*(disrate_wind(4,:)-1)'+sum(lcoee_offshorewind(:,1:2),2).*(disrate_wind(4,:)-1)'; % cost million USD 

CPPP=CP_PV_alone;
[m,n]=find(CP_PV_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_PV_alone(sub2ind(size(CP_PV_alone), m, n));
[m,n]=find(CP_onshorewind_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_onshorewind_alone(sub2ind(size(CP_onshorewind_alone), m, n));
[m,n]=find(CP_offshorewind_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_offshorewind_alone(sub2ind(size(CP_offshorewind_alone), m, n));
CP_ALL = CPPP/10^6; %TW
%% 求证一下cost_PV_alone(:,6)是否全为0
% lcoee_PV_utilize_trans_storage = cost_PV_alone(:,17).*utilize_ratio2060; % cost million USD 
phh_PV_utilize_trans_storage(:,1) = cost_PV_alone(:,1).*degrat40yr_PV(4,:)'.*utilize_ratio2060; 
% lcoee_onshorewind_utilize_trans_storage = cost_onshorewind_alone(:,17).*utilize_ratio2060; % cost million USD 
phh_onshorewind_utilize_trans_storage(:,1) = cost_onshorewind_alone(:,1).*degrat40yr_onshorewind(4,:)'.*utilize_ratio2060; 
% lcoee_offshorewind_utilize_trans_storage = cost_offshorewind_alone(:,7).*utilize_ratio2060; % cost million USD 
phh_offshorewind_utilize_trans_storage(:,1) = cost_offshorewind_alone(:,1).*degrat40yr_offshorewind(4,:)'.*utilize_ratio2060; 

phhall_all = cost_PV_alone(:,1);
[m,n]=find(cost_onshorewind_alone(:,1)~=0);
phhall_all(sub2ind(size(phhall_all), m, n))= cost_onshorewind_alone(sub2ind(size(cost_onshorewind_alone), m, n));
[m,n]=find(cost_offshorewind_alone(:,1)~=0);
phhall_all(sub2ind(size(phhall_all), m, n))= cost_offshorewind_alone(sub2ind(size(cost_offshorewind_alone), m, n));


phhall_utilize_trans_storage = phh_PV_utilize_trans_storage;
phhall_utilize_trans_storage2 = phh_PV_utilize_trans_storage;
lcoeeall_utilize_trans_storage = lcoee_PV;
lcoeeall_utilize_trans_storage2 = lcoee_PV;
% lcoeeall_utilize_trans_storage(:,6) = lcoee_PV_utilize_trans_storage;
[m,n]=find(phh_onshorewind_utilize_trans_storage~=0);
phhall_utilize_trans_storage(sub2ind(size(phhall_utilize_trans_storage), m, n))= phh_onshorewind_utilize_trans_storage(sub2ind(size(phh_onshorewind_utilize_trans_storage), m, n));
phhall_utilize_trans_storage2(sub2ind(size(phhall_utilize_trans_storage2), m, n))= phh_onshorewind_utilize_trans_storage(sub2ind(size(phh_onshorewind_utilize_trans_storage), m, n));
[m,n]=find(lcoee_onshorewind~=0);
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, n))= lcoee_onshorewind(sub2ind(size(lcoee_onshorewind), m, n));
lcoeeall_utilize_trans_storage2(sub2ind(size(lcoeeall_utilize_trans_storage2), m, n))= lcoee_onshorewind(sub2ind(size(lcoee_onshorewind), m, n));
% lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, 6*ones(size(m,1),1)))= lcoee_onshorewind_utilize_trans_storage(sub2ind(size(lcoee_onshorewind_utilize_trans_storage), m, ones(size(m,1),1)));
[m,n]=find(phh_offshorewind_utilize_trans_storage~=0);
phhall_utilize_trans_storage(sub2ind(size(phhall_utilize_trans_storage), m, n))= phh_offshorewind_utilize_trans_storage(sub2ind(size(phh_offshorewind_utilize_trans_storage), m, n));
[m,n]=find(lcoee_offshorewind~=0);
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, n))= lcoee_offshorewind(sub2ind(size(lcoee_offshorewind), m, n));
% lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, 6*ones(size(m,1),1)))= lcoee_offshorewind_utilize_trans_storage(sub2ind(size(lcoee_offshorewind_utilize_trans_storage), m, ones(size(m,1),1)));

LCOEE_all_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:11),2)+Cost_storage+Cost_mechanical1)./phhall_utilize_trans_storage(:,1)/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
% LCOEE_all_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:9),2)+Cost_storage+Cost_mechanical1)./(phhall_utilize_trans_storage(:,1)./utilize_ratio2060)/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
lcoeeall_utilize_trans_storage(:,10) = Cost_storage+Cost_mechanical1;
lcoeeall_utilize_trans_storage2(:,10) = Cost_storage+Cost_mechanical1;
lcoeeall_storage(:,1) = Cost_storage;
lcoeeall_storage(:,2) = Cost_mechanical1;

lcoeeall_utilize_trans_storage_type = lcoeeall_utilize_trans_storage(:,1:9);
lcoeeall_utilize_trans_storage_type(:,10) = Cost_storage;
lcoeeall_utilize_trans_storage_type(:,11) = Cost_mechanical1;
lcoeeall_utilize_trans_storage_type(:,12) = lcoeeall_utilize_trans_storage(:,11);
save('H:\Global PV and wind\ANS\lcoeeall_utilize_trans_storage_type2_8_2040_2s_2050_test6xzsuper.mat','lcoeeall_utilize_trans_storage_type', '-v7.3')  % billion kWh / year
% 10是battery storage，11是mechanical storage, 12是UHV cost

%%
load('H:\Global PV and wind\Data\Powerdmeand_country_coal.mat'); % TWh,行：1-192号国家，193行是全球 % 列：2000-2020
load('H:\Global PV and wind\Data\Powerdmeand_country_gas.mat'); % TWh,行：1-192号国家，193行是全球
load('H:\Global PV and wind\Data\Powerdmeand_country_oil.mat'); % TWh,行：1-192号国家，193行是全球
a = Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil;
[m,n]=find(a(:,21)==0);
clear a

r1=Powerdmeand_country_coal./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
r2=Powerdmeand_country_gas./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
r3=Powerdmeand_country_oil./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
clear Powerdmeand_country_coal
clear Powerdmeand_country_gas
clear Powerdmeand_country_oil

r1(m,21)=r1(193,21);
r2(m,21)=r2(193,21);
r3(m,21)=r3(193,21);
r = r1(:,21)+r2(:,21)+r3(:,21);

rr(:,1) = r1(1:192,21);
rr(:,2) = r2(1:192,21);
rr(:,3) = r3(1:192,21);
% 2020年各国coal gas oil发电占fossil发电的比例 
clear r
clear r1
clear r2
clear r3
% rr(35,1) = 0.9567;
% rr(35,2) = 0.9974-0.9567;
% rr(35,3) = 1-0.9974;


% [B,ixx]=sort(LCOEE_all_utilize_trans_storage,'descend');  
% optpowerunit_IX_ix = optpowerunit_IX(ixx,1); % LCOE从低到高的电厂
% powerunit_country_IX_IX_ix = powerunit_country_IX_IX(ixx,1); % LCOE从低到高的电厂
% 
% load('H:\Global PV and wind\ANS\EF_coal.mat'); % kg CO2/kWh
% load('H:\Global PV and wind\ANS\EF_gas.mat'); % kg CO2/kWh
% load('H:\Global PV and wind\ANS\EF_oil.mat'); % kg CO2/kWh
% load('H:\Global PV and wind\ANS\Price_coal.mat'); % USD/kWh, 第一列是均值，第二列std
% load('H:\Global PV and wind\ANS\Price_gas.mat'); % USD/kWh, 第一列是均值，第二列std
% load('H:\Global PV and wind\ANS\Price_oil.mat'); % USD/kWh, 第一列是均值，第二列std
% for i = 1:192
%     [m,n]=find(powerunit_country_IX_IX_ix==i);
%     rp = cumsum(optpowerunit_IX_ix(m,1))./sum(optpowerunit_IX_ix(m,1));
%     [~,Index1] = min(abs(rp-rr(i,1))); % coal
%     for jj = 1:Index1
%         ixx_2_1(m(jj),1)=1;
%         emissionfactor_1(m(jj),1)=EF_coal(i,1); % coal, kg CO2/kWh
%         Price_1(m(jj),1)=Price_coal(i,1); % coal, USD/kWh
%     end
%     [~,Index2] = min(abs(rp-rr(i,2)-rr(i,1))); % gas
%     for jj = Index1+1:Index2
%         ixx_2_1(m(jj),1)=2;
%         emissionfactor_1(m(jj),1)=EF_gas(i,1); % natural gas, kg CO2/kWh
%         Price_1(m(jj),1)=Price_gas(i,1); % natural gas, USD/kWh
%     end
%     [~,Index3] = min(abs(rp-rr(i,3)-rr(i,2)-rr(i,1))); %oil
%     for jj = Index2+1:Index3
%         ixx_2_1(m(jj),1)=3;
%         emissionfactor_1(m(jj),1)=EF_oil(i,1); % oil, kg CO2/kWh
%         Price_1(m(jj),1)=Price_oil(i,1); % oil, USD/kWh
%     end
% end
% for i = 1:size(ixx,1)
%     ixx_2(ixx(i),1) = ixx_2_1(i,1);
%     emissionfactor(ixx(i),1) = emissionfactor_1(i,1);
%     EP(ixx(i),1) = Price_1(i,1);
% end
% save('H:\Global PV and wind\ANS\emissionfactor_8_2040_2s_2050_test6xzsuper.mat','emissionfactor', '-v7.3')  
% clear ixx_2_1
% clear emissionfactor_1
% clear Price_1
% clear EF_gas
% clear EF_oil
% clear EF_coal
% clear Price_gas
% clear Price_oil
% clear Price_coal
% 
% clear optpowerunit_IX_ix
% clear powerunit_country_IX_IX_ix
% clear Cost_Aba
% clear rp
% clear ixx
% clear Index1
% clear Index2
% clear Index3
load('H:\Global PV and wind\ANS\EP_8_2040_2s_2050_test6xz.mat')  % EP
load('H:\Global PV and wind\ANS\emissionfactor_8_2040_2s_2050_test6xz.mat'); % emissionfactor


%%
CO2_PV_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_PV_alone(:,8)./CO2_C.*degrat40yr_PV(4,:)'.*utilize_ratio2060); % Mton CO2
CO2_onshorewind_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_onshorewind_alone(:,8)./CO2_C.*degrat40yr_onshorewind(4,:)'.*utilize_ratio2060); % Mton CO2
CO2_offshorewind_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_offshorewind_alone(:,5)./CO2_C.*degrat40yr_offshorewind(4,:)'.*utilize_ratio2060); % Mton CO2
% CO2_PV_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_PV_alone(:,8)./CO2_C.*degrat40yr_PV(4,:)'.*utilize_ratio2060+ cost_PV_alone(:,9)./CO2_C.*degrat40yr_PV(4,:)'+ cost_PV_alone(:,10)./CO2_C); % Mton CO2
% CO2_onshorewind_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_onshorewind_alone(:,8)./CO2_C.*degrat40yr_onshorewind(4,:)'.*utilize_ratio2060+cost_onshorewind_alone(:,9)./CO2_C.*degrat40yr_onshorewind(4,:)'+ cost_onshorewind_alone(:,10)./CO2_C); % Mton CO2
% CO2_offshorewind_all_utilize_trans_storage(:,1)= emissionfactor./EF_country.*(cost_offshorewind_alone(:,5)./CO2_C.*degrat40yr_offshorewind(4,:)'.*utilize_ratio2060); % Mton CO2

CO2_all_utilize_trans_storage = CO2_PV_all_utilize_trans_storage; % Mton CO2
[m,n]=find(CO2_onshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_onshorewind_all_utilize_trans_storage(sub2ind(size(CO2_onshorewind_all_utilize_trans_storage), m, n));
[m,n]=find(CO2_offshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_offshorewind_all_utilize_trans_storage(sub2ind(size(CO2_offshorewind_all_utilize_trans_storage), m, n));

%%
% load('H:\Global PV and wind\ANS\EP_country.mat'); % electricity price, $/kWh
CO2_all_utilize_trans_storage(CO2_all_utilize_trans_storage>0)=0;
B_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:11),2)-phhall_utilize_trans_storage(:,1).*EP*10^3)./(-CO2_all_utilize_trans_storage); %USD/t CO2
[m,n] = find(CO2_all_utilize_trans_storage==0);
B_utilize_trans_storage(m,1) = Inf;

LCOEE_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:11),2))./phhall_utilize_trans_storage(:,1)/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
save('H:\Global PV and wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','LCOEE_all_utilize_trans_storage', '-v7.3')  % billion kWh / year
cost_all_withUHVcost =  (sum(lcoeeall_utilize_trans_storage(:,1:11),2)); % million USD2019
save('H:\Global PV and wind\ANS\cost_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','cost_all_withUHVcost', '-v7.3')  % billion kWh / year
inicost_all_withUHVcost =  (sum(lcoeeall_utilize_trans_storage(:,1:8),2)); % million USD2019
save('H:\Global PV and wind\ANS\inicost_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','inicost_all_withUHVcost', '-v7.3')  % billion kWh / year
phhall_utilize_trans_storage111 = phhall_utilize_trans_storage(:,1);
% save('H:\Global PV and wind\ANS\Ph_UHV_STO_reg2_8_2040_2s_2050_test6xzsuper.mat','phhall_utilize_trans_storage111', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\Ph_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','phhall_utilize_trans_storage111', '-v7.3')  % billion kWh / year
phhall_utilize_trans_storage111_year = phhall_utilize_trans_storage111./degrat40yr_onshorewind(4,:)';
save('H:\Global PV and wind\ANS\Ph_UHV_STO_INT_year_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','phhall_utilize_trans_storage111_year', '-v7.3')  % billion kWh / year

cost = sum(lcoeeall_utilize_trans_storage(:,1:11),2); % million 2019 USD
% save('H:\Global PV and wind\ANS\costtt_UHV_STO_reg2_8_2040_2s_2050_test6xzsuper.mat','cost', '-v7.3')  % 
save('H:\Global PV and wind\ANS\costtt_UHV_STO_INT_county_all_withUHVcost2_8_2040_2s_2050_test6xzsuper.mat','cost', '-v7.3')  % million 2019 USD 

[B,IX_utilize_trans_storage]=sort(LCOEE_all_utilize_trans_storage);
CP = optpowerunit_IX(:,30)/10^6; %TW
CP_IX=CP(IX_utilize_trans_storage,1);
CP_IXs=cumsum(CP_IX);

[m,n]=find(CO2_PV_year~=0);
B_utilize_trans_storage_pv= B_utilize_trans_storage(sub2ind(size(B_utilize_trans_storage), m, n));
[m,n]=find(CO2_onshorewind_year~=0);
B_utilize_trans_storage_onshorewind= B_utilize_trans_storage(sub2ind(size(B_utilize_trans_storage), m, n));
[m,n]=find(CO2_offshorewind_year~=0);
B_utilize_trans_storage_offshorewind= B_utilize_trans_storage(sub2ind(size(B_utilize_trans_storage), m, n));

%% 成本分解
cost_all = lcoeeall_utilize_trans_storage(:,1:11)/10^6; % trillion 2019 USD
cost_all(:,12:13)  = lcoeeall_storage(:,1:2)/10^6; % trillion 2019 USD
save('H:\Global PV and wind\ANS\cost_all_withUHVcost1004_pro2_8_2040_2s_2050_test6xzsuper.mat','cost_all', '-v7.3')  % trillion 2019 USD

for cou = 1:192
    [m,n]=find(powerunit_country_IX_IX==cou);
    if size(m,1)~=1
    lcoeeall_cou(1:11,cou) = sum(lcoeeall_utilize_trans_storage(m,1:11))'; % million USD
    lcoeeall_cou(12:13,cou) = sum(lcoeeall_storage(m,1:2))'; %
    ph_cou(1,cou) = sum(phhall_utilize_trans_storage(m,1))'; %
    end
    if size(m,1)==1
    lcoeeall_cou(1:11,cou) = lcoeeall_utilize_trans_storage(m,1:11)'; % million USD
    lcoeeall_cou(12:13,cou) = lcoeeall_storage(m,1:2)'; %
    ph_cou(1,cou) = phhall_utilize_trans_storage(m,1); %
    end
end
lcoeeall_cou'/1000; % billion USD
ph_cou';
lcoeeall_cou'./ph_cou'/1000; % USD/kWh
% 
%%
CO2_PV_year_abated = emissionfactor./EF_country.*CO2_PV_year_abated;
CO2_onshorewind_year_abated = emissionfactor./EF_country.*CO2_onshorewind_year_abated;
CO2_offshorewind_year_abated = emissionfactor./EF_country.*CO2_offshorewind_year_abated;

CO2_PV_year_ls = emissionfactor./EF_country.*CO2_PV_year_ls;
CO2_onshorewind_year_ls = emissionfactor./EF_country.*CO2_onshorewind_year_ls;
CO2_offshorewind_year_ls = emissionfactor./EF_country.*CO2_offshorewind_year_ls;

CO2_PV_lcc = emissionfactor./EF_country.*CO2_PV_lcc;
CO2_onshorewind_lcc = emissionfactor./EF_country.*CO2_onshorewind_lcc;

IX=[1:1:numpowerunit]';
CO2_year_abated = CO2_PV_year_abated; % Mton CO2
[m,n]=find(CO2_onshorewind_year_abated~=0);
CO2_year_abated(sub2ind(size(CO2_year_abated), m, n))= CO2_onshorewind_year_abated(sub2ind(size(CO2_onshorewind_year_abated), m, n));
[m,n]=find(CO2_offshorewind_year_abated~=0);
CO2_year_abated(sub2ind(size(CO2_year_abated), m, n))= CO2_offshorewind_year_abated(sub2ind(size(CO2_offshorewind_year_abated), m, n));

CO2_year_ls = CO2_PV_year_ls; % Mton CO2
[m,n]=find(CO2_onshorewind_year_ls~=0);
CO2_year_ls(sub2ind(size(CO2_year_ls), m, n))= CO2_onshorewind_year_ls(sub2ind(size(CO2_onshorewind_year_ls), m, n));
[m,n]=find(CO2_offshorewind_year_ls~=0);
CO2_year_ls(sub2ind(size(CO2_year_ls), m, n))= CO2_offshorewind_year_ls(sub2ind(size(CO2_offshorewind_year_ls), m, n));

CO2_year_lcc = CO2_PV_lcc; % Mton CO2
[m,n]=find(CO2_onshorewind_lcc~=0);
CO2_year_lcc(sub2ind(size(CO2_year_lcc), m, n))= CO2_onshorewind_lcc(sub2ind(size(CO2_onshorewind_lcc), m, n));
CO2_year_cum_lcc = -cumsum(CO2_year_lcc)./10^3;

CO2_year_utilize_trans_storage = CO2_year_abated.*utilize_ratio2060+CO2_year_ls;
CO2_year_utilize_trans_storage(CO2_year_utilize_trans_storage>0)=0;
%%
CO2all_c_utilize_trans_storage = zeros(size(IX,1),1);
CO2all_c_utilize_trans_storage= CO2_year_utilize_trans_storage(sub2ind(size(CO2_year_utilize_trans_storage), IX, ones(size(IX,1),1)));

[m,n]=find(CO2_PV_year~=0);
CO2all_c_utilize_trans_storage_pv= CO2all_c_utilize_trans_storage(sub2ind(size(CO2all_c_utilize_trans_storage), m, n));
[m,n]=find(CO2_onshorewind_year~=0);
CO2all_c_utilize_trans_storage_onshorewind= CO2all_c_utilize_trans_storage(sub2ind(size(CO2all_c_utilize_trans_storage), m, n));
[m,n]=find(CO2_offshorewind_year~=0);
CO2all_c_utilize_trans_storage_offshorewind= CO2all_c_utilize_trans_storage(sub2ind(size(CO2all_c_utilize_trans_storage), m, n));


save('H:\Global PV and wind\ANS\CP_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CP', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\CO2_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage', '-v7.3')  % billion kWh / year


save('H:\Global PV and wind\ANS\CO2_PV_year_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2_PV_year', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\CO2_onshorewind_year_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2_onshorewind_year', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\CO2_offshorewind_year_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2_offshorewind_year', '-v7.3')  % billion kWh / year

% save('H:\Global PV and wind\ANS\CP_UHV_STO_reg2_8_2040_2s_2050_test6xzsuper.mat','CP', '-v7.3')  % billion kWh / year
% save('H:\Global PV and wind\ANS\CO2_UHV_STO_reg2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage', '-v7.3')  % billion kWh / year
% save('H:\Global PV and wind\ANS\B_UHV_STO_reg2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage', '-v7.3')  % billion kWh / year
% save('H:\Global PV and wind\ANS\CP_UHV_STO2_8_2040_2s_2050_test6xzsuper.mat','CP', '-v7.3')  % billion kWh / year
% save('H:\Global PV and wind\ANS\CO2_UHV_STO2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage', '-v7.3')  % billion kWh / year
% save('H:\Global PV and wind\ANS\B_UHV_STO2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage', '-v7.3')  % billion kWh / year

save('H:\Global PV and wind\ANS\CO2_UHV_STO_INT_pv_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage_pv', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\B_UHV_STO_INT_pv_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage_pv', '-v7.3')  % billion kWh / year

save('H:\Global PV and wind\ANS\CO2_UHV_STO_INT_ons_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage_onshorewind', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\B_UHV_STO_INT_ons_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage_onshorewind', '-v7.3')  % billion kWh / year

save('H:\Global PV and wind\ANS\CO2_UHV_STO_INT_off_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','CO2all_c_utilize_trans_storage_offshorewind', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\B_UHV_STO_INT_off_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','B_utilize_trans_storage_offshorewind', '-v7.3')  % billion kWh / year

B_utilize_trans_storage(B_utilize_trans_storage==-inf)=inf;
[B_utilize_trans_storage_IX,IX]=sort(B_utilize_trans_storage);
save('H:\Global PV and wind\ANS\IX_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2050_test6xzsuper.mat','IX', '-v7.3')  % 

lcoeeall_utilize_trans_storage_IX2=lcoeeall_utilize_trans_storage(IX,1:11);
phhall_utilize_trans_storage_IX2=phhall_utilize_trans_storage(IX,1);
degrat40yr_offshorewind2 = degrat40yr_offshorewind(:,IX);
phhall_all_IX2=phhall_all(IX,1);
utilize_ratio2060_IX2=utilize_ratio2060(IX,1);
CO2_year_utilize_trans_storage_IX2=CO2_year_utilize_trans_storage(IX,1);
CO2_year_ls_IX2=CO2_year_ls(IX,1);
CO2_lcc_IX2=CO2_year_lcc(IX,1);
CP_ALL_IX2=CP_ALL(IX,1); %TW
ph_plant_trans_storage = phhall_utilize_trans_storage_IX2./degrat40yr_offshorewind2(4,:)'; %TWh/year
LCOEE_plant_trans_storage = (sum(lcoeeall_utilize_trans_storage_IX2(:,1:11),2))./phhall_utilize_trans_storage_IX2/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
CO2_plant_trans_storage = -CO2_year_utilize_trans_storage_IX2;
CO2_plant_ls_year = CO2_year_ls_IX2;
CO2_plant_lcc = CO2_lcc_IX2;
% ph_plant_trans_storage_all = phhall_utilize_trans_storage_IX2./degrat40yr_offshorewind(4,:)'./utilize_ratio2060_IX2; %TWh/year
ph_plant_trans_storage_all = phhall_all_IX2; %TWh/year
ph_plant_trans_storage_all(isnan(ph_plant_trans_storage_all)==1)=0;

INNN(1,1)=0;
[~,Index1] = min(abs(B_utilize_trans_storage_IX+150));
INNN(2,1)=Index1;
[~,Index1] = min(abs(B_utilize_trans_storage_IX+100));
INNN(3,1)=Index1;
[~,Index1] = min(abs(B_utilize_trans_storage_IX+50));
INNN(4,1)=Index1;
[~,Index1] = min(abs(B_utilize_trans_storage_IX-0));
INNN(5,1)=Index1;
[~,Index1] = min(abs(B_utilize_trans_storage_IX-50));
INNN(6,1)=Index1;
[~,Index1] = min(abs(B_utilize_trans_storage_IX-100));
INNN(7,1)=Index1;
INNN(8,1)=size(B_utilize_trans_storage_IX,1);

carbon_price = 100; % $/t CO2
for i = 1:7
    table(i,1) = INNN(i+1,1)-INNN(i,1); % Number
    table(i,2) = sum(CP_ALL_IX2(INNN(i,1)+1:INNN(i+1,1),1)); % CP (TW)
    table(i,3) = sum(ph_plant_trans_storage_all(INNN(i,1)+1:INNN(i+1,1),1))./sum(CP_ALL_IX2(INNN(i,1)+1:INNN(i+1,1),1))/8760*100; % CF (%)
%     table(i,4) = sum(ph_plant_trans_storage_all(INNN(i,1)+1:INNN(i+1,1),1)); % Power potential without PUE (TWh yr-1)
%     table(i,5) = sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1:8)))/10^3; % Initial investment ((billion $))	
%     table(i,6) = sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),9)))/10^3; % O&M ((billion $))
%     table(i,7) = (sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),10)))+sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),11))))/10^3; % Transmission and storage cost  ((billion $))	
%     table(i,8) = -sum(CO2_year_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1).*degrat40yr_offshorewind2(4,INNN(i,1)+1:INNN(i+1,1))')/1000*10^9*carbon_price/10^9+sum(phhall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1).*EP(INNN(i,1)+1:INNN(i+1,1),1)); % Revenue (billion $)
    table(i,4) = sum(ph_plant_trans_storage(INNN(i,1)+1:INNN(i+1,1),1)); % available Power potential (TWh yr-1)
    table(i,5) = sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1:8)))/10^3/40; % Initial investment ((billion $ yr-1))	
    table(i,6) = sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),9)))/10^3/40; % O&M ((billion $ yr-1))
    table(i,7) = (sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),10)))+sum(sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),11))))/10^3/40; % Transmission and storage cost  ((billion $ yr-1))	
    table(i,8) = -sum(CO2_year_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1).*degrat40yr_offshorewind(4,INNN(i,1)+1:INNN(i+1,1))')/1000*10^9*carbon_price/10^9/40+sum(phhall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1).*EP(INNN(i,1)+1:INNN(i+1,1),1))/40; % Revenue (billion $)
    table(i,9) = sum(sum((sum(lcoeeall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1:11),2))))./sum(sum(phhall_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1)))/1000; % LCOE USD2019/kWh 
    table(i,10) = -sum(CO2_year_utilize_trans_storage_IX2(INNN(i,1)+1:INNN(i+1,1),1))/1000; % CO2 abatement (Gt CO2 y-1)
end
    table(8,1) = INNN(8,1)-INNN(1,1); % Number
    table(8,2) = sum(CP_ALL_IX2(:,1)); % CP (TW)
    table(8,3) = sum(ph_plant_trans_storage_all(:,1))./sum(CP_ALL_IX2(:,1))/8760*100; % CF (%)
%     table(8,4) = sum(ph_plant_trans_storage_all(:,1)); % Power potential without PUE (TWh yr-1)
%     table(8,5) = sum(sum(lcoeeall_utilize_trans_storage_IX2(:,1:8)))/10^3; % Initial investment ((billion $))	
%     table(8,6) = sum(sum(lcoeeall_utilize_trans_storage_IX2(:,9)))/10^3; % O&M ((billion $))
%     table(8,7) = (sum(sum(lcoeeall_utilize_trans_storage_IX2(:,10))))/10^3+(sum(sum(lcoeeall_utilize_trans_storage_IX2(:,11))))/10^3; % Transmission and storage cost  ((billion $))	
%     table(8,8) = -sum(CO2_year_utilize_trans_storage_IX2(:,1).*degrat40yr_offshorewind2(4,:)')/1000*10^9*carbon_price/10^9+sum(phhall_utilize_trans_storage_IX2(:,1).*EP); % Revenue (billion $)
    table(8,4) = sum(ph_plant_trans_storage(:,1)); % Power potential (TWh yr-1)
    table(8,5) = sum(sum(lcoeeall_utilize_trans_storage_IX2(:,1:8)))/10^3/40; % Initial investment ((billion $))	
    table(8,6) = sum(sum(lcoeeall_utilize_trans_storage_IX2(:,9)))/10^3/40; % O&M ((billion $))
    table(8,7) = (sum(sum(lcoeeall_utilize_trans_storage_IX2(:,10))))/10^3/40+(sum(sum(lcoeeall_utilize_trans_storage_IX2(:,11))))/10^3/40; % Transmission and storage cost  ((billion $))	
    table(8,8) = -sum(CO2_year_utilize_trans_storage_IX2(:,1).*degrat40yr_offshorewind(4,:)')/1000*10^9*carbon_price/10^9/40+sum(phhall_utilize_trans_storage_IX2(:,1).*EP)/40; % Revenue (billion $)
    table(8,9) = sum(sum((sum(lcoeeall_utilize_trans_storage_IX2(:,1:10),2))))./sum(sum(phhall_utilize_trans_storage_IX2))/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
    table(8,10) = -sum(CO2_year_utilize_trans_storage_IX2(:,1))/1000; % CO2 abatement (Gt CO2 y-1)

%%
for i = 1:size(B_utilize_trans_storage,1)
    CP11122(i,1) = CP(IX(i,1),1);
    CO2all_c_utilize_trans_storage_IX2(i,1) = CO2all_c_utilize_trans_storage(IX(i,1),1);
    B_utilize_trans_storage_IX2(i,1) = B_utilize_trans_storage(IX(i,1),1);
    storage_year_plant_IX2(i,1) = storage_year_plant(IX(i,1),1);
end
CO2all_c2_utilize_trans_storage_IX2=cumsum(CO2all_c_utilize_trans_storage_IX2);

CP_cumsum_all=cumsum(CP11122);
[m,n] =find(storage_year_plant_IX2==0);
Index_Ie(1,1) = m(end)+1; % The position to start the storage
for i = 3:3:24 % 5:5:25
    [~,Index1] = min(abs(CP_cumsum_all-i));
    Index(1,i/3) = Index1;
%     Index(1,i/5) = Index1;
end
CO2_reduction_utilize_trans_storage = -CO2all_c2_utilize_trans_storage_IX2./10^3;
CO2_reduction5_utilize_trans_storage = CO2_reduction_utilize_trans_storage(Index);
cost_CO2_reduction5_utilize_trans_storage = B_utilize_trans_storage_IX2(Index);


CO2all2_c22_utilize_trans_storage_IX2(2,1) = 0;
CO2all2_c22_utilize_trans_storage_IX2(2,1) =  CO2all_c2_utilize_trans_storage_IX2(1);

B22_utilize_trans_storage_IX2(1,1) =  B_utilize_trans_storage_IX2(1);
B22_utilize_trans_storage_IX2(2,1) =  B_utilize_trans_storage_IX2(1);

for i = 2:size(CO2all_c2_utilize_trans_storage_IX2,1)
    CO2all2_c22_utilize_trans_storage_IX2(i*2-1,1) =  CO2all_c2_utilize_trans_storage_IX2(i-1,1);
    CO2all2_c22_utilize_trans_storage_IX2(i*2,1) =  CO2all_c2_utilize_trans_storage_IX2(i,1);
    B22_utilize_trans_storage_IX2(i*2-1,1) =  B_utilize_trans_storage_IX2(i,1);
    B22_utilize_trans_storage_IX2(i*2,1) =  B_utilize_trans_storage_IX2(i,1);
end  

% PV
[B_utilize_trans_storage_IX_pv,IX_pv]=sort(B_utilize_trans_storage_pv);
for i = 1:size(CO2all_c_utilize_trans_storage_pv,1)
    CO2all_c_utilize_trans_storage_IX2_pv(i,1) = CO2all_c_utilize_trans_storage_pv(IX_pv(i,1),1);
    B_utilize_trans_storage_IX2_pv(i,1) = B_utilize_trans_storage_pv(IX_pv(i,1),1);
end
CO2all_c2_utilize_trans_storage_IX2_pv=cumsum(CO2all_c_utilize_trans_storage_IX2_pv);
CO2all2_c22_utilize_trans_storage_IX2_pv(2,1) = 0;
CO2all2_c22_utilize_trans_storage_IX2_pv(2,1) =  CO2all_c2_utilize_trans_storage_IX2_pv(1);
B22_utilize_trans_storage_IX2_pv(1,1) =  B_utilize_trans_storage_IX2_pv(1);
B22_utilize_trans_storage_IX2_pv(2,1) =  B_utilize_trans_storage_IX2_pv(1);
for i = 2:size(CO2all_c2_utilize_trans_storage_IX2_pv,1)
    CO2all2_c22_utilize_trans_storage_IX2_pv(i*2-1,1) =  CO2all_c2_utilize_trans_storage_IX2_pv(i-1,1);
    CO2all2_c22_utilize_trans_storage_IX2_pv(i*2,1) =  CO2all_c2_utilize_trans_storage_IX2_pv(i,1);
    B22_utilize_trans_storage_IX2_pv(i*2-1,1) =  B_utilize_trans_storage_IX2_pv(i,1);
    B22_utilize_trans_storage_IX2_pv(i*2,1) =  B_utilize_trans_storage_IX2_pv(i,1);
end  
% onshore wind
[B_utilize_trans_storage_IX_onshorewind,IX_onshorewind]=sort(B_utilize_trans_storage_onshorewind);
for i = 1:size(CO2all_c_utilize_trans_storage_onshorewind,1)
    CO2all_c_utilize_trans_storage_IX2_onshorewind(i,1) = CO2all_c_utilize_trans_storage_onshorewind(IX_onshorewind(i,1),1);
    B_utilize_trans_storage_IX2_onshorewind(i,1) = B_utilize_trans_storage_onshorewind(IX_onshorewind(i,1),1);
end
CO2all_c2_utilize_trans_storage_IX2_onshorewind=cumsum(CO2all_c_utilize_trans_storage_IX2_onshorewind);
CO2all2_c22_utilize_trans_storage_IX2_onshorewind(2,1) = 0;
CO2all2_c22_utilize_trans_storage_IX2_onshorewind(2,1) =  CO2all_c2_utilize_trans_storage_IX2_onshorewind(1);
B22_utilize_trans_storage_IX2_onshorewind(1,1) =  B_utilize_trans_storage_IX2_onshorewind(1);
B22_utilize_trans_storage_IX2_onshorewind(2,1) =  B_utilize_trans_storage_IX2_onshorewind(1);
for i = 2:size(CO2all_c2_utilize_trans_storage_IX2_onshorewind,1)
    CO2all2_c22_utilize_trans_storage_IX2_onshorewind(i*2-1,1) =  CO2all_c2_utilize_trans_storage_IX2_onshorewind(i-1,1);
    CO2all2_c22_utilize_trans_storage_IX2_onshorewind(i*2,1) =  CO2all_c2_utilize_trans_storage_IX2_onshorewind(i,1);
    B22_utilize_trans_storage_IX2_onshorewind(i*2-1,1) =  B_utilize_trans_storage_IX2_onshorewind(i,1);
    B22_utilize_trans_storage_IX2_onshorewind(i*2,1) =  B_utilize_trans_storage_IX2_onshorewind(i,1);
end  

% offshorewind
[B_utilize_trans_storage_IX_offshorewind,IX_offshorewind]=sort(B_utilize_trans_storage_offshorewind);
for i = 1:size(CO2all_c_utilize_trans_storage_offshorewind,1)
    CO2all_c_utilize_trans_storage_IX2_offshorewind(i,1) = CO2all_c_utilize_trans_storage_offshorewind(IX_offshorewind(i,1),1);
    B_utilize_trans_storage_IX2_offshorewind(i,1) = B_utilize_trans_storage_offshorewind(IX_offshorewind(i,1),1);
end
CO2all_c2_utilize_trans_storage_IX2_offshorewind=cumsum(CO2all_c_utilize_trans_storage_IX2_offshorewind);
CO2all2_c22_utilize_trans_storage_IX2_offshorewind(2,1) = 0;
CO2all2_c22_utilize_trans_storage_IX2_offshorewind(2,1) =  CO2all_c2_utilize_trans_storage_IX2_offshorewind(1);
B22_utilize_trans_storage_IX2_offshorewind(1,1) =  B_utilize_trans_storage_IX2_offshorewind(1);
B22_utilize_trans_storage_IX2_offshorewind(2,1) =  B_utilize_trans_storage_IX2_offshorewind(1);
for i = 2:size(CO2all_c2_utilize_trans_storage_IX2_offshorewind,1)
    CO2all2_c22_utilize_trans_storage_IX2_offshorewind(i*2-1,1) =  CO2all_c2_utilize_trans_storage_IX2_offshorewind(i-1,1);
    CO2all2_c22_utilize_trans_storage_IX2_offshorewind(i*2,1) =  CO2all_c2_utilize_trans_storage_IX2_offshorewind(i,1);
    B22_utilize_trans_storage_IX2_offshorewind(i*2-1,1) =  B_utilize_trans_storage_IX2_offshorewind(i,1);
    B22_utilize_trans_storage_IX2_offshorewind(i*2,1) =  B_utilize_trans_storage_IX2_offshorewind(i,1);
end  

y = [-20:20:200]';
[~,Index1] = min(abs(B_utilize_trans_storage_IX2-0));
CP_cumsum_all(Index1,1)
%%
save('H:\Global PV and wind\ANS\B22_utilize_trans_storage_IX22_8_2040_2s_2050_test6xzsuper.mat','B22_utilize_trans_storage_IX2', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\CO2all2_c22_utilize_trans_storage_IX22_8_2040_2s_2050_test6xzsuper.mat','CO2all2_c22_utilize_trans_storage_IX2', '-v7.3')  % billion kWh / year
CO2_cum = -CO2all_c2_utilize_trans_storage_IX2/10^3;
MAC = B_utilize_trans_storage_IX2;
CO2 = -CO2all_c_utilize_trans_storage_IX2/10^3;
save('H:\Global PV and wind\ANS\CO2_8_2040_2s_2050_test6xzsuper.mat','CO2'); %
save('H:\Global PV and wind\ANS\MAC_8_2040_2s_2050_test6xzsuper.mat','MAC'); %
save('H:\Global PV and wind\ANS\CO2_cum_8_2040_2s_2050_test6xzsuper.mat','CO2_cum'); %

save('H:\Global PV and wind\ANS\optpowerunit_IX_8_2040_2s_2050_test6xzsuper.mat','optpowerunit_IX'); %
save('H:\Global PV and wind\ANS\powerunit_country_IX_IX_8_2040_2s_2050_test6xzsuper.mat','powerunit_country_IX_IX'); %
save('H:\Global PV and wind\ANS\emissionfactor_8_2040_2s_2050_test6xzsuper.mat','emissionfactor'); %

phhall_all2 = phhall_utilize_trans_storage./degrat40yr_PV(4,:)';
% B_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:11),2)-phhall_utilize_trans_storage(:,1).*EP*10^3)./(-CO2_all_utilize_trans_storage); %USD/t CO2
save('H:\Global PV and wind\ANS\phhall_all2_8_2040_2s_2050_test6xzsuper.mat','phhall_all2', '-v7.3')  % TWh / year
save('H:\Global PV and wind\ANS\EP_8_2040_2s_2050_test6xzsuper.mat','EP', '-v7.3')  % billion kWh / year
save('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_8_2040_2s_2050_test6xzsuper.mat','CO2_year_utilize_trans_storage', '-v7.3')  % Mt CO2/year

[sum(lcoeeall_utilize_trans_storage_type) sum(-phhall_utilize_trans_storage(:,1).*EP*10^3)]

lcoeeall_component = [lcoeeall_utilize_trans_storage_type, -phhall_utilize_trans_storage(:,1).*EP*10^3];
save('H:\Global PV and wind\ANS\lcoeeall_component_8_2040_2s_2050_test6xzsuper.mat','lcoeeall_component', '-v7.3')  %


degrat = degrat40yr_PV(4,:)';
save('H:\Global PV and wind\ANS\degrat_8_2040_2s_2050_test6xzsuper.mat','degrat', '-v7.3')  %
