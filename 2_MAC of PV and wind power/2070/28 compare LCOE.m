tic
clear;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_sto_interUHVxz.mat')  % 
LCOEE(:,8) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_stoxz.mat')  % 
LCOEE(:,7) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_sto_interUHVxz.mat')  % 
LCOEE(:,6) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_interUHVxz.mat')  % 
LCOEE(:,5) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHVxz.mat')  % 
LCOEE(:,4) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_stoxz.mat')  % 
LCOEE(:,3) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_interUHVxz.mat')  % 
LCOEE(:,2) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testtxz.mat')  % 
LCOEE(:,1) = LCOEE_all_utilize_trans_storage;

LCOEEmin = min(LCOEE')';
choo = zeros(size(LCOEE,1),1);
for i = 1:8
    [m,n]=find(LCOEE(:,i)-LCOEEmin==0);
    choo(m)=i;
end
save('H:\global-PV-wind\ANS\choo_type_8_2070xz.mat','choo'); 
% 1.no UHV,sto and international UHV
% 2.no sto and international UHV
% 3.no international UHV
% 4.all

%%
tic
clear;
load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% MW
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070.mat')  % power use efficiency with global electric transport and storage
load('H:\global-PV-wind\ANS\storage_max_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070.mat')   % Maximum charging power, TWh/h  
load('H:\global-PV-wind\ANS\storage_year_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070.mat')  % Average annual electricity storage, TWh/year
numpowerunit = size(cost_trans_IX,1);
cost_trans(:,8) = cost_trans_IX;
CP_trans(:,8) = CP_trans_IX;
utilize_ratio(:,8) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,8) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,8) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% MW
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')  % power use efficiency with global electric transport and storage
load('H:\global-PV-wind\ANS\storage_max_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')   % Maximum charging power, TWh/h  
load('H:\global-PV-wind\ANS\storage_year_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')  % Average annual electricity storage, TWh/year
cost_trans(:,7) = cost_trans_IX;
CP_trans(:,7) = CP_trans_IX;
utilize_ratio(:,7) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,7) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,7) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% MW
cost_trans_IX2 = cost_trans_IX;
CP_trans_IX2 = CP_trans_IX;
load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% MW
cost_trans_IX = cost_trans_IX2 - cost_trans_IX;
CP_trans_IX = CP_trans_IX2 - CP_trans_IX;
clear cost_trans_IX2
clear CP_trans_IX2
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_nodomUHV2.mat')  % power use efficiency with global electric transport and storage
load('H:\global-PV-wind\ANS\storage_max_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070.mat')   % Maximum charging power, TWh/h  
load('H:\global-PV-wind\ANS\storage_year_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070.mat')  % Average annual electricity storage, TWh/year
cost_trans(:,6) = cost_trans_IX;
CP_trans(:,6) = CP_trans_IX;
utilize_ratio(:,6) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,6) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,6) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% MW
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_nosto2.mat')  % power use efficiency with global electric transport and storage
cost_trans(:,5) = cost_trans_IX;
CP_trans(:,5) = CP_trans_IX;
utilize_ratio(:,5) =utilize_ratio2060_UHV_storage_inter;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX

load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% MW
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_nointerUHVsto2.mat','utilize_ratio2060_UHV_storage_inter')
storage_max_plant_UHV_storage_inter = zeros(numpowerunit,1);
storage_inter_year_plant_UHV_storage = zeros(numpowerunit,1);
cost_trans(:,4) = cost_trans_IX;
CP_trans(:,4) = CP_trans_IX;
utilize_ratio(:,4) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,4) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,4) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

cost_trans_IX = zeros(numpowerunit,1);
CP_trans_IX = zeros(numpowerunit,1);
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_noUHV2.mat')  % power use efficiency with global electric transport and storage
load('H:\global-PV-wind\ANS\storage_max_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')   % Maximum charging power, TWh/h  
load('H:\global-PV-wind\ANS\storage_year_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')  % Average annual electricity storage, TWh/year
cost_trans(:,3) = cost_trans_IX;
CP_trans(:,3) = CP_trans_IX;
utilize_ratio(:,3) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,3) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,3) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070.mat')% MW
cost_trans_IX2 = cost_trans_IX;
CP_trans_IX2 = CP_trans_IX;
load('H:\global-PV-wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV.mat')% MW
cost_trans_IX = cost_trans_IX2 - cost_trans_IX;
CP_trans_IX = CP_trans_IX2 - CP_trans_IX;
clear cost_trans_IX2
clear CP_trans_IX2
load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_nodomUHVsto2.mat')  % power use efficiency with global electric transport and storage
storage_max_plant_UHV_storage_inter = zeros(numpowerunit,1);
storage_inter_year_plant_UHV_storage = zeros(numpowerunit,1);
cost_trans(:,2) = cost_trans_IX;
CP_trans(:,2) = CP_trans_IX;
utilize_ratio(:,2) =utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,2) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,2) = storage_inter_year_plant_UHV_storage;
clear utilize_ratio2060_UHV_storage_inter
clear cost_trans_IX
clear CP_trans_IX
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all_battery_mechanical_pro2_8_2070_noUHVsto2.mat','utilize_ratio2060_UHV_storage_inter')
utilize_ratio(:,1) =utilize_ratio2060_UHV_storage_inter;
clear utilize_ratio2060_UHV_storage_inter

load('H:\global-PV-wind\ANS\choo_type_8_2070xz.mat'); % choo
for i = 1:8
    [m,n]=find(choo==i);
    cost_trans_IX(m,1) = cost_trans(m,i);
    CP_trans_IX(m,1) = CP_trans(m,i);
    utilize_ratio2060_UHV_storage_inter(m,1) = utilize_ratio(m,i);
    storage_max_plant_UHV_storage_inter(m,1) = storage_max_plant(m,i);
    storage_inter_year_plant_UHV_storage(m,1) = storage_year_plant(m,i);
end
clear utilize_ratio
clear cost_trans
clear CP_trans
clear storage_max_plant
clear storage_year_plant

save('H:\global-PV-wind\ANS\utilize_ratio_8_2070_xzxz.mat','utilize_ratio2060_UHV_storage_inter'); 
save('H:\global-PV-wind\ANS\cost_trans_8_2070_xzxz.mat','cost_trans_IX'); 
save('H:\global-PV-wind\ANS\CP_trans_8_2070_xzxz.mat','CP_trans_IX'); 
save('H:\global-PV-wind\ANS\storage_max_plant_8_2070_xzxz.mat','storage_max_plant_UHV_storage_inter'); 
save('H:\global-PV-wind\ANS\storage_year_plant_8_2070_xzxz.mat','storage_inter_year_plant_UHV_storage'); 




%%
tic
clear;
load('H:\global-PV-wind\Data\r_module_hardwarecost_pv.mat');
% Module占代码中Moduleprice2的比例
% 硬件成本占代码中Moduleprice2的比例
r_module_hardwarecost_ons = 0.704960836;
r_fixednomodule_hardwarecost_ons = 0.7493-r_module_hardwarecost_ons;

load('H:\global-PV-wind\ANS\cho_county_all_pro2_8_2070.mat'); %
load('H:\global-PV-wind\ANS\countryy_IX_county_all_pro2_8_2070.mat')
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];
load('H:\global-PV-wind\ANS\cost_trans_8_2070_xzxz.mat')% million $和optpowerunit_IX顺序一样
load('H:\global-PV-wind\ANS\CP_trans_8_2070_xzxz.mat')% MW

load('H:\global-PV-wind\ANS\utilize_ratio_8_2070_xzxz.mat')  % power use efficiency with global electric transport and storage
utilize_ratio2060 =utilize_ratio2060_UHV_storage_inter;
utilize_ratio2060(find(isnan(utilize_ratio2060)==1))=0;
clear utilize_ratio2060_UHV_storage_inter
numpowerunit = size(utilize_ratio2060,1);

load('H:\global-PV-wind\ANS\storage_max_plant_8_2070_xzxz.mat')   % Maximum charging power, TWh/h
load('H:\global-PV-wind\ANS\storage_year_plant_8_2070_xzxz.mat')  % Average annual electricity storage, TWh/year
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

load('H:\global-PV-wind\Data\discount_country.mat')
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
Ip_storageall=zeros(numpowerunit,1);
Ie_storageall=zeros(numpowerunit,1);
Ip_mechanical=zeros(numpowerunit,1);
cost_trans=zeros(numpowerunit,1);
CP_trans=zeros(numpowerunit,1);
Cost_m_ini1 = zeros(numpowerunit,1); %
Cost_m_OM1 = zeros(numpowerunit,1); %
utilize_ratio2060_2=zeros(numpowerunit,1);
Cost_mechanical1(mmm,1) = Cost_mechanical(mmm,1);
Cost_m_ini1(mmm,1) = Cost_m_ini(mmm,1); %million $
Cost_m_OM1(mmm,1) = Cost_m_OM(mmm,1); %million $
Ip_mechanical(mmm,1) = Ip(mmm,1); %kW
Ip_storageall(mmm,1) = Ip(mmm,1); %kW
Ie_storageall(mmm,1) = Ie(mmm,1); %kWh
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1);
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);


% load('H:\global-PV-wind\ANS\utilize_ratio2060_UHV_storage_inter_county_all.mat')  % power use efficiency with global electric transport and storage
% utilize_ratio2060 =utilize_ratio2060_UHV_storage_inter;
% clear utilize_ratio2060_UHV_storage_inter
%
% load('H:\global-PV-wind\ANS\storage_max_plant_UHV_storage_inter_county_all.mat')   % Maximum charging power, TWh/h
% load('H:\global-PV-wind\ANS\storage_year_plant_UHV_storage_inter_county_all.mat')  % Average annual electricity storage, TWh/year
% storage_max_plant = storage_max_plant_UHV_storage_inter;
% storage_year_plant = storage_inter_year_plant_UHV_storage;
% clear storage_max_plant_UHV_storage_inter
% clear storage_year_plant_UHV_storage_inter
% load('H:\global-PV-wind\ANS\cost_trans_IX_all_battery_mechanical_all.mat')% million $和optpowerunit_IX顺序一样
% load('H:\global-PV-wind\ANS\CP_trans_IX_all_battery_mechanical_all.mat')% MW

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
load('H:\global-PV-wind\ANS\unitmin_global_IX_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070xz.mat');
unitmin_CHA = unitmin;
unitmin(unitmin==1) = 11;
unitmin(unitmin==2) = 7;
unitmin(unitmin==3) = 8;
unitmin(unitmin==4) = 9;
unitmin(unitmin==5) = 10;
unitmin(unitmin==6) = 11;
rIp=[595.73,440,374.45,374.45-(374.45-327.22)/2, 327.22,327.22-(327.22-280.78)/2, 280.78, 280.78-(280.78-234.343)/2, 234.343, 234.343-(280.78-234.343)/2,234.343-(280.78-234.343)]./595.73;
rIe=[345,242,198,186, 174, 161, 149, 149-(149-124)/2, 124, 124-(149-124)/2, 124-(149-124)]./345;
for i = 1:11
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
Ip_storageall(mmm,1) = Ip(mmm,1); %kW
Ie_storageall(mmm,1) = Ie(mmm,1); %kWh
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1);
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);
% clear Inc
clear Cost_storage
Cost_storage=Cost_storage1;
clear Cost_storage1

clear Ip_storage

Cost_sto_ini = Cost_m_ini1;
Cost_sto_OM = Cost_m_OM1;
[m,nnn]=find(cho==-1);
Cost_sto_ini(m) = Cost_storage1_ini(m); % million $/y
Cost_sto_OM(m) = Cost_storage1_OM1(m); % million $/y

%
[mmm,nnn]=find(cho==0);
utilize_ratio2060_2(mmm,1) = utilize_ratio2060(mmm,1);
cost_trans(mmm,1) = cost_trans_IX(mmm,1);
CP_trans(mmm,1) = CP_trans_IX(mmm,1);
Ip_storageall(mmm,1) = Ip(mmm,1); %kW
Ie_storageall(mmm,1) = Ie(mmm,1); %kWh
clear utilize_ratio2060
clear CP_trans_IX
clear cost_trans_IX
% utilize_ratio2060_trans_plant_alone = utilize_ratio2060_2;
% utilize_ratio2060_trans_storage_plant_alone = utilize_ratio2060_2;
utilize_ratio2060 = utilize_ratio2060_2;
clear utilize_ratio2060_2

%
LR_PV2 = 0.18; %0.18; % 0.37; %  0.37 PV module
CP0_PV = sum(sum(CP_trans))/10; % MW
learnprice_trans(1)=1;
for i = 2:10
    [m,n]=find(unitmin_CHA<=i-1);
    learnprice_trans(i,1) = ((sum(CP_trans(m))+CP0_PV)./CP0_PV).^(log(1- LR_PV2)/log(2));
end
learnprice_trans(11,1) = ((sum(CP_trans)+CP0_PV)./CP0_PV).^(log(1- LR_PV2)/log(2));
learnprice_trans(learnprice_trans>1)=1;


CO2_C=0.2727;
lifetime_power=25;
lifetime_inverter=10; % renewed per 10 years
OMratio_PV = 0.01;
OMratio_wind = 0.03;
discountinverter=zeros(192,1);

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
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind;

% %
% load('H:\global-PV-wind\ANS\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
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

load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
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

% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\global-PV-wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\global-PV-wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\global-PV-wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

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


%
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

%
load('H:\global-PV-wind\Data\Stocks2018_PumpedHydro.mat'); % "A" grade (GWh)	"B" grade (GWh)	"C" grade (GWh)	"D" grade (GWh)	"E" grade (GWh)
Stocks2018_PumpedHydro = Stocks2018_PumpedHydro/1000; % TWh

for ccc = 1:192
    [mmm,nnn]=find(cho==1 & powerunit_country_IX_IX==ccc);
    aaaa(ccc,1) = sum(Ie_storageall(mmm))/10^9; % TWh
end
[ccou,n] = find(Stocks2018_PumpedHydro(:,1)-aaaa<0);
% Stocks2018_PumpedHydro(ccou,1)-aaaa(ccou,1)
choo_stoxz = cho*0;
for ccc = 1:size(ccou,1)
    [mmm,nnn]=find(cho==1 & powerunit_country_IX_IX==ccou(ccc));
    a = cumsum(Ie_storageall(mmm))/10^9; % TWh
    b = Stocks2018_PumpedHydro(ccou(ccc),1);
    [mmm2,nnn2]=find(a>b);
    cho(mmm(mmm2)) = -1;
end
save('H:\global-PV-wind\ANS\cho_county_all_pro2_8_2070_xz.mat','cho'); %
% % 0是无需储能，1是机械储能，-1是化学储能
% save('H:\global-PV-wind\ANS\Ie_storageall_2070CCCCCC.mat','Ie_storageall'); %
% save('H:\global-PV-wind\ANS\Ip_storageall_2070CCCCCC.mat','Ip_storageall'); %

