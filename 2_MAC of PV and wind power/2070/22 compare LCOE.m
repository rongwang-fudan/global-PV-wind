tic
clear;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_sto_interUHV.mat')  % 
LCOEE(:,8) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_sto.mat')  % 
LCOEE(:,7) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_sto_interUHV.mat')  % 
LCOEE(:,6) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV_interUHV.mat')  % 
LCOEE(:,5) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_UHV.mat')  % 
LCOEE(:,4) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_sto.mat')  % 
LCOEE(:,3) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt_interUHV.mat')  % 
LCOEE(:,2) = LCOEE_all_utilize_trans_storage;
load('H:\global-PV-wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_testt.mat')  % 
LCOEE(:,1) = LCOEE_all_utilize_trans_storage;

LCOEEmin = min(LCOEE')';
choo = zeros(size(LCOEE,1),1);
for i = 1:8
    [m,n]=find(LCOEE(:,i)-LCOEEmin==0);
    choo(m)=i;
end
save('H:\global-PV-wind\ANS\choo_type_8_2070.mat','choo'); 
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

load('H:\global-PV-wind\ANS\choo_type_8_2070.mat'); % choo
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

save('H:\global-PV-wind\ANS\utilize_ratio_8_2070_xz.mat','utilize_ratio2060_UHV_storage_inter'); 
save('H:\global-PV-wind\ANS\cost_trans_8_2070_xz.mat','cost_trans_IX'); 
save('H:\global-PV-wind\ANS\CP_trans_8_2070_xz.mat','CP_trans_IX'); 
save('H:\global-PV-wind\ANS\storage_max_plant_8_2070_xz.mat','storage_max_plant_UHV_storage_inter'); 
save('H:\global-PV-wind\ANS\storage_year_plant_8_2070_xz.mat','storage_inter_year_plant_UHV_storage'); 


