tic
clear;
load('H:\Global PV and wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_testt_UHVxz.mat')  % 
LCOEE(:,2) = LCOEE_all_utilize_trans_storage;
load('H:\Global PV and wind\ANS\LCOE_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_testtxz.mat')  % 
LCOEE(:,1) = LCOEE_all_utilize_trans_storage;
clear LCOEE_all_utilize_trans_storage

LCOEEmin = min(LCOEE')';
choo = zeros(size(LCOEE,1),1);
for i = 1:2
    [m,n]=find(LCOEE(:,i)-LCOEEmin==0);
    choo(m)=i;
end
save('H:\Global PV and wind\ANS\choo_type_8_2040xz_CaseB.mat','choo'); 
% 1.no UHV,sto and international UHV
% 2.no sto and international UHV
% 3.no international UHV
% 4.all

%%
tic
clear;
load('H:\Global PV and wind\ANS\cost_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV_nodomSto.mat')% million $和optpowerunit_IX顺序一样
load('H:\Global PV and wind\ANS\CP_trans_IX_battery_mechanical_all_pro2_8_2070_nointerUHV_nodomSto.mat')% MW
load('H:\Global PV and wind\ANS\utilize_ratio2060_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV_nodomSto.mat')  % power use efficiency with global electric transport and storage
load('H:\Global PV and wind\ANS\storage_max_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV_nodomSto.mat')   % Maximum charging power, TWh/h  
load('H:\Global PV and wind\ANS\storage_year_plant_UHV_storage_inter_county_battery_mechanical_all_pro2_8_2070_nointerUHV_nodomSto.mat')  % Average annual electricity storage, TWh/year
cost_trans(:,2) = cost_trans_IX;
CP_trans(:,2) = CP_trans_IX;
utilize_ratio(:,2) = utilize_ratio2060_UHV_storage_inter;
storage_max_plant(:,2) = storage_max_plant_UHV_storage_inter;
storage_year_plant(:,2) = storage_inter_year_plant_UHV_storage;

load('H:\Global PV and wind\ANS\utilize_ratio2060_county_battery_mechanical_all_pro2_8_2070.mat')  % power use efficiency with global electric transport and storage
utilize_ratio(:,1) = utilize_ratio2060;
clear utilize_ratio2060
clear cost_trans_IX
clear CP_trans_IX
clear utilize_ratio2060_UHV_storage_inter
clear storage_max_plant_UHV_storage_inter
clear storage_inter_year_plant_UHV_storage

load('H:\Global PV and wind\ANS\choo_type_8_2040xz_CaseB.mat'); % choo
for i = 1:2
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
clear storage_max_planta
clear storage_year_plant

save('H:\Global PV and wind\ANS\utilize_ratio_8_2040_xzxz_CaseB.mat','utilize_ratio2060_UHV_storage_inter'); 
save('H:\Global PV and wind\ANS\cost_trans_8_2040_xzxz_CaseB.mat','cost_trans_IX'); 
save('H:\Global PV and wind\ANS\CP_trans_8_2040_xzxz_CaseB.mat','CP_trans_IX'); 
save('H:\Global PV and wind\ANS\storage_max_plant_8_2040_xzxz_CaseB.mat','storage_max_plant_UHV_storage_inter'); 
save('H:\Global PV and wind\ANS\storage_year_plant_8_2040_xzxz_CaseB.mat','storage_inter_year_plant_UHV_storage'); 


