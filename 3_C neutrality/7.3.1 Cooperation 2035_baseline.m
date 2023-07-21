%% 
tic
clear;

%% 筛选出2040年电厂
load('H:\Global PV and wind\ANS\powerunit_country_IX_IX_8_2040_2s_2060s_test6xz_baseline.mat'); % powerunit_country_IX_IX
% load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_7_2040_2s_2020s_test6.mat')  % B_utilize_trans_storage
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2060s_test6xz_baseline.mat')  % B_utilize_trans_storage
B_utilize_trans_storage_nodiffuse = B_utilize_trans_storage;
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2020s_test6xz.mat')  % B_utilize_trans_storage
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_192countrya_2040_8.mat')  % B_utilize_trans_storage_192
load('H:\Global PV and wind\ANS\EP_8_2040_2s_2060s_test6xz_baseline.mat')  % EP，USD/kWh
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_8_2040_2s_2060s_test6xz_baseline.mat')  % CO2_year_utilize_trans_storage, Mt CO2/year
load('H:\Global PV and wind\ANS\phhall_all2_8_2040_2s_2060s_test6xz_baseline.mat')
load('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat'); %%
load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat')
load('H:\Global PV and wind\ANS\optpowerunit_IX_8_2040_2s_2060s_test6xz_baseline.mat'); % optpowerunit_IX
% Plant_ID_IX
% 但这个不是按照原始LCOE排序的
optpowerunit_IX_IX = optpowerunit_IX(Plant_ID_IX,:);
powerunit_country_IX_IX = powerunit_country_IX_IX(Plant_ID_IX,:);
B_utilize_trans_storage_nodiffuse = B_utilize_trans_storage_nodiffuse(Plant_ID_IX,:);
B_utilize_trans_storage = B_utilize_trans_storage(Plant_ID_IX,:);
B_utilize_trans_storage_192 = B_utilize_trans_storage_192(Plant_ID_IX,:);
EP = EP(Plant_ID_IX,:);
CO2_year_utilize_trans_storage = CO2_year_utilize_trans_storage(Plant_ID_IX,:);
phhall_all2 = phhall_all2(Plant_ID_IX,:);
unitmin = unitmin(Plant_ID_IX,:);
load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_8_2070.mat')
etrans = etrans_cou1_num1(Plant_ID_IX,:,:);
clear Plant_ID_IX
clear etrans_cou1_num1
EF = CO2_year_utilize_trans_storage./phhall_all2; % 单位发电量的减排量

etrans1 = reshape(sum(etrans,2),[size(unitmin,1),192]);
clear etrans

load('H:\Global PV and wind\ANS\Ph_PVwind5_ef_cou_p_2040_8_2s_6_2a.mat')  % TWh/yr, Ph_PVwind5_ef_cou
load('H:\Global PV and wind\ANS\m_2040_8_2s_6_2a.mat')  % m,由于没有带来更大的收益在2021-2030年被剔除的电厂
Ph_PVwind5_ef_cou_p(m,:,:) = [];
Ph_PVwind5_ef_cou_p2040 = Ph_PVwind5_ef_cou_p(:,:,3);
etrans2 = min(etrans1,Ph_PVwind5_ef_cou_p2040);
clear etrans1
Ph_PVwind5_ef_cou_p2040 = etrans2;
r = sum(Ph_PVwind5_ef_cou_p2040,2)./phhall_all2;
r(find(isnan(r)==1))=0;
phhall_all2 = sum(Ph_PVwind5_ef_cou_p2040,2);
CO2_year_utilize_trans_storage = r.*CO2_year_utilize_trans_storage;
etrans_cou1_num = zeros(size(phhall_all2,1),192,192);
for i = 1:size(phhall_all2,1)
    cc = powerunit_country_IX_IX(i);
    etrans_cou1_num(i,cc,:) = Ph_PVwind5_ef_cou_p2040(i,:);
end
clear Ph_PVwind5_ef_cou_p

numpowerunit_all = size(unitmin,1);
[m,n]=find(unitmin<=3);
m21to40 = m;
[m,n]=find(unitmin>3);
optpowerunit_IX_IX(m,:) = [];
powerunit_country_IX_IX(m,:) = [];
B_utilize_trans_storage_nodiffuse(m,:) = [];
B_utilize_trans_storage(m,:) = [];
B_utilize_trans_storage_192(m,:) = [];
EP(m,:) = [];
CO2_year_utilize_trans_storage(m,:) = [];
phhall_all2(m,:) = [];
unitmin(m,:) = [];
etrans_cou1_num(m,:,:) = [];
numpowerunit = size(unitmin,1);
etrans_t = phhall_all2;
r_etrans_cou1_num = etrans_cou1_num*0;
for i = 1:numpowerunit
    if phhall_all2(i)~=0
    r_etrans_cou1_num(i,:,:) = etrans_cou1_num(i,:,:)./phhall_all2(i);
    end
end
Inv_nodiffuse = B_utilize_trans_storage_nodiffuse.*(-CO2_year_utilize_trans_storage); % million USD/yr
Inv_nodiffuse(find(isnan(Inv_nodiffuse)==1))=0;
Inv_nodiffuse_total = B_utilize_trans_storage_nodiffuse.*(-CO2_year_utilize_trans_storage)+phhall_all2.*EP*10^3; % million USD/yr
Inv_nodiffuse_total(find(isnan(Inv_nodiffuse_total)==1))=0;
MAC_CCS = 130.85;

%%
load('H:\Global PV and wind\ANS\mineral_limit_dis_2020s_8a.mat')
mineral_limit_20sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2020s_8b.mat')
mineral_limit_20sb = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2030s_8a.mat')
mineral_limit_30sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2030s_8b.mat')
mineral_limit_30sb = mineral_limit;
clear mineral_limit
% mineral_limit, thousand metric tons
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
mineral_limit_20sa(:,7)=0;
mineral_limit_20sb(:,7)=0;
mineral_limit_30sa(:,7)=0;
mineral_limit_30sb(:,7)=0;

load('H:\Global PV and wind\Data\mineral_type.mat');
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% 数值：1是仅PV用，2是仅wind用，3是PV和wind都用得到
mineral_limit_20s2a = mineral_limit_20sa./sum(mineral_limit_20sa);
mineral_limit_20s2a(find(isnan(mineral_limit_20s2a)==1))=0;
mineral_limit_20s2b = mineral_limit_20sb./sum(mineral_limit_20sb);
mineral_limit_20s2b(find(isnan(mineral_limit_20s2b)==1))=0;
mineral_limit_30s2a = mineral_limit_30sa./sum(mineral_limit_30sa);
mineral_limit_30s2a(find(isnan(mineral_limit_30s2a)==1))=0;
mineral_limit_30s2b = mineral_limit_30sb./sum(mineral_limit_30sb);
mineral_limit_30s2b(find(isnan(mineral_limit_30s2b)==1))=0;
[m,n]=find(mineral_type~=2); % PV
r_mineral_20sa(:,1) = sum(mineral_limit_20s2a(:,n),2)./sum(sum(mineral_limit_20s2a(:,n)))
r_mineral_20sb(:,1) = sum(mineral_limit_20s2b(:,n),2)./sum(sum(mineral_limit_20s2b(:,n)))
r_mineral_30sa(:,1) = sum(mineral_limit_30s2a(:,n),2)./sum(sum(mineral_limit_30s2a(:,n)))
r_mineral_30sb(:,1) = sum(mineral_limit_30s2b(:,n),2)./sum(sum(mineral_limit_30s2b(:,n)))
[m,n]=find(mineral_type~=1); % Wind
r_mineral_20sa(:,2) = sum(mineral_limit_20s2a(:,n),2)./sum(sum(mineral_limit_20s2a(:,n)))
r_mineral_20sb(:,2) = sum(mineral_limit_20s2b(:,n),2)./sum(sum(mineral_limit_20s2b(:,n)))
r_mineral_30sa(:,2) = sum(mineral_limit_30s2a(:,n),2)./sum(sum(mineral_limit_30s2a(:,n)))
r_mineral_30sb(:,2) = sum(mineral_limit_30s2b(:,n),2)./sum(sum(mineral_limit_30s2b(:,n)))
r_mineral_pv = [r_mineral_20sa(:,1),r_mineral_20sb(:,1),r_mineral_30sa(:,1),r_mineral_30sb(:,1)];
r_mineral_wind = [r_mineral_20sa(:,2),r_mineral_20sb(:,2),r_mineral_30sa(:,2),r_mineral_30sb(:,2)];

Owncountry_ID2_all = zeros(numpowerunit,8);
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2020s_8a.mat')
% Owncountry_ID2, mineral from own country is 0
[m,n]=find(unitmin==1);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2020s_8b.mat')
[m,n]=find(unitmin==2);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2030s_8a.mat')
[m,n]=find(unitmin==3);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
Owncountry_ID2 = Owncountry_ID2_all;
clear Owncountry_ID2_all

%% 各国自给自足，无需国家之间mineral和power transportion和technology diffussion
power_mineral_cou = zeros(192,2);
EPpower_mineral_cou = zeros(192,1);
EFpower_mineral_cou = zeros(192,1);

plant_use = zeros(size(B_utilize_trans_storage,1),1);
inv = zeros(1,1);
inv_diffuse = zeros(1,1);
power_mineral = zeros(1,1);
asw = zeros(192,1);
nn_c = 1;

r_power_coop = zeros(numpowerunit,1); % (i,j)表示两国合作导致的i电厂的有用电量比例
power_coop = zeros(numpowerunit,1); 
inv_cou = zeros(192,1);
   
inv_total_cou = zeros(192,1);
for i = 1:192
    j = i;
    %
    [m,n]=find(powerunit_country_IX_IX==i & sum(Owncountry_ID2,2)==0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS );% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,1) = r_power_coop(m,1)+r_etrans_cou1_num(m,i,i);
    power_coop(m,1) = power_coop(m,1)+etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = 1;
    %
%     [m,n]=find(powerunit_country_IX_IX==i & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,1);
    r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
    power_coop(m,1) = power_coop(m,1)+r_min*etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = plant_use(m,1)+r_min;
    
    for ttt = 2:3
    [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_pv(i,ttt);
    r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
    power_coop(m,1) = power_coop(m,1)+r_min*etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = plant_use(m,1)+r_min;
    end

    
    %
%     [m,n]=find(powerunit_country_IX_IX==i & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,2);
    r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
    power_coop(m,1) = power_coop(m,1)+r_min*etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = plant_use(m,1)+r_min;
    
    for ttt = 2:3
    [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_wind(i,ttt);
    r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
    power_coop(m,1) = power_coop(m,1)+r_min*etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = plant_use(m,1)+r_min;
    end
    
end
save('H:\Global PV and wind\ANS\power_mineral_cou2035_2040baseline.mat','power_mineral_cou','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2035_2040baseline.mat','inv_total_cou','-v7.3')  % million $,total cost
save('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2035_2040baseline.mat','inv_cou','-v7.3')  % TWh/year
nn_c = nn_c+1;

r_power_coop_notech_nopowertrans = r_power_coop;
Inv_pvwind_mineral(1,1) = inv;% million USD

r_power_coop_notech_nopowertrans2035 = zeros(numpowerunit_all,1);
power_coop_notech_nopowertrans2035 = zeros(numpowerunit_all,1);
r_power_coop_notech_nopowertrans2035(m21to40,:) = r_power_coop_notech_nopowertrans;
power_coop_notech_nopowertrans2035(m21to40,:) = power_coop;
power_mineral2035 = power_mineral;
Inv_pvwind_mineral2035 = Inv_pvwind_mineral;
save('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2035_2040baseline.mat','power_coop_notech_nopowertrans2035','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2035_2040baseline.mat','r_power_coop_notech_nopowertrans2035','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_2035_2040baseline.mat','power_mineral2035','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2035_2040baseline.mat','Inv_pvwind_mineral2035','-v7.3')  %
