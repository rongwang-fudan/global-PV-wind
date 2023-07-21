%%
tic
clear;

%% 筛选出2040年电厂
load('H:\Global PV and wind\ANS\powerunit_country_IX_IX_8_2040_2s_2070_test6xz.mat'); % powerunit_country_IX_IX
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_7_2040_2s_2070_test6.mat')  % B_utilize_trans_storage
B_utilize_trans_storage_nodiffuse = B_utilize_trans_storage;
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2070_test6xz.mat')  % B_utilize_trans_storage
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_192countrya_2040_8_2070.mat')  % B_utilize_trans_storage_192
% num*192,用各国的学习率计算出来的MAC
load('H:\Global PV and wind\ANS\EP_8_2040_2s_2070_test6xz.mat')  % EP，USD/kWh
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_8_2040_2s_2070_test6xz.mat')  % CO2_year_utilize_trans_storage, Mt CO2/year
load('H:\Global PV and wind\ANS\phhall_all2_8_2040_2s_2070_test6xz.mat')
load('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat'); %%
load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat')
load('H:\Global PV and wind\ANS\optpowerunit_IX_8_2040_2s_2070_test6xz.mat'); % optpowerunit_IX
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
clear Plant_ID_IX
EF = CO2_year_utilize_trans_storage./phhall_all2; % 单位发电量的减排量
save('H:\Global PV and wind\ANS\EF_2070_2040.mat','EF','-v7.3')  % t CO2/kWh
save('H:\Global PV and wind\ANS\unitmin_2070_2040.mat','unitmin','-v7.3')  % t CO2/kWh

load('H:\Global PV and wind\ANS\Ph_PVwind5_ef_cou_p_2040_8_2s_6_2a.mat')  % TWh/yr, Ph_PVwind5_ef_cou
load('H:\Global PV and wind\ANS\m_2040_8_2s_6_2a.mat')  % m,由于没有带来更大的收益在2021-2030年被剔除的电厂
Ph_PVwind5_ef_cou_p(m,:,:) = [];
Ph_PVwind5_ef_cou_p2040 = Ph_PVwind5_ef_cou_p(:,:,10);
r = sum(Ph_PVwind5_ef_cou_p2040,2)./phhall_all2;
phhall_all2 = sum(Ph_PVwind5_ef_cou_p(:,:,10),2);
CO2_year_utilize_trans_storage = r.*CO2_year_utilize_trans_storage;
etrans_cou1_num = zeros(size(phhall_all2,1),192,192);
for i = 1:size(phhall_all2,1)
    cc = powerunit_country_IX_IX(i);
    etrans_cou1_num(i,cc,:) = Ph_PVwind5_ef_cou_p2040(i,:);
end
clear Ph_PVwind5_ef_cou_p

numpowerunit_all = size(unitmin,1);
[m,n]=find(unitmin<=10);
m21to40 = m;
[m,n]=find(unitmin>10);
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
Inv_nodiffuse_total = B_utilize_trans_storage_nodiffuse.*(-CO2_year_utilize_trans_storage)+phhall_all2.*EP*10^3; % million USD/yr
Inv_nodiffuse_total(find(isnan(Inv_nodiffuse_total)==1))=0;
Inv = B_utilize_trans_storage.*(-CO2_year_utilize_trans_storage); % million USD/yr
Inv192 = B_utilize_trans_storage_192.*(-CO2_year_utilize_trans_storage); % million USD/yr

% Inv_nodiffuse = B_utilize_trans_storage_nodiffuse.*(-CO2_year_utilize_trans_storage)+phhall_all2.*EP*10^3; % million USD/yr
% Inv = B_utilize_trans_storage.*(-CO2_year_utilize_trans_storage)+phhall_all2.*EP*10^3; % million USD/yr
% Inv192 = B_utilize_trans_storage_192.*(-CO2_year_utilize_trans_storage)+phhall_all2.*EP*10^3; % million USD/yr
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
% 1	North America
% 2	South & Central America
% 3	Europe
% 4	Middle East
% 5	North Africa
% 6	Tropical Africa
% 7	Former Soviet Union
% 8	East Asia
% 9	South Asia
% 10 Southeast Asia
% 11 Pacific Developed region
region_ID11 = region_ID;
region_ID(region_ID==1)=101; % 1 North America
region_ID(region_ID==2)=102; % 2 South & Central America
region_ID(region_ID==3)=103; % 3 Europe
region_ID(region_ID==4)=104; % 4 Middle East
region_ID(region_ID==5)=105; % 5 Africa
region_ID(region_ID==6)=105; % 5 Africa
region_ID(region_ID==7)=106; % 6 Former Soviet Union
region_ID(region_ID==8)=107; % 7 Asia
region_ID(region_ID==9)=107; % 7 Asia
region_ID(region_ID==10)=107; % 7 Asia
region_ID(region_ID==11)=108; % 8 Pacific Developed region
% 1.North America; 2.South & Central America; 3.Europe; 4.Middle East
% 5.Africa; 6.Former Soviet Union; 7.Asia; 8.Pacific Developed region
region_ID = region_ID-100;
powerunit_reg_IX_IX = powerunit_country_IX_IX*0;
for i= 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    powerunit_reg_IX_IX(m,1) = region_ID(i);
end

MAC_CCS = 130.85;
[m,n]=find(B_utilize_trans_storage>MAC_CCS);
B_utilize_trans_storage(B_utilize_trans_storage>MAC_CCS)=MAC_CCS-0.01;
for reg = 1:8
    [m1,n]=find(region_ID==reg);
    [m2,n]=find(powerunit_reg_IX_IX==reg);
    B_min(m2,1) = min(B_utilize_trans_storage_192(m2,m1)')';
end
% [m,n]=find(B_min>MAC_CCS);
for i = 1:size(m,1)
    if B_utilize_trans_storage_192(m(i),powerunit_country_IX_IX(m(i)))==B_min(m(i))
        aaa(i) = 1;
        B_utilize_trans_storage_192(m(i),powerunit_country_IX_IX(m(i))) = MAC_CCS-0.01;
        B_utilize_trans_storage_nodiffuse(m(i),1) = MAC_CCS-0.01;
        asd(i,1) = powerunit_country_IX_IX(m(i));
    else
        [m1,n]=find(region_ID==region_ID(powerunit_country_IX_IX(m(i))));
        [mmm,nnn]=find(B_utilize_trans_storage_192(m(i),m1)==B_min(m(i)));
        B_utilize_trans_storage_192(m(i),m1(nnn(1))) = MAC_CCS-0.01;
        asd(i,1) = nnn(1);
    end
end

for i = 1:8
    [m,n]=find(region_ID==i);
    [m2,n]=find(powerunit_reg_IX_IX==i);
    AAA(i,1) = sum((min(B_utilize_trans_storage_192(m2,m)'))'-B_utilize_trans_storage(m2));
end
save('H:\Global PV and wind\ANS\B_utilize_trans_storage_nodiffuse_2070_2040.mat','B_utilize_trans_storage_nodiffuse','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\B_utilize_trans_storage_2070_2040.mat','B_utilize_trans_storage','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_2070_2040.mat','CO2_year_utilize_trans_storage','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\phhall_all2_2070_2040.mat','phhall_all2','-v7.3')  % TWh/year

%%
load('H:\Global PV and wind\ANS\mineral_limit_dis_2020s_8a.mat')
mineral_limit_20sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2020s_8b.mat')
mineral_limit_20sb = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2030s_8a.mat')
mineral_limit_30sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2030s_8b.mat')
mineral_limit_30sb = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2040s_8a.mat')
mineral_limit_40sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2040s_8b.mat')
mineral_limit_40sb = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2050s_8a.mat')
mineral_limit_50sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2050s_8b.mat')
mineral_limit_50sb = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2060s_8a.mat')
mineral_limit_60sa = mineral_limit;
load('H:\Global PV and wind\ANS\mineral_limit_dis_2060s_8b.mat')
mineral_limit_60sb = mineral_limit;
clear mineral_limit
% mineral_limit, thousand metric tons
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
mineral_limit_20sa(:,7)=0;
mineral_limit_20sb(:,7)=0;
mineral_limit_30sa(:,7)=0;
mineral_limit_30sb(:,7)=0;
mineral_limit_40sa(:,7)=0;
mineral_limit_40sb(:,7)=0;
mineral_limit_50sa(:,7)=0;
mineral_limit_50sb(:,7)=0;
mineral_limit_60sa(:,7)=0;
mineral_limit_60sb(:,7)=0;

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
mineral_limit_40s2a = mineral_limit_40sa./sum(mineral_limit_40sa);
mineral_limit_40s2a(find(isnan(mineral_limit_40s2a)==1))=0;
mineral_limit_40s2b = mineral_limit_40sb./sum(mineral_limit_40sb);
mineral_limit_40s2b(find(isnan(mineral_limit_40s2b)==1))=0;
mineral_limit_50s2a = mineral_limit_50sa./sum(mineral_limit_50sa);
mineral_limit_50s2a(find(isnan(mineral_limit_50s2a)==1))=0;
mineral_limit_50s2b = mineral_limit_50sb./sum(mineral_limit_50sb);
mineral_limit_50s2b(find(isnan(mineral_limit_50s2b)==1))=0;
mineral_limit_60s2a = mineral_limit_60sa./sum(mineral_limit_60sa);
mineral_limit_60s2a(find(isnan(mineral_limit_60s2a)==1))=0;
mineral_limit_60s2b = mineral_limit_60sb./sum(mineral_limit_60sb);
mineral_limit_60s2b(find(isnan(mineral_limit_60s2b)==1))=0;
[m,n]=find(mineral_type~=2); % PV
r_mineral_20sa(:,1) = sum(mineral_limit_20s2a(:,n),2)./sum(sum(mineral_limit_20s2a(:,n)))
r_mineral_20sb(:,1) = sum(mineral_limit_20s2b(:,n),2)./sum(sum(mineral_limit_20s2b(:,n)))
r_mineral_30sa(:,1) = sum(mineral_limit_30s2a(:,n),2)./sum(sum(mineral_limit_30s2a(:,n)))
r_mineral_30sb(:,1) = sum(mineral_limit_30s2b(:,n),2)./sum(sum(mineral_limit_30s2b(:,n)))
r_mineral_40sa(:,1) = sum(mineral_limit_40s2a(:,n),2)./sum(sum(mineral_limit_40s2a(:,n)))
r_mineral_40sb(:,1) = sum(mineral_limit_40s2b(:,n),2)./sum(sum(mineral_limit_40s2b(:,n)))
r_mineral_50sa(:,1) = sum(mineral_limit_50s2a(:,n),2)./sum(sum(mineral_limit_50s2a(:,n)))
r_mineral_50sb(:,1) = sum(mineral_limit_50s2b(:,n),2)./sum(sum(mineral_limit_50s2b(:,n)))
r_mineral_60sa(:,1) = sum(mineral_limit_60s2a(:,n),2)./sum(sum(mineral_limit_60s2a(:,n)))
r_mineral_60sb(:,1) = sum(mineral_limit_60s2b(:,n),2)./sum(sum(mineral_limit_60s2b(:,n)))
[m,n]=find(mineral_type~=1); % Wind
r_mineral_20sa(:,2) = sum(mineral_limit_20s2a(:,n),2)./sum(sum(mineral_limit_20s2a(:,n)))
r_mineral_20sb(:,2) = sum(mineral_limit_20s2b(:,n),2)./sum(sum(mineral_limit_20s2b(:,n)))
r_mineral_30sa(:,2) = sum(mineral_limit_30s2a(:,n),2)./sum(sum(mineral_limit_30s2a(:,n)))
r_mineral_30sb(:,2) = sum(mineral_limit_30s2b(:,n),2)./sum(sum(mineral_limit_30s2b(:,n)))
r_mineral_40sa(:,2) = sum(mineral_limit_40s2a(:,n),2)./sum(sum(mineral_limit_40s2a(:,n)))
r_mineral_40sb(:,2) = sum(mineral_limit_40s2b(:,n),2)./sum(sum(mineral_limit_40s2b(:,n)))
r_mineral_50sa(:,2) = sum(mineral_limit_50s2a(:,n),2)./sum(sum(mineral_limit_50s2a(:,n)))
r_mineral_50sb(:,2) = sum(mineral_limit_50s2b(:,n),2)./sum(sum(mineral_limit_50s2b(:,n)))
r_mineral_60sa(:,2) = sum(mineral_limit_60s2a(:,n),2)./sum(sum(mineral_limit_60s2a(:,n)))
r_mineral_60sb(:,2) = sum(mineral_limit_60s2b(:,n),2)./sum(sum(mineral_limit_60s2b(:,n)))
r_mineral_pv = [r_mineral_20sa(:,1),r_mineral_20sb(:,1),r_mineral_30sa(:,1),r_mineral_30sb(:,1),r_mineral_40sa(:,1),r_mineral_40sb(:,1),r_mineral_50sa(:,1),r_mineral_50sb(:,1),r_mineral_60sa(:,1),r_mineral_60sb(:,1)];
r_mineral_wind = [r_mineral_20sa(:,2),r_mineral_20sb(:,2),r_mineral_30sa(:,2),r_mineral_30sb(:,2),r_mineral_40sa(:,2),r_mineral_40sb(:,2),r_mineral_50sa(:,2),r_mineral_50sb(:,2),r_mineral_60sa(:,2),r_mineral_60sb(:,2)];


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
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2030s_8b.mat')
[m,n]=find(unitmin==4);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2040s_8a.mat')
[m,n]=find(unitmin==5);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2040s_8b.mat')
[m,n]=find(unitmin==6);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2050s_8a.mat')
[m,n]=find(unitmin==7);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2050s_8b.mat')
[m,n]=find(unitmin==8);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2060s_8a.mat')
[m,n]=find(unitmin==9);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
load('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2060s_8b.mat')
[m,n]=find(unitmin==10);
Owncountry_ID2_all(m,:) = Owncountry_ID2;
Owncountry_ID2 = Owncountry_ID2_all;
clear Owncountry_ID2_all

%% 各国自给自足，无需国家之间mineral和power transportion和technology diffussion
power_mineral_cou = zeros(192,40000);
EPpower_mineral_cou = zeros(192,40000);
EFpower_mineral_cou = zeros(192,40000);

plant_use = zeros(size(B_utilize_trans_storage,1),1);
inv = zeros(1,1);
inv_diffuse = zeros(1,1);
power_mineral = zeros(40000,1);
asw = zeros(192,1);
nn_c = 1;

r_power_coop = zeros(numpowerunit,40000); % (i,j)表示两国合作导致的i电厂的有用电量比例
inv_cou = zeros(192,1);
inv_total_cou = zeros(192,1);
for i = 1:192
    j = i;
    %
    [m,n]=find(powerunit_country_IX_IX==i & sum(Owncountry_ID2,2)==0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS );% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,1) = r_power_coop(m,1)+r_etrans_cou1_num(m,i,i);
    power_mineral(1,1) = power_mineral(1,1)+sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    power_mineral_cou(i,2) = power_mineral_cou(i,2)+sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
    inv_cou(i,1) = inv_cou(i,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    inv_total_cou(i,1) = inv_total_cou(i,1)+sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
    inv = inv+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    asw(i,1) = asw(i,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
    plant_use(m,1) = 1;
    %
    %     [m,n]=find(powerunit_country_IX_IX==i & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    for ttt = 1:10
        [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
        r_min = r_mineral_pv(i,ttt);
        r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
        power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
        power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
        inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
        inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        plant_use(m,1) = plant_use(m,1)+r_min;
        
        [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
        r_min = r_mineral_wind(i,ttt);
        r_power_coop(m,1) = r_power_coop(m,1)+r_min*r_etrans_cou1_num(m,i,i);
        power_mineral(1,1) = power_mineral(1,1)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
        power_mineral_cou(i,2) = power_mineral_cou(i,2)+r_min*sum(etrans_cou1_num(m,i,i)); % total power without mineral and power trade
        inv_cou(i,1) = inv_cou(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        inv_total_cou(i,1) = inv_total_cou(i,1)+r_min*sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,i));
        inv = inv+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        asw(i,1) = asw(i,1)+r_min*sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,i));
        plant_use(m,1) = plant_use(m,1)+r_min;
    end
end
save('H:\Global PV and wind\ANS\power_mineral_cou2070_2040.mat','power_mineral_cou','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2070_2040.mat','inv_total_cou','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2070_2040.mat','inv_cou','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\asw_2070_2040.mat','asw','-v7.3')  % million $
nn_c = nn_c+1;
power_mineral_powertrans = power_mineral;
power_mineral_powertrans_tech = power_mineral;
plant_use2 = plant_use;

Inv_pvwind_mineral = zeros(40000,1);
Inv_pvwind_mineral(1,1) = inv;% million USD
Inv_pvwind_mineral_powertrans = Inv_pvwind_mineral;
Inv_pvwind_mineral_powertrans_tech = Inv_pvwind_mineral;

Inv_pvwind_mineral_total = zeros(40000,1);
Inv_pvwind_mineral_total(1,1) = sum(inv_total_cou);% million USD
Inv_pvwind_mineral_powertrans_total = Inv_pvwind_mineral_total;

power_mineral_powertrans_cou = power_mineral_cou;
power_mineral_powertrans_tech_cou = power_mineral_cou;

power_mineral_cou2a = power_mineral_cou;
power_mineral_powertrans_cou2a = power_mineral_cou;
power_mineral_powertrans_tech_cou2a = power_mineral_cou;

r_power_coop_notech = r_power_coop;
r_power_coop_notech_nopowertrans = r_power_coop;
r_power_coop_notech_nopower_nomineraltrans2040 = zeros(numpowerunit_all,1);
r_power_coop_notech_nopower_nomineraltrans2040(m21to40,:) = r_power_coop_notech_nopowertrans(:,1);
save('H:\Global PV and wind\ANS\r_power_coop_notech_nopower_nomineraltrans2070.mat','r_power_coop_notech_nopower_nomineraltrans2040','-v7.3')  % TWh/year

%% 单向建立联系，
% i to j:
% i can transport mineral, power and technology to j
for i = 1:192
    for j = 1:192
        if i~=j
            %             asd = nn_c+1;
            asd = nn_c;
            % only mineral cooperation, country i to country j
            [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            r_min = r_mineral_20sa(i,1);
            power_mineral(nn_c,1) = r_min*sum(etrans_cou1_num(m,j,j)); % total power without mineral and power trade
            plant_use(m,1) = plant_use(m,1)+r_min;
            Inv_pvwind_mineral(asd,1) = sum(r_min.*Inv_nodiffuse(m).*r_etrans_cou1_num(m,j,j));
            Inv_pvwind_mineral_total(asd,1) = sum(r_min.*Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,j,j));
            r_power_coop_notech_nopowertrans(m,nn_c) = r_power_coop_notech_nopowertrans(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
            
            [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            r_min = r_mineral_20sa(i,2);
            power_mineral(nn_c,1) = power_mineral(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j)); % total power without mineral and power trade
            plant_use(m,1) = plant_use(m,1)+r_min;
            Inv_pvwind_mineral(asd,1) = Inv_pvwind_mineral(asd,1)+sum(r_min.*Inv_nodiffuse(m).*r_etrans_cou1_num(m,j,j));
            Inv_pvwind_mineral_total(asd,1) = Inv_pvwind_mineral_total(asd,1)+sum(r_min.*Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,j,j));
            r_power_coop_notech_nopowertrans(m,nn_c) = r_power_coop_notech_nopowertrans(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
            
            for ttt = 2:10
                [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                r_min = r_mineral_pv(i,ttt);
                power_mineral(nn_c,1) = power_mineral(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j)); % total power without mineral and power trade
                plant_use(m,1) = plant_use(m,1)+r_min;
                Inv_pvwind_mineral(asd,1) =  Inv_pvwind_mineral(asd,1)+sum(r_min.*Inv_nodiffuse(m).*r_etrans_cou1_num(m,j,j));
                Inv_pvwind_mineral_total(asd,1) =  Inv_pvwind_mineral_total(asd,1)+sum(r_min.*Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,j,j));
                r_power_coop_notech_nopowertrans(m,nn_c) = r_power_coop_notech_nopowertrans(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
                
                [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                r_min = r_mineral_wind(i,ttt);
                power_mineral(nn_c,1) = power_mineral(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j)); % total power without mineral and power trade
                plant_use(m,1) = plant_use(m,1)+r_min;
                Inv_pvwind_mineral(asd,1) =  Inv_pvwind_mineral(asd,1)+sum(r_min.*Inv_nodiffuse(m).*r_etrans_cou1_num(m,j,j));
                Inv_pvwind_mineral_total(asd,1) =  Inv_pvwind_mineral_total(asd,1)+sum(r_min.*Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,j,j));
                r_power_coop_notech_nopowertrans(m,nn_c) = r_power_coop_notech_nopowertrans(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
            end
            
            
            % mineral and power transportation cooperation, country i to country j
            [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            r_min = r_mineral_20sa(i,1);
            power_mineral_powertrans(nn_c,1) = r_min*sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m,i,j));
            Inv_pvwind_mineral_powertrans(asd,1) = sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse(m1).*r_etrans_cou1_num(m1,j,j));
            Inv_pvwind_mineral_powertrans_total(asd,1) = sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse_total(m1).*r_etrans_cou1_num(m1,j,j));
            r_power_coop_notech(m,nn_c) = r_power_coop_notech(m,nn_c)+r_etrans_cou1_num(m,i,j);
            r_power_coop_notech(m1,nn_c) = r_power_coop_notech(m1,nn_c)+r_min*r_etrans_cou1_num(m1,j,j);
            
            [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
            r_min = r_mineral_20sa(i,2);
            power_mineral_powertrans(nn_c,1) = power_mineral_powertrans(nn_c,1)+r_min*sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m,i,j));
            Inv_pvwind_mineral_powertrans(asd,1) = Inv_pvwind_mineral_powertrans(asd,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse(m1).*r_etrans_cou1_num(m1,j,j));
            Inv_pvwind_mineral_powertrans_total(asd,1) = Inv_pvwind_mineral_powertrans_total(asd,1)+sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse_total(m1).*r_etrans_cou1_num(m1,j,j));
            r_power_coop_notech(m,nn_c) = r_power_coop_notech(m,nn_c)+r_etrans_cou1_num(m,i,j);
            r_power_coop_notech(m1,nn_c) = r_power_coop_notech(m1,nn_c)+r_min*r_etrans_cou1_num(m1,j,j);
            
            for ttt = 2:10
                [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                r_min = r_mineral_pv(i,ttt);
                power_mineral_powertrans(nn_c,1) = power_mineral_powertrans(nn_c,1)+r_min*sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m,i,j));
                Inv_pvwind_mineral_powertrans(asd,1) = Inv_pvwind_mineral_powertrans(asd,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse(m1).*r_etrans_cou1_num(m1,j,j));
                Inv_pvwind_mineral_powertrans_total(asd,1) = Inv_pvwind_mineral_powertrans_total(asd,1)+sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse_total(m1).*r_etrans_cou1_num(m1,j,j));
                r_power_coop_notech(m,nn_c) = r_power_coop_notech(m,nn_c)+r_etrans_cou1_num(m,i,j);
                r_power_coop_notech(m1,nn_c) = r_power_coop_notech(m1,nn_c)+r_min*r_etrans_cou1_num(m1,j,j);
                
                [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & sum(Owncountry_ID2,2)~=0 & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                [m,n]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & B_utilize_trans_storage_nodiffuse <=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
                r_min = r_mineral_wind(i,ttt);
                power_mineral_powertrans(nn_c,1) = power_mineral_powertrans(nn_c,1)+r_min*sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m,i,j));
                Inv_pvwind_mineral_powertrans(asd,1) = Inv_pvwind_mineral_powertrans(asd,1)+sum(Inv_nodiffuse(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse(m1).*r_etrans_cou1_num(m1,j,j));
                Inv_pvwind_mineral_powertrans_total(asd,1) = Inv_pvwind_mineral_powertrans_total(asd,1)+sum(Inv_nodiffuse_total(m).*r_etrans_cou1_num(m,i,j))+sum(r_min.*Inv_nodiffuse_total(m1).*r_etrans_cou1_num(m1,j,j));
                r_power_coop_notech(m,nn_c) = r_power_coop_notech(m,nn_c)+r_etrans_cou1_num(m,i,j);
                r_power_coop_notech(m1,nn_c) = r_power_coop_notech(m1,nn_c)+r_min*r_etrans_cou1_num(m1,j,j);
            end
            
            index_ij(nn_c,1) = i;
            index_ij(nn_c,2) = j;
            nn_c = nn_c+1;
        end
    end
    i
end
save('H:\Global PV and wind\ANS\power_mineraltrans_cou2070_2040.mat','power_mineral_cou','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_cou2070_2040.mat','power_mineral_powertrans_cou','-v7.3')  % TWh/year
r_power_coop_notech_nopowertrans(:,nn_c:end)=[];
Inv_pvwind_mineral(nn_c:end)=[];
power_mineral(nn_c:end)=[];
r_power_coop_notech_nopowertrans2040 = zeros(numpowerunit_all,nn_c-1);
r_power_coop_notech_nopowertrans2040(m21to40,:) = r_power_coop_notech_nopowertrans;
power_mineral2040 = power_mineral;
Inv_pvwind_mineral2040 = Inv_pvwind_mineral;
Inv_pvwind_mineral_total(nn_c:end)=[];
Inv_pvwind_mineral2040_total = Inv_pvwind_mineral_total;
save('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2070.mat','r_power_coop_notech_nopowertrans2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_2070.mat','power_mineral2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2070.mat','Inv_pvwind_mineral2040','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2070_total.mat','Inv_pvwind_mineral2040_total','-v7.3')  %

r_power_coop_notech(:,nn_c:end)=[];
Inv_pvwind_mineral_powertrans(nn_c:end)=[];
power_mineral_powertrans(nn_c:end)=[];
r_power_coop_notech2040 = zeros(numpowerunit_all,nn_c-1);
r_power_coop_notech2040(m21to40,:) = r_power_coop_notech;
power_mineral_powertrans2040 = power_mineral_powertrans;
Inv_pvwind_mineral_powertrans2040 = Inv_pvwind_mineral_powertrans;
Inv_pvwind_mineral_powertrans_total(nn_c:end)=[];
Inv_pvwind_mineral_powertrans2040_total = Inv_pvwind_mineral_powertrans_total;
save('H:\Global PV and wind\ANS\r_power_coop_notech_2070.mat','r_power_coop_notech2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_2070.mat','power_mineral_powertrans2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_2070.mat','Inv_pvwind_mineral_powertrans2040','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_2070_total.mat','Inv_pvwind_mineral_powertrans2040_total','-v7.3')  %

[B,IX2] = sort(power_mineral_powertrans(2:nn_c-1),'descend');
index_ij_IX = index_ij(IX2+1,1:2);
save('H:\Global PV and wind\ANS\index_ij_2070.mat','index_ij')  % TWh/year
save('H:\Global PV and wind\ANS\index_ij_IX_2070.mat','index_ij_IX')  %

%% With mineral trade, power trans and technology diffussion
ind_cou = zeros(5000,192);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);% find power plant situated in ocuntry i and without mineral trade
    ind_cou(1:size(m,1),i) = m;
end

nn_c = 2;
Inv192min = Inv_nodiffuse;
Bmin = B_utilize_trans_storage_nodiffuse;

load('H:\Global PV and wind\ANS\Inn_2070.mat')
Inn_2055_2040_ixx = flip(Inn_2070(:,1:2));
Inn_2055_2040_num = Inn_2070(:,3:size(Inn_2070,2));
for i = 1:size(Inn_2055_2040_ixx,1)
    [m,n]=find(index_ij_IX(:,1)==Inn_2055_2040_ixx(i,1) & index_ij_IX(:,2)==Inn_2055_2040_ixx(i,2));
    index_ij_IX(m,:) = [];
    index_ij_IX = [Inn_2055_2040_ixx(i,1:2);index_ij_IX];
end
save('H:\Global PV and wind\ANS\index_ij_IX_2070xz.mat','index_ij_IX')  %

for i2 = 1:1:size(Inn_2055_2040_ixx,1)
    i = index_ij_IX(i2,1);
    j = index_ij_IX(i2,2);
    asd = nn_c;
    
%     if Inn_2055_2040_num(i2,4+1)==0
%         asd1 = 1;
%     else if Inn_2055_2040_num(i2,4*1+1)~=0
%             asd1 = 2;
%         end
%     end    
    asd1 = size(Inn_2055_2040_num,2)/3;
    Invtest = ones(size(Inv,1),size(Inv,2))*10^10;
    for nnnum = 1:asd1
        [m,n]=find(optpowerunit_IX_IX(:,35)==Inn_2055_2040_num(i2,(nnnum-1)*3+1)  & powerunit_reg_IX_IX==Inn_2055_2040_num(i2,(nnnum-1)*3+2) & unitmin==Inn_2055_2040_num(i2,(nnnum-1)*3+3));
        Inv192(m,j) = Inv(m);
        Invtest(m,1) =  Inv(m);
        B_utilize_trans_storage_192(m,j) = B_utilize_trans_storage(m);
    end
    Inv192min1 = Inv192min;
    % mineral, power transportation and technology diffusion cooperation, country i to country j
    mm = ind_cou(:,j);
    mm(mm==0)=[];
    Inv2(mm,1) = min(Inv192(mm,j),Inv192(mm,i));
    B2(mm,1) = min(B_utilize_trans_storage_192(mm,j),B_utilize_trans_storage_192(mm,i));
    asw1 = sum(Inv192min(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,j,j));
    [mvvv,nvvv]=find(powerunit_country_IX_IX==j);
    %     Inv192min(mvvv) = min(Inv192(mvvv,i),Inv192min(mvvv));
    Inv192min(mvvv) = min(Inv2(mvvv,1),Inv192min(mvvv));
    Bmin(mvvv) = min(B2(mvvv,1),Bmin(mvvv));
    asw(j,1) = sum(Inv192min(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,j,j));
    xxx = asw(j,1)-asw1;
    xxx(xxx>0)=0;
    xxxall(i2) = xxx;
    
    [m1,n1]=find(index_ij_IX(1:i2-1,1)==j);
    if size(m1,1)>0
    for iyu = 1:size(m1,1)
    mm = ind_cou(:,index_ij_IX(m1(iyu),2));
    mm(mm==0)=[];
    Inv2(mm,1) = min(Inv192(mm,index_ij_IX(m1(iyu),2)),Invtest(mm,1));
    asw1 = sum(Inv192min1(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,index_ij_IX(m1(iyu),2),index_ij_IX(m1(iyu),2)));
    [mvvv,nvvv]=find(powerunit_country_IX_IX==index_ij_IX(m1(iyu),2));
    Inv192min(mvvv) = min(Inv2(mvvv,1),Inv192min(mvvv));
    asw(index_ij_IX(m1(iyu),2),1) = sum(Inv192min(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,index_ij_IX(m1(iyu),2),index_ij_IX(m1(iyu),2)));
    aqa = asw(index_ij_IX(m1(iyu),2),1)-asw1;
    aqa(aqa>0)=0;
    xxx = xxx+aqa;
    end
    end
    
    
    
    %     [m,n]=find(powerunit_country_IX_IX==j & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0 & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,1);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse>MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = xxx+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,2);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse>MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    %     [m,n]=find(powerunit_country_IX_IX==j & unitmin==2 & plant_use2<1 & sum(Owncountry_ID2,2)~=0 & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    for ttt = 2:1:10
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_pv(i,ttt);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse >MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_wind(i,ttt);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse >MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    end
    
    nn_c = nn_c+1;
    i2
end

for i2 = size(Inn_2055_2040_ixx,1)+1:1:size(index_ij_IX,1)
    i = index_ij_IX(i2,1);
    j = index_ij_IX(i2,2);
    %     asd = nn_c+1;
    asd = nn_c;
    xxx = 0;
    if region_ID(i) == region_ID(j)
        % mineral, power transportation and technology diffusion cooperation, country i to country j
        mm = ind_cou(:,j);
        mm(mm==0)=[];
        asw1 = sum(Inv192min(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,j,j));
        [mvvv,nvvv]=find(powerunit_country_IX_IX==j);
        Inv192min(mvvv) = min(Inv192(mvvv,i),Inv192min(mvvv));
        Bmin(mvvv) = min(B_utilize_trans_storage_192(mvvv,i),Bmin(mvvv));
        asw(j,1) = sum(Inv192min(mm,1).*plant_use2(mm,1).*r_etrans_cou1_num(mm,j,j));
        xxx = asw(j,1)-asw1;
        xxx(xxx>0)=0;
        xxxall(i2,1) = xxx;
    end
    
    %     [m,n]=find(powerunit_country_IX_IX==j & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0 & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,1);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse>MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==1 & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = xxx+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_20sa(i,2);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse>MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==1 & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    %     [m,n]=find(powerunit_country_IX_IX==j & unitmin==2 & plant_use2<1 & sum(Owncountry_ID2,2)~=0 & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    for ttt = 2:1:10
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_pv(i,ttt);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse >MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)==1 & unitmin==ttt & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    
    [m,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)~=0);% find power plant situated in ocuntry i and without mineral trade
    r_min = r_mineral_wind(i,ttt);
    plant_use2(m,1) = plant_use2(m,1)+r_min;
    [m1,n]=find(powerunit_country_IX_IX==j & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & plant_use2<1 & sum(Owncountry_ID2,2)==0 &  B_utilize_trans_storage_nodiffuse >MAC_CCS & min(B_utilize_trans_storage_192(:,i),B_utilize_trans_storage_192(:,j))<=MAC_CCS & region_ID(i) == region_ID(j));% find power plant situated in ocuntry i and without mineral trade
    plant_use2(m1,1) = 1;
    [m2,n2]=find(powerunit_country_IX_IX==i & optpowerunit_IX_IX(:,35)~=1 & unitmin==ttt & B_utilize_trans_storage<=MAC_CCS);% find power plant situated in ocuntry i and without mineral trade
    r_power_coop(m,nn_c) = r_power_coop(m,nn_c)+r_min*r_etrans_cou1_num(m,j,j);
    r_power_coop(m1,nn_c) = r_power_coop(m1,nn_c)+r_etrans_cou1_num(m1,j,j);
    r_power_coop(m2,nn_c) = r_power_coop(m2,nn_c)+r_etrans_cou1_num(m2,i,j);
    power_mineral_powertrans_tech(nn_c,1) = power_mineral_powertrans_tech(nn_c,1)+r_min*sum(etrans_cou1_num(m,j,j))+sum(etrans_cou1_num(m1,j,j))+sum(etrans_cou1_num(m2,i,j)); % total power without mineral and power trade
    Inv_pvwind_mineral_powertrans_tech(asd,1) = Inv_pvwind_mineral_powertrans_tech(asd,1)+sum(Inv(m2).*r_etrans_cou1_num(m2,i,j))+sum(Inv192min(m1,1).*r_etrans_cou1_num(m1,j,j))+sum(r_min.*Inv192min(m,1).*r_etrans_cou1_num(m,j,j));
    end
    
    nn_c = nn_c+1;
    i2
end
r_power_coop(:,nn_c:end)=[];
Inv_pvwind_mineral_powertrans_tech(nn_c:end)=[];
power_mineral_powertrans_tech(nn_c:end)=[];
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_cou2070_2040.mat','power_mineral_powertrans_tech_cou','-v7.3')  % TWh/year

% a1 = sum(r_power_coop,2);
% [m,n]=find(phhall_all2==0);
% a1(m)=1;
% sum(phhall_all2)
% sum(a1.*phhall_all2)

r_power_coop2040 = zeros(numpowerunit_all,nn_c-1);
r_power_coop2040(m21to40,:) = r_power_coop;
power_mineral_powertrans_tech2040 = power_mineral_powertrans_tech;
Inv_pvwind_mineral_powertrans_tech2040 = Inv_pvwind_mineral_powertrans_tech;

save('H:\Global PV and wind\ANS\r_power_coop_2070.mat','r_power_coop2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070.mat','power_mineral_powertrans_tech2040','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070.mat','Inv_pvwind_mineral_powertrans_tech2040','-v7.3')  %

%%
% sum(Inv_pvwind_mineral_powertrans_tech2030)-sum(xxxall)
sum(Inv)
sum(Inv_nodiffuse)
sum(Inv192min)
sum(Inv_pvwind_mineral_powertrans_tech2040)

