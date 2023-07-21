%% 2040 C neutrality year, optimal path
tic
clear;
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2060s_test6.mat')  % B_utilize_trans_storage
load('H:\Global PV and wind\ANS\unitmin2040_8_2s.mat'); %% rIp=[595.73, 374.45, 327.22,280.78]./595.73;
load('H:\Global PV and wind\ANS\optpowerunit_IX_8_2040_2s_2060s_test6.mat'); %optpowerunit_IX

load('H:\Global PV and wind\Data\mineral_production2021.mat');  %  thousand metric tons/year
mineral_production_all = sum(mineral_production2021*1000*5)/10^6; % Mt/10years
load('H:\Global PV and wind\Data\mineral_CP.mat');
% 9*5
% 9行 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others
% 1:3列为PV, onshorewind和offshore wind的Minerals consumption (kg MW-1)
% 4列：Production (ton/year)
% 5列：Global reserves (ton) 全球储量，其中silicon的储量未知，假设是充足的
for j = 1:10
    [m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_pv(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,1)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==2 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_ons(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,2)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==3 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_off(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,3)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
end
mineral_8 = (mineral_pv+mineral_ons+mineral_off)/10^6; % Mton
mineral_8(:,9) = [];
r_mineral_8 = mineral_8./repmat(mineral_production_all,10,1);
[m,n]=find(r_mineral_8>1);


% silicon
[m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==4 & B_utilize_trans_storage<=130.85);
mineral_pv5 = optpowerunit_IX(m,30)*mineral_CP(4,1)'/1000/10^6; % MW*(kg MW-1)/1000→Mton, 消耗的量
[mmm,nnn]=find(flip(cumsum(flip(mineral_pv5)))>mineral_production_all(4));
unitmin(m(mmm),1)=3;
clear mineral_pv
clear mineral_ons
clear mineral_off

for j = 1:10
    [m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_pv(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,1)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==2 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_ons(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,2)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==3 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_off(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,3)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
end
mineral_8 = (mineral_pv+mineral_ons+mineral_off)/10^6; % Mton
mineral_8(:,9) = [];
r_mineral_8 = mineral_8./repmat(mineral_production_all,10,1);
[m,n]=find(r_mineral_8>1);

% silicon
[m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==3 & B_utilize_trans_storage<=130.85);
mineral_pv5 = optpowerunit_IX(m,30)*mineral_CP(4,1)'/1000/10^6; % MW*(kg MW-1)/1000→Mton, 消耗的量
[mmm,nnn]=find(flip(cumsum(flip(mineral_pv5)))>mineral_production_all(4));
unitmin(m(mmm),1)=2;
clear mineral_pv
clear mineral_ons
clear mineral_off

for j = 1:10
    [m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_pv(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,1)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==2 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_ons(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,2)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==3 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_off(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,3)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
end
mineral_8 = (mineral_pv+mineral_ons+mineral_off)/10^6; % Mton
mineral_8(:,9) = [];
r_mineral_8 = mineral_8./repmat(mineral_production_all,10,1);
[m,n]=find(r_mineral_8>1);

% silicon
[m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==2 & B_utilize_trans_storage<=130.85);
mineral_pv5 = optpowerunit_IX(m,30)*mineral_CP(4,1)'/1000/10^6; % MW*(kg MW-1)/1000→Mton, 消耗的量
[mmm,nnn]=find(flip(cumsum(flip(mineral_pv5)))>mineral_production_all(4));
unitmin(m(mmm),1)=1;
clear mineral_pv
clear mineral_ons
clear mineral_off

for j = 1:10
    [m,n]=find(optpowerunit_IX(:,35)==1 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_pv(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,1)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==2 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_ons(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,2)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
    [m,n]=find(optpowerunit_IX(:,35)==3 & unitmin==j & B_utilize_trans_storage<=130.85);
    mineral_off(j,:) = sum(optpowerunit_IX(m,30))*mineral_CP(:,3)'/1000; % MW*(kg MW-1)/1000→ton, 消耗的量
end
mineral_8 = (mineral_pv+mineral_ons+mineral_off)/10^6; % Mton
mineral_8(:,9) = [];
r_mineral_8 = mineral_8./repmat(mineral_production_all,10,1);
[m,n]=find(r_mineral_8>1);

save('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat','unitmin');

