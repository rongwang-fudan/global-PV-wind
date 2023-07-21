tic
clear;
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
% load('H:\world\code\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
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

clear index_mineral_off_time2
clear index_mineral_ons_time2
clear index_mineral_pv_time2

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

load('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat');
optpowerunit_IX(:,41) = unitmin;

load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat') % Plant_ID_IX
optpowerunit_IX = optpowerunit_IX(Plant_ID_IX,:);
powerunit_IX = powerunit_IX(Plant_ID_IX,:);
lines_IX_IX = lines_IX_IX(Plant_ID_IX,:);
powerunit_num_IX_IX = powerunit_num_IX_IX(Plant_ID_IX,:);
powerunit_country_IX_IX = powerunit_country_IX_IX(Plant_ID_IX,:);
unitmin = unitmin(Plant_ID_IX,:);
numpowerunit = size(optpowerunit_IX,1);
clear B
clear off_pro_IX
clear off_pro_IX_ori
clear pos
clear powerunit_IX_offshorewind
clear powerunit_IX_onshorewind
clear powerunit_num_IX_onshorewind
clear powerunit_num_IX_onshorewind_ori
clear powerunit_num_IX_PV
clear powerunit_num_IX_PV_ori
[m,n]=find(unitmin~=6);
optpowerunit_IX(m,:) = [];
powerunit_IX(m,:) = [];
lines_IX_IX(m,:) = [];
powerunit_num_IX_IX(m,:) = [];
powerunit_country_IX_IX(m,:) = [];
unitmin(m,:) = [];
numpowerunit = size(optpowerunit_IX,1);
% 
% [B,IX]=sort(optpowerunit_IX(:,20),1);
% optpowerunit_IX = optpowerunit_IX(IX,:);
% powerunit_IX = powerunit_IX(IX,:);
% lines_IX_IX = lines_IX_IX(IX,:);
% powerunit_num_IX_IX = powerunit_num_IX_IX(IX,:);
% powerunit_country_IX_IX = powerunit_country_IX_IX(IX,:);
% unitmin = unitmin(IX,:);
optpowerunit_IX(:,42) = [1:1:size(optpowerunit_IX,1)]';


%%
load('H:\Global PV and wind\Data\mineral_production2021.mat');  %  thousand metric tons/year
mineral_production2021 = mineral_production2021*5;
load('H:\Global PV and wind\Data\mineral_reserve2021.mat');  %  thousand metric tons
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% mineral_production2021的 Molybdenum 年产量（第7列）假设是足够的
% mineral_reserve2021 Molybdenum 和 silicon 储量（第4、7列）假设是足够的
load('H:\Global PV and wind\Data\mineral_CP.mat');
% 9*5
% 9行 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others
% 1:3列为PV, onshorewind和offshore wind的Minerals consumption (kg MW-1)
% 4列：Production (ton/year)
% 5列：Global reserves (ton) 全球储量，其中silicon的储量未知，假设是充足的
load('H:\Global PV and wind\Data\mineral_type.mat');
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% 数值：1是仅PV用，2是仅wind用，3是PV和wind都用得到
mineral_trans_plant1 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant2 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant3 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant4 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant5 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant6 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans_plant8 = zeros(size(optpowerunit_IX,1),192,192);
mineral_trans8 = zeros(192,192,8);
mineral_country = ones(192,8); % 如果此国家此类mineral供不应求，则对应序号为0
Owncountry_ID2 = ones(size(optpowerunit_IX,1),8);
optpowerunit_id = [1:1:size(optpowerunit_IX,1)]';
for iii = 1:8
    optpowerunit_country_final = [];
    for country = 1:1:192
        [ma,n]=find(powerunit_country_IX_IX == country);
        optpowerunit_country = optpowerunit_IX(ma,:);
        optpowerunit_ida = optpowerunit_id(ma,:);
        [B,IX]=sort(optpowerunit_country(:,20),1);
        optpowerunit_country_IX=optpowerunit_country(IX,:); % lat lon
        optpowerunit_ida_IX=optpowerunit_ida(IX,:); 
        clear optpowerunit_ida
        clear optpowerunit_country
        clear mineral
        mineral = zeros(size(ma,1),8);
        for i = iii %:8
            [m,n] = find(optpowerunit_country_IX(:,35)==1); % PV
            mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,1)/1000; % ton
            [m,n] = find(optpowerunit_country_IX(:,35)==2); % onshorewind
            mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,2)/1000; % ton
            [m,n] = find(optpowerunit_country_IX(:,35)==3); % offshorewind
            mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,3)/1000; % ton
        end
        aaa = mineral_production2021(country,:)*1000;
        bbb = mineral_reserve2021(country,:)*1000;
        for i = 1:size(optpowerunit_country_IX,1)
            if min(aaa(:,iii)-mineral(i,iii)) >= 0 & min(bbb(:,iii)-mineral(i,iii)) >= 0
                optpowerunit_country_final = [optpowerunit_country_final; optpowerunit_country_IX(i,:)];
                mineral_reserve2021(country,iii) = mineral_reserve2021(country,iii)-mineral(i,iii)/1000;
                mineral_production2021(country,iii) = mineral_production2021(country,iii)-mineral(i,iii)/1000;
                aaa(:,iii) = aaa(:,iii) - mineral(i,iii);
                bbb(:,iii) = bbb(:,iii) - mineral(i,iii);
                mineral_trans8(country,country,iii) = mineral_trans8(country,country,iii)+mineral(i,iii);
                if iii == 1
                    mineral_trans_plant1(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant1(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 2
                    mineral_trans_plant2(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant2(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 3
                    mineral_trans_plant3(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant3(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 4
                    mineral_trans_plant4(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant4(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 5
                    mineral_trans_plant5(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant5(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 6
                    mineral_trans_plant6(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant6(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
                if iii == 8
                    mineral_trans_plant8(optpowerunit_ida_IX(i),country,country) = mineral_trans_plant8(optpowerunit_ida_IX(i),country,country)+sum(mineral(i,iii));
                end
            else
                mineral_country(country,iii) = 0;
                %                 mineral_reserve2021(country,iii) = 0;
                %                 mineral_production2021(country,iii) = 0;
                %                 aaa(:,iii)= 0;
                %                 bbb(:,iii) = 0;
            end
        end
    end
    Owncountry_ID = unique(optpowerunit_country_final(:,42));
    Owncountry_ID2(Owncountry_ID,iii) = 0; % mineral from own country is 0
    [m,n]=find(mineral_country(:,iii)==0);
    mineral_reserve2021(m,iii) = 0;
    mineral_production2021(m,iii) = 0;
    iii
end
save('H:\Global PV and wind\ANS\Owncountry_ID_type_pro_2040s_8b.mat', 'Owncountry_ID2', '-v7.3')
% mineral from own country is 0
save('H:\Global PV and wind\ANS\mineral_reserve2021_others_2040s_8b.mat', 'mineral_reserve2021', '-v7.3')  %  thousand metric tons
save('H:\Global PV and wind\ANS\mineral_production2021_others_2040s_8b.mat', 'mineral_production2021', '-v7.3')  %  thousand metric tons/year
% 各国自用之后剩余可提供给他国使用的原材料
% 9列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others

%% 各电厂消耗的各类原材料的量
for i = 1:8
    [m,n] = find(optpowerunit_IX(:,35)==1); % PV
    mineral_plant(m,i) = optpowerunit_IX(m,30).* mineral_CP(i,1)/1000; % ton
    [m,n] = find(optpowerunit_IX(:,35)==2); % onshorewind
    mineral_plant(m,i) = optpowerunit_IX(m,30).* mineral_CP(i,2)/1000; % ton
    [m,n] = find(optpowerunit_IX(:,35)==3); % offshorewind
    mineral_plant(m,i) = optpowerunit_IX(m,30).* mineral_CP(i,3)/1000; % ton
end
save('H:\Global PV and wind\ANS\mineral_plant_2040s_8b.mat', 'mineral_plant', '-v7.3')
% mineral consumption, ton
% 9列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others
mineral_trans8_self = mineral_trans8;

%% 计算国家之间的mineral流动方向
aaa = mineral_production2021*1000;
bbb = mineral_reserve2021*1000;
mineral_limit = min(aaa,bbb);
save('H:\Global PV and wind\ANS\mineral_limit_dis_2040s_8b.mat', 'mineral_limit', '-v7.3')
 %  thousand metric tons
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
for iii = 1:1:8
    [m,n]=find(Owncountry_ID2(:,iii)~=0); % 来自于其他国家
    for j = 1:size(m,1)
        cou = powerunit_country_IX_IX(m(j));
        mineral_trans8(:,cou,iii) = mineral_trans8(:,cou,iii)+mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii);
        mineral_limit(:,iii) = mineral_limit(:,iii)-mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii);
        if iii == 1
            mineral_trans_plant1(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant1(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 2
            mineral_trans_plant2(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant2(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 3
            mineral_trans_plant3(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant3(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 4
            mineral_trans_plant4(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant4(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 5
            mineral_trans_plant5(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant5(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 6
            mineral_trans_plant6(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant6(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
        if iii == 8
            mineral_trans_plant8(optpowerunit_id(m(j)),:,cou) = mineral_trans_plant8(optpowerunit_id(m(j)),:,cou)+(mineral_plant(m(j),iii)/sum(mineral_limit(:,iii)).*mineral_limit(:,iii))';
        end
    end
    iii
end
save('H:\Global PV and wind\ANS\mineral_trans8_2040s_8b.mat', 'mineral_trans8', '-v7.3')
% ton，各国之间原材料运输
% 9列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others
save('H:\Global PV and wind\ANS\mineral_trans_plant1_2040s_8b.mat', 'mineral_trans_plant1', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant2_2040s_8b.mat', 'mineral_trans_plant2', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant3_2040s_8b.mat', 'mineral_trans_plant3', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant4_2040s_8b.mat', 'mineral_trans_plant4', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant5_2040s_8b.mat', 'mineral_trans_plant5', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant6_2040s_8b.mat', 'mineral_trans_plant6', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_plant8_2040s_8b.mat', 'mineral_trans_plant8', '-v7.3')
% ton，各国之间原材料运输
% 8 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 8 Rare earths,
sum(sum(sum(mineral_trans8)))-sum(sum(mineral_trans8(:,:,7)))-sum(sum(sum(mineral_trans_plant1)))-sum(sum(sum(mineral_trans_plant2)))-sum(sum(sum(mineral_trans_plant3)))-sum(sum(sum(mineral_trans_plant4)))-sum(sum(sum(mineral_trans_plant5)))-sum(sum(sum(mineral_trans_plant6)))--sum(sum(sum(mineral_trans_plant8)))

%%
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

mineral_trans_reg11_8 = zeros(11,11,8);
for j = 1:1:8
    minerala = mineral_trans8(:,:,j);
    for i = 1:192
        for i2 = 1:192
            mineral_trans_reg11_8(region_ID11(i),region_ID11(i2),j) = mineral_trans_reg11_8(region_ID11(i),region_ID11(i2),j)+minerala(i,i2);
        end
    end
end
mineral_trans_reg5_8 = zeros(8,8,8);
for j = 1:1:8
    minerala = mineral_trans8(:,:,j);
    for i = 1:192
        for i2 = 1:192
            mineral_trans_reg5_8(region_ID(i),region_ID(i2),j) = mineral_trans_reg5_8(region_ID(i),region_ID(i2),j)+minerala(i,i2);
        end
    end
end
save('H:\Global PV and wind\ANS\mineral_trans_reg11_8_2040s_8b.mat', 'mineral_trans_reg11_8', '-v7.3')
save('H:\Global PV and wind\ANS\mineral_trans_reg5_8_2040s_8b.mat', 'mineral_trans_reg5_8', '-v7.3')
% 2040年前间的电厂，2040年碳中和，optimal path

mineral_trans_reg5_8_noself = zeros(8,8,8);
for j = 1:1:8
    minerala = mineral_trans8(:,:,j);
    for i = 1:192
        for i2 = 1:192
            if i ~= i2
            mineral_trans_reg5_8_noself(region_ID(i),region_ID(i2),j) = mineral_trans_reg5_8_noself(region_ID(i),region_ID(i2),j)+minerala(i,i2);
            end
        end
    end
end
save('H:\Global PV and wind\ANS\mineral_trans_reg5_8_2040s_8_noselfb.mat', 'mineral_trans_reg5_8_noself', '-v7.3')

