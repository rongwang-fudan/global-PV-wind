tic
clear;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_IX_100GW_3_2_all2_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV;

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
% 1.power plant ID;2 country;3 pro;4 county
optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind;

%
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

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
optpowerunit(:,41)=[1:1:size(optpowerunit,1)]';
clear optpowerunit_PV
clear optpowerunit_onshorewind
clear optpowerunit_offshorewind
lines_IX_offshorewind(:,12:18)=0;
lines_IX = [lines_IX_PV;lines_IX_onshorewind;lines_IX_offshorewind];
clear lines_IX_PV
clear lines_IX_onshorewind
clear lines_IX_offshorewind
powerunit_country_IX = [powerunit_num_IX_PV_ori;powerunit_num_IX_onshorewind_ori;off_pro_IX_ori];
clear powerunit_num_IX_PV_ori
clear powerunit_num_IX_onshorewind_ori
clear off_pro_IX_ori
powerunit_num_IX = [powerunit_num_IX_PV(:,5);powerunit_num_IX_onshorewind(:,5);off_pro_IX(:,3)];
clear powerunit_num_IX_PV
clear powerunit_num_IX_onshorewind
clear off_pro_IX
powerunit_pro_IX=lines_IX(:,10);
% [B,IX]=sort(optpowerunit(:,20),1);
% optpowerunit_IX(:,1:40)=optpowerunit(IX,1:40); % lat lon

optpowerunit1 = optpowerunit;
powerunit_country_IX1 = powerunit_country_IX;
powerunit_num_IX1 = powerunit_num_IX;
powerunit_pro_IX1=powerunit_pro_IX;

%%
load('H:\global-PV-wind\ANS\mineral_production_pro.mat');  %  thousand metric tons/year
load('H:\global-PV-wind\ANS\mineral_reserve_pro.mat');  %  thousand metric tons
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% mineral_production2021的 Molybdenum 年产量（第7列）假设是足够的
% mineral_reserve2021 Molybdenum 和 silicon 储量（第4、7列）假设是足够的
mineral_production_pro_all = mineral_production_pro*1000*50; % t/50years
mineral_reserve_pro_all = mineral_reserve_pro*1000; % t/50years
load('H:\global-PV-wind\Data\mineral_production2021.mat');  %  thousand metric tons/year
load('H:\global-PV-wind\Data\mineral_reserve2021.mat');  %  thousand metric tons
mineral_production_cou_all = mineral_production2021*1000*50; % t/50years
mineral_reserve_cou_all = mineral_reserve2021*1000; % t/50years
load('H:\global-PV-wind\Data\mineral_CP.mat');
% 9*5
% 9行 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths, 9 Others
% 1:3列为PV, onshorewind和offshore wind的Minerals consumption (kg MW-1)
% 4列：Production (ton/year)
% 5列：Global reserves (ton) 全球储量，其中silicon的储量未知，假设是充足的
load('H:\global-PV-wind\Data\mineral_type.mat');
% 8列 为各类矿物：1 Copper, 2 Zinc, 3 Nickel, 4 Silicon,
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% 数值：1是仅PV用，2是仅wind用，3是PV和wind都用得到
load('H:\global-PV-wind\Data\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
optpowerunit_country_final = [];
powerunit_pro_IX_final = [];
powerunit_country_IX_final = [];
for country = 1:1:192
    [m1,n]=find(ID_pro(:,3)==country);
    for pro = 1:size(m1,1)
        [mp,np]=find(powerunit_pro_IX==ID_pro(m1(pro),1));
        optpowerunit_pro = optpowerunit(mp,:);
        powerunit_pro_IX2 = powerunit_pro_IX(mp,:);
        powerunit_country_IX2 = powerunit_country_IX(mp,:);
        [B,IX]=sort(optpowerunit_pro(:,20),1);
        optpowerunit_pro_IX=optpowerunit_pro(IX,:); % lat lon
        powerunit_pro_IX_IX = powerunit_pro_IX2(IX,:);
        powerunit_country_IX_IX = powerunit_country_IX2(IX,:);
        clear optpowerunit_pro
        clear powerunit_pro_IX2
        clear powerunit_country_IX2
        clear mineral
        for i = 1:8
            [m,n] = find(optpowerunit_pro_IX(:,35)==1); % PV
            mineral(m,i) = optpowerunit_pro_IX(m,30).* mineral_CP(i,1)/1000; % ton
            [m,n] = find(optpowerunit_pro_IX(:,35)==2); % onshorewind
            mineral(m,i) = optpowerunit_pro_IX(m,30).* mineral_CP(i,2)/1000; % ton
            [m,n] = find(optpowerunit_pro_IX(:,35)==3); % offshorewind
            mineral(m,i) = optpowerunit_pro_IX(m,30).* mineral_CP(i,3)/1000; % ton
        end
        aaa = mineral_production_pro_all(ID_pro(m1(pro),1)+1,:);
        bbb = mineral_reserve_pro_all(ID_pro(m1(pro),1)+1,:);
        for i = 1:size(optpowerunit_pro_IX,1)
            if min(aaa-mineral(i,:)) >= 0 & min(bbb-mineral(i,:)) >= 0
                optpowerunit_country_final = [optpowerunit_country_final; optpowerunit_pro_IX(i,:)];
                powerunit_pro_IX_final = [powerunit_pro_IX_final; powerunit_pro_IX_IX(i,:)];
                powerunit_country_IX_final = [powerunit_country_IX_final; powerunit_country_IX_IX(i,:)];
                aaa = aaa - mineral(i,:);
                bbb = bbb - mineral(i,:);
                mineral_production_cou_all(country,:) = mineral_production_cou_all(country,:)- mineral(i,:); % t/40years
                mineral_reserve_cou_all(country,:) = mineral_reserve_cou_all(country,:)- mineral(i,:);; % t/40years
            end
        end
    end
    [is,pos]=ismember(optpowerunit_country_final(:,41),optpowerunit(:,41));
    pos(pos==0)=[];
    optpowerunit(pos,:)=[];
    powerunit_country_IX(pos,:)=[];
    powerunit_num_IX(pos,:)=[];
    powerunit_pro_IX(pos,:)=[];
    [m,n]=find(powerunit_country_IX == country);
    optpowerunit_country = optpowerunit(m,:);
    powerunit_pro_IX2 = powerunit_pro_IX(m,:);
    powerunit_country_IX2 = powerunit_country_IX(m,:);
    [B,IX]=sort(optpowerunit_country(:,20),1);
    optpowerunit_country_IX=optpowerunit_country(IX,:); % lat lon
    powerunit_pro_IX_IX = powerunit_pro_IX2(IX,:);
    powerunit_country_IX_IX = powerunit_country_IX2(IX,:);
    clear powerunit_pro_IX2
    clear powerunit_country_IX2
    clear optpowerunit_country
    clear mineral
    for i = 1:8
        [m,n] = find(optpowerunit_country_IX(:,35)==1); % PV
        mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,1)/1000; % ton
        [m,n] = find(optpowerunit_country_IX(:,35)==2); % onshorewind
        mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,2)/1000; % ton
        [m,n] = find(optpowerunit_country_IX(:,35)==3); % offshorewind
        mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,3)/1000; % ton
    end
    aaa = mineral_production_cou_all(country,:);
    bbb = mineral_reserve_cou_all(country,:);
    for i = 1:size(optpowerunit_country_IX,1)
        if min(aaa-mineral(i,:)) >= 0 & min(bbb-mineral(i,:)) >= 0
            optpowerunit_country_final = [optpowerunit_country_final; optpowerunit_country_IX(i,:)];
            powerunit_pro_IX_final = [powerunit_pro_IX_final; powerunit_pro_IX_IX(i,:)];
            powerunit_country_IX_final = [powerunit_country_IX_final; powerunit_country_IX_IX(i,:)];
            aaa = aaa - mineral(i,:);
            bbb = bbb - mineral(i,:);
            mineral_production_cou_all(country,:) = mineral_production_cou_all(country,:)- mineral(i,:); % t/40years
            mineral_reserve_cou_all(country,:) = mineral_reserve_cou_all(country,:)- mineral(i,:);; % t/40years
        end
    end
    country
end

%%
[is,pos]=ismember(optpowerunit_country_final(:,41),optpowerunit(:,41));
pos(pos==0)=[];
optpowerunit(pos,:)=[];
powerunit_country_IX(pos,:)=[];
powerunit_pro_IX(pos,:)=[];
powerunit_num_IX(pos,:)=[];
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
for i = 1:1:192
    [m,n]=find(powerunit_country_IX_final==i & optpowerunit_country_final(:,35)==1);
    P_cou_pv(i,1) = sum(optpowerunit_country_final(m,1));
    CP_cou_pv(i,1) = sum(optpowerunit_country_final(m,30))/1000; % GW
    
    [m,n]=find(powerunit_country_IX_final==i & optpowerunit_country_final(:,35)~=1);
    P_cou_wind(i,2) = sum(optpowerunit_country_final(m,1));
    CP_cou_wind(i,2) = sum(optpowerunit_country_final(m,30))/1000; % GW
end

for i = 1:1:size(mineral_reserve_pro,1)
    [m,n]=find(powerunit_pro_IX_final==i-1 & optpowerunit_country_final(:,35)==1);
    if ~isempty(m)
        P_pro_pv(i,1) = sum(optpowerunit_country_final(m,1));
        CP_pro_pv(i,1) = sum(optpowerunit_country_final(m,30))/1000; % GW
    end
    if ~isempty(m)
        [m,n]=find(powerunit_pro_IX_final==i-1 & optpowerunit_country_final(:,35)~=1);
        P_pro_wind(i,2) = sum(optpowerunit_country_final(m,1));
        CP_pro_wind(i,2) = sum(optpowerunit_country_final(m,30))/1000; % GW
    end
end

%%
clear mineral
for i = 1:8
    [m,n] = find(optpowerunit_country_final(:,35)==1); % PV
    mineral(m,i) = optpowerunit_country_final(m,30).* mineral_CP(i,1)/1000; % ton
    [m,n] = find(optpowerunit_country_final(:,35)==2); % onshorewind
    mineral(m,i) = optpowerunit_country_final(m,30).* mineral_CP(i,2)/1000; % ton
    [m,n] = find(optpowerunit_country_final(:,35)==3); % offshorewind
    mineral(m,i) = optpowerunit_country_final(m,30).* mineral_CP(i,3)/1000; % ton
end
mineral_total = sum(mineral);


%%
[B,IX]=sort(optpowerunit(:,20),1);
optpowerunit_country_IX=optpowerunit(IX,:); % lat lon
clear optpowerunit
clear mineral
for i = 1:8
    [m,n] = find(optpowerunit_country_IX(:,35)==1); % PV
    mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,1)/1000; % ton
    [m,n] = find(optpowerunit_country_IX(:,35)==2); % onshorewind
    mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,2)/1000; % ton
    [m,n] = find(optpowerunit_country_IX(:,35)==3); % offshorewind
    mineral(m,i) = optpowerunit_country_IX(m,30).* mineral_CP(i,3)/1000; % ton
end
aaa = sum(mineral_production_cou_all);
bbb = sum(mineral_reserve_cou_all);
for i = 1:size(optpowerunit_country_IX,1)
    if min(aaa-mineral(i,:)) >= 0 & min(bbb-mineral(i,:)) >= 0
        optpowerunit_country_final = [optpowerunit_country_final; optpowerunit_country_IX(i,:)];
        aaa = aaa - mineral(i,:);
        bbb = bbb - mineral(i,:);
    end
end

%%
[B,IX]=sort(optpowerunit1(:,20),1);
optpowerunit_IX(:,:)=optpowerunit1(IX,:); % lat lon
[is,m_mineral]=ismember(optpowerunit_country_final(:,41),optpowerunit_IX(:,41));
save('H:\global-PV-wind\ANS\m_mineral_county0811_2_all_CNfirst_5%_inilow_pro.mat','m_mineral')

[m,n] = find(optpowerunit_country_final(:,35)==1); % PV
index_mineral_pv = optpowerunit_country_final(m,40);
[m,n] = find(optpowerunit_country_final(:,35)==2); % ons
index_mineral_ons = optpowerunit_country_final(m,40);
[m,n] = find(optpowerunit_country_final(:,35)==3); % off
index_mineral_off = optpowerunit_country_final(m,40);
save('H:\global-PV-wind\ANS\index_mineral_pv_county0811_2_all_CNfirst_5%_inilow_pro.mat','index_mineral_pv') % 按照成本排序后保留的PV电厂原始序号
save('H:\global-PV-wind\ANS\index_mineral_ons_county0811_2_all_CNfirst_5%_inilow_pro.mat','index_mineral_ons') % 按照成本排序后保留的onshorewind电厂原始序号
save('H:\global-PV-wind\ANS\index_mineral_off_county0811_2_all_CNfirst_5%_inilow_pro.mat','index_mineral_off') % 按照成本排序后保留的offshorewind电厂原始序号
