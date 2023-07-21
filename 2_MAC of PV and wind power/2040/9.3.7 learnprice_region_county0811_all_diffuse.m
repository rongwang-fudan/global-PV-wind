tic
clear;
load('H:\Global PV and wind\Data\initialcost_ratio_country_0111low_off.mat')  % 1.PV 2.Onshorewind的成本与China的对比
initialcost_ratio_country_off = initialcost_ratio_country;
load('H:\Global PV and wind\Data\initialcost_ratio_country_0111low.mat')  % 1.PV 2.Onshorewind的成本与China的对比
initialcost_ratio_country2 = [initialcost_ratio_country initialcost_ratio_country_off];
clear initialcost_ratio_country
clear initialcost_ratio_country_off
initialcost_ratio_country =initialcost_ratio_country2;
load('H:\Global PV and wind\Data\Country_classify.mat')
% 1:developed country; 2: developing country; 3: Least Developed Countries
Country_classify(Country_classify>=2)=2;
% load('H:\Global PV and wind\ANS\CP2020_country.mat'); % GW 2020
% 1 Nuclear, 2	Hydroelectricity,3 Tide and wave, 4 Hydroelectric pumped storage, 5 Geothermal,
% 6 Solar, 7 Wind, 8 Biomass and waste
load('H:\Global PV and wind\Data\CP_PV_ons_off2020.mat'); % MW 2020
% 1 Solar photovoltaic, 2 Onshore wind energy, 3 Offshore wind energy
CP_max = CP_PV_ons_off2020(35,:);  % MW 2020

% 实际各国learning rate，当为nan和inf时和全球统一
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
load('H:\Global PV and wind\ANS\unitmin2040_8_2s.mat'); %
optpowerunit_IX(:,41) = unitmin;
[M,N]=find(optpowerunit_IX(:,35)==1);
powerunit_num_IX_PV_ori = powerunit_country_IX_IX(M);
optpowerunit_PV = optpowerunit_IX(M,:);
[M,N]=find(optpowerunit_IX(:,35)==2);
powerunit_num_IX_onshorewind_ori = powerunit_country_IX_IX(M);
optpowerunit_onshorewind = optpowerunit_IX(M,:);
[M,N]=find(optpowerunit_IX(:,35)==3);
off_pro_IX_ori = powerunit_country_IX_IX(M);
optpowerunit_offshorewind = optpowerunit_IX(M,:);

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

learnprice_region = ones(192,6);
learnprice_region(:,1:2) = learnprice_region(:,1:2).*initialcost_ratio_country(:,1);
learnprice_region(:,3:4) = learnprice_region(:,3:4).*initialcost_ratio_country(:,2);
learnprice_region(:,5:6) = learnprice_region(:,5:6).*initialcost_ratio_country(:,3);
LR = ones(11,1);
for time = 1:1:10
    CP_PV_ons_off2020_0 = CP_PV_ons_off2020;
    for country = 1:1:192
        [id_m_pv,n]=find(powerunit_num_IX_PV_ori(:,1)==country & optpowerunit_PV(:,41)==time);
        [id_m_ons,n]=find(powerunit_num_IX_onshorewind_ori(:,1)==country & optpowerunit_onshorewind(:,41)==time);
        [id_m_off,n]=find(off_pro_IX_ori(:,1)==country & optpowerunit_offshorewind(:,41)==time);
        num_plant(country,1) = size(id_m_pv,1); % PV power plant number in each country
        num_plant(country,2) = size(id_m_ons,1); % onshore wind power plant number in each country
        num_plant(country,3) = size(id_m_off,1); % offshore wind power plant number in each country
        
        optpowerunit = [optpowerunit_PV(id_m_pv,:);optpowerunit_onshorewind(id_m_ons,:);optpowerunit_offshorewind(id_m_off,:)];
        if ~isempty(optpowerunit)
            [B,IX]=sort(optpowerunit(:,20),1);
            numpowerunit = size(optpowerunit,1);
            for i=1:numpowerunit
                i2=IX(i);
                powerunit_IX(i,1)=i2;
                optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
            end
            
            power_country(country,1) = sum(optpowerunit_IX(:,1));
            
            % PV
            if Country_classify(country)==1 % developed country 发达国家
                LR_PV1 = 0.18; % PV module
                LR_PV2 = 0.1643; % BOS
            else
                LR_PV1 = 0.26; %0.18; % 0.37; %  0.37 PV module
                LR_PV2 = 0.25; % BOS
            end
            CP0_PV = CP_PV_ons_off2020(country,1); %  cumulative capacity potential in China in 2020, MW
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
            if CP0_PV<CP_max(1)
                CP0_PV = CP_max(1);
            end
            r_pv = CP_PV./CP0_PV;
            r_pv(r_pv<1)=1;
            CP_PV_ons_off2020(country,1) = CP_PV(end);
            
            for i = 1:size(r_pv,1)
                LR_PV(i,1) = LR_PV1;
                learnprice_PV(i,1) = (r_pv(i,1)).^(log(1- LR_PV(i,1))/log(2));
                LR_PV_2(i,1) = LR_PV2;
                learnprice_PV_2(i,1) = (r_pv(i,1)).^(log(1- LR_PV_2(i,1))/log(2));
            end
            clear LR_PV
            clear LR_PV_2
            learnprice_PV(find(isnan(learnprice_PV)==1))=1;
            learnprice_PV_2(find(isnan(learnprice_PV_2)==1))=1;
            learnprice_region(country,1)=learnprice_PV(end)*learnprice_region(country,1);
            learnprice_region(country,2)=learnprice_PV_2(end)*learnprice_region(country,2);
            for i = 1:size(CP_PV,1)-1
                CP_PV22(i,1) = CP_PV(i+1,1)-CP_PV(1,1);
                CP_PV22(i,2) = (CP_PV(i+1,1)-CP_PV(1,1))/(CP_PV(end,1)-CP_PV(1,1));
            end
            cost_PV = zeros(numpowerunit,19);
            if optpowerunit_IX(1,35)==1
                cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
            end
            for i = 2: numpowerunit
                if optpowerunit_IX(i,35)==1
                    cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
                else
                    cost_PV(i,11:19) =  cost_PV(i-1,11:19);
                end
            end
            clear CP_PV
            clear CP0_PV
            
            
            % onshorewind
            if Country_classify(country)==1 % developed country 发达国家
                LR_onshorewind1 = 0.073; %0.18; % 0.37; %  0.37 PV module
                LR_onshorewind2= 0.1643;
            else
                LR_onshorewind1 = 0.122; %0.18; % 0.37; %  0.37 PV module
                LR_onshorewind2= 0.25;
            end
            CP0_onshorewind = CP_PV_ons_off2020(country,2); % 278324; % MW
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
            if CP0_onshorewind<CP_max(2)
                CP0_onshorewind = CP_max(2);
            end
            r_onshorewind = CP_onshorewind./CP0_onshorewind;
            r_onshorewind(r_onshorewind<1)=1;
            CP_PV_ons_off2020(country,2) = CP_onshorewind(end);
            
            for i = 1:size(r_onshorewind,1)
                LR_onshorewind(i,1) = LR_onshorewind1;
                learnprice_onshorewind(i,1) = (r_onshorewind(i,1)).^(log(1- LR_onshorewind(i,1))/log(2));
                LR_onshorewind_2(i,1) = LR_onshorewind2;
                learnprice_onshorewind_2(i,1) = (r_onshorewind(i,1)).^(log(1- LR_onshorewind_2(i,1))/log(2));
            end
            clear LR_onshorewind
            clear LR_onshorewind_2
            learnprice_onshorewind(find(isnan(learnprice_onshorewind)==1))=1;
            learnprice_onshorewind_2(find(isnan(learnprice_onshorewind_2)==1))=1;
            learnprice_region(country,3)=learnprice_onshorewind(end)*learnprice_region(country,3);
            learnprice_region(country,4)=learnprice_onshorewind_2(end)*learnprice_region(country,4);
            for i = 1:size(CP_onshorewind,1)-1
                CP_onshorewind22(i,1) = CP_onshorewind(i+1,1)-CP_onshorewind(1,1);
                CP_onshorewind22(i,2) = (CP_onshorewind(i+1,1)-CP_onshorewind(1,1))/(CP_onshorewind(end,1)-CP_onshorewind(1,1));
            end
            cost_onshorewind = zeros(numpowerunit,19);
            if optpowerunit_IX(1,35)==2
                cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
            end
            for i = 2: numpowerunit
                if optpowerunit_IX(i,35)==2
                    cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
                else
                    cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
                end
            end
            clear CP_onshorewind
            clear CP0_onshorewind
            
            % offhorewind
            if Country_classify(country)==1 % developed country 发达国家
                LR_offshorewind1 = 0.073; %0.18; % 0.37; %  0.37 PV module
                LR_offshorewind2= 0.1643;
            else
                LR_offshorewind1 = 0.122; %0.18; % 0.37; %  0.37 PV module
                LR_offshorewind2= 0.25;
            end
            CP0_offshorewind = CP_PV_ons_off2020(country,3); % 8990;% 11130; % MW
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
            if CP0_offshorewind<CP_max(3)
                CP0_offshorewind = CP_max(3);
            end
            r_offshorewind = CP_offshorewind./CP0_offshorewind;
            r_offshorewind(r_offshorewind<1)=1;
            CP_PV_ons_off2020(country,3) = CP_offshorewind(end);
            
            for i = 1:size(r_offshorewind,1)
                LR_offshorewind(i,1) = LR_offshorewind1;
                learnprice_offshorewind(i,1) = (r_offshorewind(i,1)).^(log(1- LR_offshorewind(i,1))/log(2));
                LR_offshorewind_2(i,1) = LR_offshorewind2;
                learnprice_offshorewind_2(i,1) = (r_offshorewind(i,1)).^(log(1- LR_offshorewind_2(i,1))/log(2));
            end
            clear LR_offshorewind
            clear LR_offshorewind_2
            learnprice_offshorewind(find(isnan(learnprice_offshorewind)==1))=1;
            learnprice_offshorewind_2(find(isnan(learnprice_offshorewind_2)==1))=1;
            learnprice_region(country,5)=learnprice_offshorewind(end)*learnprice_region(country,5);
            learnprice_region(country,6)=learnprice_offshorewind_2(end)*learnprice_region(country,6);
            
            cost_offshorewind = zeros(numpowerunit,10);
            cost_storage_offshorewind = zeros(numpowerunit,1);
            if optpowerunit_IX(1,35)==3
                cost_offshorewind(1,6:10) =  optpowerunit_IX(1,6:10);
            end
            for i = 2: numpowerunit
                if optpowerunit_IX(i,35)==3
                    cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10) + optpowerunit_IX(i,6:10);
                else
                    cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10);
                end
            end
        end
        clear CP_offshorewind
        clear CP0_offshorewind
        
        clear powerunit_IX;
        clear optpowerunit_IX;
        clear learnprice_PV;
        clear learnprice_PV_2;
        clear CP_PV22;
        clear learnprice_onshorewind;
        clear learnprice_onshorewind_2;
        clear CP_onshorewind22;
        clear CP_offshorewind22;
        clear learnprice_offshorewind;
        clear learnprice_offshorewind_2;
        clear P_ratio_3type
        clear unit_PV;
        clear unit_onshorewind;
        clear unit_offshorewind;
        clear penex;
        clear pene_real;
        clear pene;
        clear penex_real;
        clear penemin_real;
        clear P_ratio_3type;clear pe_ratio
        
        clear pene_reg;
        clear costmax1;
        clear costmin1
        clear penemax1
        clear penemin1
        country
    end
    
    %%
    
    %     for i = 1:6
    %         [m,n]=find(learnprice_region(:,i)>min(learnprice_region(:,i)));
    %         learnprice_region(m,i) = min(learnprice_region(:,i));
    %     end
    learnprice_region_d192 = learnprice_region;
    for i = 1:6
        for reg = 1:8
            [mmm,nnn]=find(region_ID==reg);
            [m2,n2]=find(learnprice_region(:,i)==min(learnprice_region(mmm,i)) & region_ID==reg);
            if size(m2,1)==1
                index_cou(reg,i,time) = m2;
            else
                index_cou(reg,i,time) = 0;
            end
            [m,n]=find(learnprice_region(mmm,i)>min(learnprice_region(mmm,i)));
            learnprice_region(mmm,i) = min(learnprice_region(mmm,i));
            learnprice_region2(mmm,i,time+1) = min(learnprice_region(mmm,i));
        end
    end
    learnprice_region_d192(:,1) = learnprice_region_d192(:,1)./initialcost_ratio_country(:,1);
    learnprice_region_d192(:,2) = learnprice_region_d192(:,2)./initialcost_ratio_country(:,1);
    learnprice_region_d192(:,3) = learnprice_region_d192(:,3)./initialcost_ratio_country(:,2);
    learnprice_region_d192(:,4) = learnprice_region_d192(:,4)./initialcost_ratio_country(:,2);
    learnprice_region_d192(:,5) = learnprice_region_d192(:,5)./initialcost_ratio_country(:,3);
    learnprice_region_d192(:,6) = learnprice_region_d192(:,6)./initialcost_ratio_country(:,3);
    filename = strcat('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_pro',num2str(time+1),'2_192.mat');
    save(filename, 'learnprice_region_d192','-v7.3');
    time
end
clear learnprice_region
learnprice_region1 = learnprice_region2(:,:,1:11);
% save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro_8_2040.mat','learnprice_region'); % 10 regions * 6
% 第1，3列表示的是成本
% 其余是学习率

% learnprice_region = ones(192,6);
learnprice_region(:,1) = learnprice_region1(:,1,2)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,2)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,2)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,2)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,2)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,2)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro22_8_2040.mat','learnprice_region'); % 10 regions * 6

learnprice_region(:,1) = learnprice_region1(:,1,3)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,3)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,3)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,3)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,3)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,3)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro32_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,4)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,4)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,4)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,4)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,4)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,4)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro42_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,5)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,5)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,5)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,5)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,5)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,5)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro52_8_2040.mat','learnprice_region'); %


learnprice_region(:,1) = learnprice_region1(:,1,6)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,6)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,6)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,6)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,6)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,6)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro62_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,7)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,7)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,7)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,7)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,7)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,7)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro72_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,8)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,8)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,8)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,8)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,8)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,8)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro82_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,9)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,9)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,9)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,9)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,9)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,9)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro92_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,10)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,10)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,10)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,10)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,10)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,10)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro102_8_2040.mat','learnprice_region'); %

learnprice_region(:,1) = learnprice_region1(:,1,11)./initialcost_ratio_country(:,1);
learnprice_region(:,2) = learnprice_region1(:,2,11)./initialcost_ratio_country(:,1);
learnprice_region(:,3) = learnprice_region1(:,3,11)./initialcost_ratio_country(:,2);
learnprice_region(:,4) = learnprice_region1(:,4,11)./initialcost_ratio_country(:,2);
learnprice_region(:,5) = learnprice_region1(:,5,11)./initialcost_ratio_country(:,3);
learnprice_region(:,6) = learnprice_region1(:,6,11)./initialcost_ratio_country(:,3);
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_2_diffuse_pro112_8_2040.mat','learnprice_region'); %


