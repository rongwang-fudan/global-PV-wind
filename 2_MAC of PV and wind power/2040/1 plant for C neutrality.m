% 实际各国learning rate，当为nan和inf时和全球统一
tic
clear;
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

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
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

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_8_2070.mat')
etrans_t = sum(sum(etrans_cou1_num,3),2);
load('H:\global-PV-wind\ANS\powerunit_country_IX_IX_8_2070_2s_2060s_test6.mat'); % powerunit_country_IX_IX
for i = 1:size(etrans_t,1)
    etrans_self(i,1) = etrans_cou1_num(i,powerunit_country_IX_IX(i),powerunit_country_IX_IX(i));
    a = reshape(etrans_cou1_num(i,:,:),[192,192]);
    [m,n]=find(a~=0);
    nn(i,1) = size(m,1);
end
[Plant_ID_others,n]=find(etrans_self./etrans_t<0.05); 
% 认为这些电厂是为了给别的国家输送电力而建设
etrans_t(Plant_ID_others,1) = 0;
ID_others = zeros(size(etrans_t,1),1);
ID_others(Plant_ID_others,1) = 1; % 为别的国家而建的电厂的值为1
load('H:\global-PV-wind\ANS\unitmin_global_IX_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat');
% load('H:\global-PV-wind\ANS\unitmin_global_IX_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070xz.mat');

optpowerunit_IX_self = optpowerunit_IX;
optpowerunit_IX_self(Plant_ID_others,:) = 0;
for j = 1:10
%     [m,n]=find(unitmin<=j & B_utilize_trans_storage5(:,j)<=130.85);
    [m,n]=find(unitmin<=j);
    Ph_(:,j) = reshape(sum(sum(etrans_cou1_num(m,:,:),2)),[192,1]);
end

etrans_cou1_num2 = etrans_cou1_num;
etrans_con = reshape(sum(etrans_cou1_num2,2),[size(unitmin,1),192]); % plant num*192,耗电量

%% For country which achieve C neutrality before 2040
load('H:\global-PV-wind\Data\P2030_others_xz.mat') 
% 1.Biomass; 2.Geothermal; 3.Hydro; 4.Nuclear; 5.Ocean;
load('H:\global-PV-wind\Data\P2019_others_xz.mat')
P2025_others = P2030_others-(P2030_others-P2019_others)/11*5;
% load('H:\global-PV-wind\ANS\Ph_country2025_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2025_IIASA_C1toC8_median.mat');
PD2025 = Ph_country2025_IIASA_2deg;
PD(:,1) = PD2025;
P_others(:,1) = sum(P2025_others,2);
% load('H:\global-PV-wind\ANS\Ph_country2030_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2030_IIASA_C1toC8_median.mat');
PD2030 = Ph_country2030_IIASA_2deg;
PD(:,2) = PD2030;
load('H:\global-PV-wind\Data\P2030_others_xz.mat') 
P_others(:,2) = sum(P2030_others,2);
% load('H:\global-PV-wind\ANS\Ph_country2040_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2040_IIASA_C1toC8_median.mat');
PD2040 = Ph_country2040_IIASA_2deg;
PD(:,4) = PD2040;
load('H:\global-PV-wind\Data\P2040_others_xz.mat') 
P_others(:,4) = sum(P2040_others,2);
% load('H:\global-PV-wind\ANS\Ph_country2050_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2050_IIASA_C1toC8_median.mat');
PD2050 = Ph_country2050_IIASA_2deg;
PD(:,6) = PD2050;
load('H:\global-PV-wind\Data\P2050_others_xz.mat') 
P_others(:,6) = sum(P2050_others,2);
% load('H:\global-PV-wind\ANS\Ph_country2060_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2060_IIASA_C1toC8_median.mat');
PD2060 = Ph_country2060_IIASA_2deg;
PD(:,8) = PD2060;
PD(PD<0)=0;
load('H:\global-PV-wind\Data\P2060_others_xz.mat')  % P2060_others
P_others(:,8) = sum(P2060_others,2);
% load('H:\global-PV-wind\ANS\Ph_country2070_IIASA_2deg_max.mat');
load('H:\global-PV-wind\Data\Ph_country2070_IIASA_C1toC8_median.mat');
PD2070 = Ph_country2070_IIASA_2deg;
PD(:,10) = PD2070;
load('H:\global-PV-wind\Data\P2070_others_xz.mat')  % P2070_others
P_others(:,10) = sum(P2070_others,2);

for i = [3 5 7 9]
    PD(:,i) = (PD(:,i-1)+PD(:,i+1))/2;
    P_others(:,i) = (P_others(:,i-1)+P_others(:,i+1))/2;
end

% 2065年
[m,n]=find(Ph_(:,9)+P_others(:,9)<PD(:,9)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,9)-(P_others(:,9));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>9 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=9;
    else
        [m3,n]=find(unitmin>9 & etrans_con(:,m(i))~=0);
        unitmin(m3)=9;
%         P_others(m(i),4) = PD(m(i),4)-sum(etrans_con(:,m(i)));
    end  
end


% 2060年
[m,n]=find(Ph_(:,8)+P_others(:,8)<PD(:,8)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,8)-(P_others(:,8));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>8 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=8;
    else
        [m3,n]=find(unitmin>8 & etrans_con(:,m(i))~=0);
        unitmin(m3)=8;
%         P_others(m(i),4) = PD(m(i),4)-sum(etrans_con(:,m(i)));
    end  
end


% 2055年
[m,n]=find(Ph_(:,7)+P_others(:,7)<PD(:,7)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,7)-(P_others(:,7));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>7 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=7;
    else
        [m3,n]=find(unitmin>7 & etrans_con(:,m(i))~=0);
        unitmin(m3)=7;
%         P_others(m(i),4) = PD(m(i),4)-sum(etrans_con(:,m(i)));
    end  
end

% 2050年
[m,n]=find(Ph_(:,6)+P_others(:,6)<PD(:,6)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,6)-(P_others(:,6));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>6 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=6;
    else
        [m3,n]=find(unitmin>6 & etrans_con(:,m(i))~=0);
        unitmin(m3)=6;
%         P_others(m(i),4) = PD(m(i),4)-sum(etrans_con(:,m(i)));
    end  
end

% 2045年
[m,n]=find(Ph_(:,5)+P_others(:,5)<PD(:,5)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,5)-(P_others(:,5));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>5 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=5;
    else
        [m3,n]=find(unitmin>5 & etrans_con(:,m(i))~=0);
        unitmin(m3)=5;
    end  
end

% 2040年
[m,n]=find(Ph_(:,4)+P_others(:,4)<PD(:,4)); % 不满足需电量的国家,但是2060年的可以全部满足
country_needadjust = m;
B = PD(:,4)-(P_others(:,4));
for i = 1:size(m,1)
    [m4,n4]=find(cumsum(etrans_con(:,m(i)))>=B(m(i)));
    if ~isempty(m4)
        [m3,n]=find(unitmin>4 & etrans_con(:,m(i))~=0);
        m3(m3>m4(1)) = [];
        unitmin(m3)=4;
    else
        [m3,n]=find(unitmin>4 & etrans_con(:,m(i))~=0);
        unitmin(m3)=4;
    end  
end


save('H:\global-PV-wind\ANS\P_others2040_8_2s.mat','P_others'); %
save('H:\global-PV-wind\ANS\unitmin2040_8_2s_1.mat','unitmin'); %
save('H:\global-PV-wind\ANS\country_needadjust2040_8_2s.mat','country_needadjust'); %

