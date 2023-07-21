tic
clear;
load('H:\Global PV and wind\Data\GDP_Country_2020_Nominal.mat')  % billion
% D:\FUdan\solar\world\GDP_Global\GDP by country.xlsx
load('H:\Global PV and wind\ANS\Ph_PVwind5_ef_cou_p_2040_8_2s_6_2a.mat')  % TWh/yr, Ph_PVwind5_ef_cou
load('H:\Global PV and wind\ANS\m_2040_8_2s_6_2a.mat')  % m,由于没有带来更大的收益在2021-2030年被剔除的电厂
Ph_PVwind5_ef_cou_p(m,:,:) = [];
for i = 1:10
    Ph_PW(:,i) = sum(sum(Ph_PVwind5_ef_cou_p(:,:,i),3),2);
end
load('H:\Global PV and wind\ANS\m_2040_8_2s_6_2a.mat')  % m,由于没有带来更大的收益在2021-2030年被剔除的电厂
m_cho = m;
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_8_2s_6_2a.mat')  % million $,Inv_CCS2070
load('H:\Global PV and wind\ANS\Inv_others2040_8_2s_6_2a.mat')  % million $,Inv_others
load('H:\Global PV and wind\ANS\Inv_others_e2040_8_2s_6_2a.mat')  % million $,% million $,抵消的fossil fuel的成本
Inv_others_e = -Inv_others_e;
Inv_CCS0 = Inv_CCS2070_c(1,:);
Inv_others0 = Inv_others(1,:);
Inv_others_e0 = Inv_others_e(1,:);
Inv_others_c0 = Inv_others0+Inv_others_e0;
clear Inv_others0
clear Inv_others_e0
Inv_CCS_dif = diff(Inv_CCS2070_c);
Inv_others_dif = diff(Inv_others);
Inv_others_e_dif = diff(Inv_others_e);
clear Inv_CCS2070_c
clear Inv_others
clear Inv_others_e
Inv_CCS_dif(m_cho,:) = []; % 增加此电厂导致CCS净成本减少量
Inv_others_dif(m_cho,:) = [];
Inv_others_e_dif(m_cho,:) = [];
Inv_others_c_dif = Inv_others_dif+Inv_others_e_dif; % 增加此电厂导致other renewbales净成本减少量
clear Inv_others_dif
clear Inv_others_e_dif

% Cost_CCS = (Inv_CCS0+sum(Inv_CCS_dif))'/10^6; % trillion $
% Cost_others = (Inv_others0+sum(Inv_others_dif))'/10^6; % trillion $
% Rev_others = (Inv_others_e0+sum(Inv_others_e_dif))'/10^6; % trillion $

load('H:\Global PV and wind\ANS\lrr_CCS2040_8_2s_6_2a.mat')  % lrr_CCS
MAC_CCS = 130.85*lrr_CCS; % MAC of CCS in 2020s,2030s,2040s,2050s,2060s
clear lrr_CCS
% load('H:\Global PV and wind\ANS\EF_2030_2040CaseA.mat')  % EF, t CO2/kWh
load('H:\Global PV and wind\ANS\EF_2030_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
load('H:\Global PV and wind\ANS\unitmin_2030_2040.mat')  % 
MAC_CCS_p = unitmin*0;
for i = 1:10
    [m,n]=find(unitmin==i);
    MAC_CCS_p(m,1) =MAC_CCS(i);
end
clear MAC_CCS

load('H:\Global PV and wind\ANS\powerunit_country_IX_IX_8_2040_2s_2070_test6xz_CaseA.mat'); % powerunit_country_IX_IX
load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat')
powerunit_country_IX_IX = powerunit_country_IX_IX(Plant_ID_IX,:);

%% 2025, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2025_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2025
load('H:\Global PV and wind\ANS\r_power_coop_2025xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2025_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,1).*sum(r_power_coop2025xz,2)-power_coop_notech_nopowertrans2025;
clear power_coop_notech_nopowertrans2025
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2025a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2025 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2025(i) = sum(Inv_PWtoCCS2025a(m));
end
clear Inv_PWtoCCS2025a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2025_2040CaseA.mat')  % inv_cou
Inv_pvwind2025 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2025_2040CaseA.mat')  % inv_cou
Inv_pvwind2025_total = inv_total_cou;
clear inv_cou
clear r_power_coop2025xz
clear r_surplus

%% 2030, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2030_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2030
load('H:\Global PV and wind\ANS\r_power_coop_2030xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2030_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,2).*sum(r_power_coop2030xz,2)-power_coop_notech_nopowertrans2030;
clear power_coop_notech_nopowertrans2030
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2030a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2030 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2030(i) = sum(Inv_PWtoCCS2030a(m));
end
clear Inv_PWtoCCS2030a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2030_2040CaseA.mat')  % inv_cou
Inv_pvwind2030 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2030_2040CaseA.mat')  % inv_cou
Inv_pvwind2030_total = inv_total_cou;
clear inv_cou
clear r_power_coop2030xz
clear r_surplus


%% 2035, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2035_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2035
load('H:\Global PV and wind\ANS\r_power_coop_2035xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2035_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,3).*sum(r_power_coop2035xz,2)-power_coop_notech_nopowertrans2035;
clear power_coop_notech_nopowertrans2035
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2035a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2035 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2035(i) = sum(Inv_PWtoCCS2035a(m));
end
clear Inv_PWtoCCS2035a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2035_2040CaseA.mat')  % inv_cou
Inv_pvwind2035 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2035_2040CaseA.mat')  % inv_cou
Inv_pvwind2035_total = inv_total_cou;
clear inv_cou
clear r_power_coop2035xz
clear r_surplus

%% 2040, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2040_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2040
load('H:\Global PV and wind\ANS\r_power_coop_2040xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2040_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,4).*sum(r_power_coop2040xz,2)-power_coop_notech_nopowertrans2040;
clear power_coop_notech_nopowertrans2040
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2040a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2040 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2040(i) = sum(Inv_PWtoCCS2040a(m));
end
clear Inv_PWtoCCS2040a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2040_2040CaseA.mat')  % inv_cou
Inv_pvwind2040 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2040_2040CaseA.mat')  % inv_cou
Inv_pvwind2040_total = inv_total_cou;
clear inv_cou
clear r_power_coop2040xz
clear r_surplus

%% 2045, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2045_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2045
load('H:\Global PV and wind\ANS\r_power_coop_2045xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2045_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,5).*sum(r_power_coop2045xz,2)-power_coop_notech_nopowertrans2045;
clear power_coop_notech_nopowertrans2045
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2045a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2045 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2045(i) = sum(Inv_PWtoCCS2045a(m));
end
clear Inv_PWtoCCS2045a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2045_2040CaseA.mat')  % inv_cou
Inv_pvwind2045 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2045_2040CaseA.mat')  % inv_cou
Inv_pvwind2045_total = inv_total_cou;
clear inv_cou
clear r_power_coop2045xz
clear r_surplus

%% 2050, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2050_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2050
load('H:\Global PV and wind\ANS\r_power_coop_2050xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2050_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,6).*sum(r_power_coop2050xz,2)-power_coop_notech_nopowertrans2050;
clear power_coop_notech_nopowertrans2050
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2050a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2050 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2050(i) = sum(Inv_PWtoCCS2050a(m));
end
clear Inv_PWtoCCS2050a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2050_2040CaseA.mat')  % inv_cou
Inv_pvwind2050 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2050_2040CaseA.mat')  % inv_cou
Inv_pvwind2050_total = inv_total_cou;
clear inv_cou
clear r_power_coop2050xz
clear r_surplus

%% 2055, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2055_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2055
load('H:\Global PV and wind\ANS\r_power_coop_2055xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2055_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,7).*sum(r_power_coop2055xz,2)-power_coop_notech_nopowertrans2055;
clear power_coop_notech_nopowertrans2055
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2055a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2055 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2055(i) = sum(Inv_PWtoCCS2055a(m));
end
clear Inv_PWtoCCS2055a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2055_2040CaseA.mat')  % inv_cou
Inv_pvwind2055 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2055_2040CaseA.mat')  % inv_cou
Inv_pvwind2055_total = inv_total_cou;
clear inv_cou
clear r_power_coop2055xz
clear r_surplus


%% 2060, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2060_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2060
load('H:\Global PV and wind\ANS\r_power_coop_2060xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2060_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,8).*sum(r_power_coop2060xz,2)-power_coop_notech_nopowertrans2060;
clear power_coop_notech_nopowertrans2060
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2060a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2060 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2060(i) = sum(Inv_PWtoCCS2060a(m));
end
clear Inv_PWtoCCS2060a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2060_2040CaseA.mat')  % inv_cou
Inv_pvwind2060 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2060_2040CaseA.mat')  % inv_cou
Inv_pvwind2060_total = inv_total_cou;
clear inv_cou
clear r_power_coop2060xz
clear r_surplus

%% 2065, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2065_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2065
load('H:\Global PV and wind\ANS\r_power_coop_2065xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2065_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,9).*sum(r_power_coop2065xz,2)-power_coop_notech_nopowertrans2065;
clear power_coop_notech_nopowertrans2065
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2065a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2065 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2065(i) = sum(Inv_PWtoCCS2065a(m));
end
clear Inv_PWtoCCS2065a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2065_2040CaseA.mat')  % inv_cou
Inv_pvwind2065 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2065_2040CaseA.mat')  % inv_cou
Inv_pvwind2065_total = inv_total_cou;
clear inv_cou
clear r_power_coop2065xz
clear r_surplus

%% 2070, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\power_coop_notech_nopowertrans2070_2040CaseA.mat')  % TWh/year,power_coop_notech_nopowertrans2070
load('H:\Global PV and wind\ANS\r_power_coop_2070xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2070_2040.mat')  % EF, t CO2/kWh
EF(find(isnan(EF)==1))=0;
p_surplus = Ph_PW(:,10).*sum(r_power_coop2070xz,2)-power_coop_notech_nopowertrans2070;
clear power_coop_notech_nopowertrans2070
p_surplus(p_surplus<0)=0;
Inv_PWtoCCS2070a = p_surplus.*(-EF).*MAC_CCS_p;
Inv_PWtoCCS2070 = zeros(192,1);
for i = 1:192
    [m,n]=find(powerunit_country_IX_IX==i);
    Inv_PWtoCCS2070(i) = sum(Inv_PWtoCCS2070a(m));
end
clear Inv_PWtoCCS2070a
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2070_2040CaseA.mat')  % inv_cou
Inv_pvwind2070 = inv_cou;
load('H:\Global PV and wind\ANS\inv_total_cou_PW_noconn_2070_2040CaseA.mat')  % inv_cou
Inv_pvwind2070_total = inv_total_cou;
clear inv_cou
clear r_power_coop2070xz
clear r_surplus

%%
discount = 0.03;
rr2(1,1) = 1/(1+discount)^(25-20);
rr2(2,1) = 1/(1+discount)^(30-20);
rr2(3,1) = 1/(1+discount)^(35-20);
rr2(4,1) = 1/(1+discount)^(40-20);
rr2(5,1) = 1/(1+discount)^(45-20);
rr2(6,1) = 1/(1+discount)^(50-20);
rr2(7,1) = 1/(1+discount)^(55-20);
rr2(8,1) = 1/(1+discount)^(60-20);
rr2(9,1) = 1/(1+discount)^(65-20);
rr2(10,1) = 1/(1+discount)^(70-20);
Inv_PWall_dis = [Inv_pvwind2025 Inv_pvwind2030 Inv_pvwind2035 Inv_pvwind2040 Inv_pvwind2045 Inv_pvwind2050 Inv_pvwind2055 Inv_pvwind2060 Inv_pvwind2065 Inv_pvwind2070].*repmat(rr2',[size(Inv_pvwind2070,1),1]);
Inv_PWtoCCS_dis = [Inv_PWtoCCS2025 Inv_PWtoCCS2030 Inv_PWtoCCS2035 Inv_PWtoCCS2040 Inv_PWtoCCS2045 Inv_PWtoCCS2050 Inv_PWtoCCS2055 Inv_PWtoCCS2060 Inv_PWtoCCS2065 Inv_PWtoCCS2070].*repmat(rr2',[size(Inv_pvwind2070,1),1]);
Inv_PWall_total_dis = [Inv_pvwind2025_total Inv_pvwind2030_total Inv_pvwind2035_total Inv_pvwind2040_total Inv_pvwind2045_total Inv_pvwind2050_total Inv_pvwind2055_total Inv_pvwind2060_total Inv_pvwind2065_total Inv_pvwind2070_total].*repmat(rr2',[size(Inv_pvwind2070_total,1),1]);

Inv_PWall_dis2 = [zeros(size(Inv_PWall_dis,1),1) Inv_PWall_dis];
Inv_PWall_total_dis2 = [zeros(size(Inv_PWall_total_dis,1),1) Inv_PWall_total_dis];
Inv_PWtoCCS_dis2 = [zeros(size(Inv_PWtoCCS_dis,1),1) Inv_PWtoCCS_dis];
for i = 1:10
    Inv_PWall_dis2a(:,i) = sum(Inv_PWall_dis2(:,i:i+1),2)/2;
    Inv_PWall_total_dis2a(:,i) = sum(Inv_PWall_total_dis2(:,i:i+1),2)/2;
    Inv_PWtoCCS_dis2a(:,i) = sum(Inv_PWtoCCS_dis2(:,i:i+1),2)/2;
end
% sum(Inv_CCSall_dis2a)+sum(Inv_othersall_dis2a)+sum(Inv_PWall_dis2a)+sum(Inv_PWtoCCS_dis2a)+Inv_CCS0_dis2a+Inv_others_c0_dis2a

save('H:\Global PV and wind\ANS\Inv_PWall_dis2a_2_2040CaseA_cou.mat','Inv_PWall_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_2_2040CaseA_cou.mat','Inv_PWall_total_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_2_2040CaseA_cou.mat','Inv_PWtoCCS_dis2a','-v7.3')  %

aa4 = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
sum(aa4(:,1))
