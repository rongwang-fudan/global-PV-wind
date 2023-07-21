tic
clear;
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
load('H:\Global PV and wind\ANS\unitmin_2030_2040.mat')  % 
MAC_CCS_p = unitmin*0;
for i = 1:10
    [m,n]=find(unitmin==i);
    MAC_CCS_p(m,1) =MAC_CCS(i);
end
clear MAC_CCS

%% 2025
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2025_2040.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2025xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2025_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2025xz-r_power_coop_notech_nopowertrans2030);
clear r_power_coop_notech_nopowertrans2030
Inv_PWtoCCS2025 = sum(repmat(Ph_PW(:,1),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2025_2040.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2025 = Inv_pvwind_mineral2030;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2025_2040_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2025_total = Inv_pvwind_mineral2030_total;
clear Inv_pvwind_mineral2030
clear r_power_coop2025xz
clear r_surplus

%% 2030, technology diffusion exert no influence
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2030_2040.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2030xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2030_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2030xz-r_power_coop_notech_nopowertrans2030);
clear r_power_coop_notech_nopowertrans2030
Inv_PWtoCCS2030 = sum(repmat(Ph_PW(:,2),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2030_2040.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2030 = Inv_pvwind_mineral2030;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2030_2040_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2030_total = Inv_pvwind_mineral2030_total;
clear Inv_pvwind_mineral2030
clear r_power_coop2030xz
clear r_surplus

%% 2035
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2035_2040.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2035xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2035_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2035xz-r_power_coop_notech_nopowertrans2030);
clear r_power_coop_notech_nopowertrans2030
Inv_PWtoCCS2035 = sum(repmat(Ph_PW(:,3),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2035_2040.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2035 = Inv_pvwind_mineral2030;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2035_2040_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2035_total = Inv_pvwind_mineral2030_total;
clear Inv_pvwind_mineral2030
clear r_power_coop2035xz
clear r_surplus

%% 2040
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2040.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2040xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2040_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2040xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2040 = sum(repmat(Ph_PW(:,4),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2040.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2040 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2040_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2040_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
clear r_power_coop2040xz
clear r_surplus


%% 2045
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2045.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2045xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2045_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2045xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2045 = sum(repmat(Ph_PW(:,5),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2045.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2045 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2045_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2045_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
clear r_power_coop2045xz
clear r_surplus

%% 2050
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2050.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2050xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2050_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2050xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2050 = sum(repmat(Ph_PW(:,6),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2050.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2050 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2050_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2050_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
clear r_power_coop2050xz
clear r_surplus

%% 2055
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2055.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2055xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2055_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2055xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2055 = sum(repmat(Ph_PW(:,7),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2055.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2055 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2055_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2055_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2055
clear r_power_coop2055xz
clear r_surplus

%% 2060
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2060.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2060xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2060_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2060xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2060 = sum(repmat(Ph_PW(:,8),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2060.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2060 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2060_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2060_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
clear r_power_coop2060xz
clear r_surplus


%% 2065
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2065.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2065xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2065_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2065xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2065 = sum(repmat(Ph_PW(:,9),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2065.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2065 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2065_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2065_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
clear r_power_coop2065xz
clear r_surplus

%% 2070
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2070.mat')  % r_power_coop_notech2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\r_power_coop_2070xz.mat')  % r_power_coop2040 % plant*connection(0-max)
load('H:\Global PV and wind\ANS\EF_2070_2040.mat')  % EF, t CO2/kWh
r_surplus = (r_power_coop2070xz-r_power_coop_notech_nopowertrans2040);
clear r_power_coop_notech_nopowertrans2040
Inv_PWtoCCS2070 = sum(repmat(Ph_PW(:,10),[1,size(r_surplus,2)]).*r_surplus.*(-repmat(EF,[1,size(r_surplus,2)])).*(repmat(MAC_CCS_p,[1,size(r_surplus,2)])));
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2070.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2070 = Inv_pvwind_mineral2040;
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_2070_total.mat')  % Inv_pvwind_mineral_powertrans2040
Inv_pvwind2070_total = Inv_pvwind_mineral2040_total;
clear Inv_pvwind_mineral2040
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
Inv_PWtoCCS_dis = [Inv_PWtoCCS2025' Inv_PWtoCCS2030' Inv_PWtoCCS2035' Inv_PWtoCCS2040' Inv_PWtoCCS2045' Inv_PWtoCCS2050' Inv_PWtoCCS2055' Inv_PWtoCCS2060' Inv_PWtoCCS2065' Inv_PWtoCCS2070'].*repmat(rr2',[size(Inv_pvwind2070,1),1]);
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

save('H:\Global PV and wind\ANS\Inv_PWall_dis2a_6_2040xz.mat','Inv_PWall_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_6_2040xz.mat','Inv_PWall_total_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_6_2040xz.mat','Inv_PWtoCCS_dis2a','-v7.3')  %

