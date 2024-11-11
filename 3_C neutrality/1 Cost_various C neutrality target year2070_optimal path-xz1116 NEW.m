tic
clear;
% load('H:\global-PV-wind\ANS\P_others2040_8_2s.mat'); %
% P_others_xz = P_others;
MAC_others = [47.09,7.3,2.25,17.01,200];
% marginalcost_World.xlsx
% 1.Biomass; 2.Geothermal; 3.Hydro; 4.Nuclear; 5.Ocean;
IXX = [3,2,4,1,5];
MAC_others = MAC_others(IXX);

%%
load('H:\global-PV-wind\Data\Powerdmeand_country_coal.mat'); % TWh,行：1-192号国家，193行是全球 % 列：2000-2020
load('H:\global-PV-wind\Data\Powerdmeand_country_gas.mat'); % TWh,行：1-192号国家，193行是全球
load('H:\global-PV-wind\Data\Powerdmeand_country_oil.mat'); % TWh,行：1-192号国家，193行是全球
a = Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil;
[m,n]=find(a(:,21)==0);
clear a

r1=Powerdmeand_country_coal./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
r2=Powerdmeand_country_gas./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
r3=Powerdmeand_country_oil./(Powerdmeand_country_coal+Powerdmeand_country_gas+Powerdmeand_country_oil);
clear Powerdmeand_country_coal
clear Powerdmeand_country_gas
clear Powerdmeand_country_oil

r1(m,21)=r1(193,21);
r2(m,21)=r2(193,21);
r3(m,21)=r3(193,21);
r = r1(:,21)+r2(:,21)+r3(:,21);

rr(:,1) = r1(1:192,21);
rr(:,2) = r2(1:192,21);
rr(:,3) = r3(1:192,21);
% 2020年各国coal gas oil发电占fossil发电的比例
clear r
clear r1
clear r2
clear r3
load('H:\global-PV-wind\Data\EF_coal.mat'); % kg CO2/kWh
load('H:\global-PV-wind\Data\EF_gas.mat'); % kg CO2/kWh
load('H:\global-PV-wind\Data\EF_oil.mat'); % kg CO2/kWh
EF = [EF_coal EF_gas EF_oil];
clear EF_coal
clear EF_gas
clear EF_oil
EF_mean = sum(rr.*EF,2);
load('H:\global-PV-wind\Data\Price_coal.mat'); % USD/kWh, 第一列是 均值，第二列 std
load('H:\global-PV-wind\Data\Price_gas.mat'); % USD/kWh, 第一列是 均值，第二列 std
load('H:\global-PV-wind\Data\Price_oil.mat'); % USD/kWh, 第一列是 均值，第二列 std
Price = [Price_coal(:,1) Price_gas(:,1) Price_oil(:,1)];
EP_mean = sum(rr.*Price,2);
clear Price_coal
clear Price_gas
clear Price_oil
clear Price
save('H:\global-PV-wind\ANS\EF_mean_country.mat','EF_mean','-v7.3')  %  kg CO2/kWh
save('H:\global-PV-wind\ANS\EP_mean_country.mat','EP_mean','-v7.3')  %  USD/kWh

%% C neutrality target year of 2040
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2020s_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,1) = B_utilize_trans_storage; % 2025年MAC
B_utilize_trans_storage5(:,2) = B_utilize_trans_storage; % 2030年MAC
B_utilize_trans_storage5(:,3) = B_utilize_trans_storage; % 2035年MAC
B_utilize_trans_storage5(:,4) = B_utilize_trans_storage; % 2040年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2045_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,5) = B_utilize_trans_storage; % 2045年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2050_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,6) = B_utilize_trans_storage; % 2050年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2055_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,7) = B_utilize_trans_storage; % 2060年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2060_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,8) = B_utilize_trans_storage; % 2060年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2065_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,9) = B_utilize_trans_storage; % 2060年MAC
load('H:\global-PV-wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2070_2s_2070_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,10) = B_utilize_trans_storage; % 2060年MAC
load('H:\global-PV-wind\ANS\unitmin_global_IX_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070xz.mat');
load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_8_2070.mat')% etrans_cou1_num, 实际有效发电量, TWh/year
etrans_t = sum(sum(etrans_cou1_num,3),2);
load('H:\global-PV-wind\ANS\powerunit_country_IX_IX_8_2070_2s_2070_test6xz.mat'); % powerunit_country_IX_IX

load('H:\global-PV-wind\ANS\EP_8_2070_2s_2070_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2070_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2070_test6xz.mat')
CO2_year_utilize_trans_storage(:,10) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,10) = emissionfactor ; %
EP5(:,10) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2065_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2065_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2065_test6xz.mat')
CO2_year_utilize_trans_storage(:,9) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,9) = emissionfactor ; %
EP5(:,9) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2060_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2060_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2060_test6xz.mat')
CO2_year_utilize_trans_storage(:,8) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,8) = emissionfactor ; %
EP5(:,8) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2055_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2055_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2055_test6xz.mat')
CO2_year_utilize_trans_storage(:,7) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,7) = emissionfactor ; %
EP5(:,7) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2050_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2050_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2050_test6xz.mat')
CO2_year_utilize_trans_storage(:,6) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,6) = emissionfactor ; %
EP5(:,6) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2045_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2045_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2045_test6xz.mat')
CO2_year_utilize_trans_storage(:,5) = -phhall_all2.*emissionfactor ; % Mt CO2/y
EF5(:,5) = emissionfactor ; %
EP5(:,5) = EP ; %
load('H:\global-PV-wind\ANS\EP_8_2070_2s_2020s_test6xz.mat')  % EP，USD/kWh
load('H:\global-PV-wind\ANS\emissionfactor_8_2070_2s_2020s_test6xz.mat')
load('H:\global-PV-wind\ANS\phhall_all2_8_2070_2s_2020s_test6xz.mat')
CO2_year_utilize_trans_storage(:,1:4) = repmat(-phhall_all2.*emissionfactor,[1,4]); % Mt CO2/y
EF5(:,1:4) = repmat(emissionfactor,[1,4]) ; %
EP5(:,1:4) = repmat(EP,[1,4]); %
clear emissionfactor
clear EP

Inv = B_utilize_trans_storage5.*(-CO2_year_utilize_trans_storage)+repmat(phhall_all2,1,10).*EP5*10^3; % million USD/yr
Inv_noe = B_utilize_trans_storage5.*(-CO2_year_utilize_trans_storage); % million USD/yr
Inv(find(isnan(Inv)==1))=0;
Inv_noe(find(isnan(Inv_noe)==1))=0;
[B,IX] = sort(B_utilize_trans_storage5(:,10));
phhall_all2_IX70 = phhall_all2(IX);
B_utilize_trans_storage_IX70 = B_utilize_trans_storage5(IX,:);
unitmin_IX70 = unitmin(IX);
EP5_IX70 = EP5(IX,:);
EF5_IX70 = EF5(IX,:);
Inv_IX70 = Inv(IX,:);
Inv_noe_IX70 = Inv_noe(IX,:);
powerunit_country_IX70 = powerunit_country_IX_IX(IX);
CO2_year_utilize_trans_storage_IX70 = CO2_year_utilize_trans_storage(IX,:);
etrans_t_IX70 =etrans_t(IX);
etrans_cou1_num_IX70 =etrans_cou1_num(IX,:,:);
clear phhall_all2
clear B_utilize_trans_storage
clear Inv
clear powerunit_country_IX_IX
clear CO2_year_utilize_trans_storage
clear etrans_t
clear etrans_cou1_num
clear unitmin
[m,n]=find(B_utilize_trans_storage_IX70(:,10)==Inf);
phhall_all2_IX70(m) = [];
B_utilize_trans_storage_IX70(m,:) = [];
unitmin_IX70(m) = [];
EP5_IX70(m,:) = [];
EF5_IX70(m,:) = [];
Inv_IX70(m,:) = [];
Inv_noe_IX70(m,:) = [];
powerunit_country_IX70(m) = [];
CO2_year_utilize_trans_storage_IX70(m,:) = [];
etrans_t_IX70(m) = [];
etrans_cou1_num_IX70(m,:,:) = [];
save('H:\global-PV-wind\ANS\unitmin_IX70_2070_8_2s_6_2a_all.mat','unitmin_IX70','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\powerunit_country_IX70_2070_8_2s_6_2a_all.mat','powerunit_country_IX70','-v7.3')  % TWh/yr
pcon_2070C = reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]);
save('H:\global-PV-wind\ANS\pcon_2070C.mat','pcon_2070C','-v7.3')  % TWh/yr

%% 2070
pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2070_IIASA_C1toC8_median.mat');
r = pcon./repmat(Ph_country2070_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2070_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,10) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,10) = p;
r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,10).*r;
Inv_IX70_2=Inv_IX70(:,10).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,10).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,10) = [0;phhall_all2_IX70_2];

Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,10) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];

r = Ph_cou_per./repmat(sum(Ph_cou_per),192,1);
r(find(isnan(r)==1))=0;
CO2_pvwind_cou = repmat([0;-(CO2_year_utilize_trans_storage_IX70_2)/1000]',192,1).*r;
CO2a =cumsum(CO2_pvwind_cou,2);

% other renewable
load('H:\global-PV-wind\Data\P2070_others_xz.mat')  % P2070_others
% 1.Biomass; 2.Geothermal; 3.Hydro; 4.Nuclear; 5.Ocean;
IXX = [3,2,4,1,5];
P2070_others = P2070_others(:,IXX);
% 3.Hydro; 2.Geothermal; 4.Nuclear; 1.Biomass; 5.Ocean;
P_cum = (cumsum(P2070_others'))';
others(:,1) = min(sum(P2070_others,2),Ph_country2070_IIASA_2deg);
% 没有PV和wind时的other renewable的发电量
[m,n] = find(sum(P2070_others,2)-others(:,1)==0);
others_type(m,:) = P2070_others(m,:);
[m,n] = find(sum(P2070_others,2)-others(:,1)>0);
if ~isempty(m)
    for i = 1:size(m,1)
        [m2,n2] = find(P_cum(m(i),:)<=Ph_country2070_IIASA_2deg(m(i)));
        others_type(m(i),n2) = P2070_others(m(i),n2);
        others_type(m(i),n2(end)+1) = Ph_country2070_IIASA_2deg(m(i))-P_cum(m(i),n2(end));
    end
end
p_others2070(1,10) = sum(sum(others_type));
save('H:\global-PV-wind\ANS\others_type0_2070Cneutrality.mat','others_type','-v7.3')  % TWh/yr

Aba_others_cou2 = others_type.*repmat(EF_mean,1,5)/1000; % Gt CO2/yr
Aba_others_cou2_2070_0 = sum(Aba_others_cou2,2);
Inv_others2070_t0(:,1) = sum(repmat(MAC_others,192,1).*Aba_others_cou2*1000,2)+sum(others_type,2).*EP_mean*10^3; % million USD
Inv_others2070_e0(:,1) = sum(others_type,2).*EP_mean*10^3; % million USD

Inv_others2070(1,10) = 0; % million USD, change value
Inv_others2070_e(1,10) = 0; % million USD
Inv_others2070_c(:,1) = zeros(192,1); % million USD
Inv_others2070_c_e(:,1) = zeros(192,1); % million USD
p_country = zeros(192,1);
others_type_0123 = others_type;
for i = 1:1:size(B_utilize_trans_storage_IX70,1)
    p_country(:,1) = (pcon(i,:))';
    [cou_all,n]=find((p(i,:))'~=0);
    [cou_allno,n]=find((p(i,:))'==0);
    if isempty(cou_all)
        others(:,i+1) = others(:,i);
        Inv_others2070(i+1,10) = 0;
        Inv_others2070_e(i+1,10) = 0;
        p_others2070(i+1,10) = sum(sum(others(:,i+1)));
        Inv_others2070_c(:,i+1) = 0; % million USD
        Inv_others2070_c_e(:,i+1) = 0; % million USD
    else
        others(cou_allno,i+1) = others(cou_allno,i);
        Inv_others2070_c(cou_allno,i+1) = 0; % million USD
        Inv_others2070_c_e(cou_allno,i+1) = 0; % million USD
        for jj = 1:size(cou_all,1)
            cou = cou_all(jj);
            if p_country(cou,1)+others(cou,i)<=Ph_country2070_IIASA_2deg(cou,1)
                others(cou,i+1) = others(cou,i);
                Inv_others2070_c(cou,i+1) = 0; % million USD
                Inv_others2070_c_e(cou,i+1) = 0; % million USD
            else
                others(cou,i+1) = (Ph_country2070_IIASA_2deg(cou,1)-p_country(cou,1));
                [m2,n2] = find(P_cum(cou,:)<=others(cou,i+1));
                others_type0 = others_type(cou,:);
                if ~isempty(m2)
                    others_type(cou,n2) = P2070_others(cou,n2);
                    others_type(cou,n2(end)+1) = P_cum(cou,n2(end)+1)-others(cou,i+1);
                    [m3,n3] = find(P_cum(cou,:)>others(cou,i+1));
                    others_type(cou,n3) = 0;
                else
                    others_type(cou,1) = others(cou,i+1);
                    others_type(cou,2:5) = 0;
                end
                Aba_others_cou2 = (others_type(cou,:)-others_type0).*EF5_IX70(i,10)/1000; % Gt CO2/yr
                Inv_others2070_c(cou,i+1) = sum(MAC_others.*Aba_others_cou2*1000,2)+sum(others_type(cou,:)-others_type0).*EP5_IX70(i,10)*10^3; % million USD
                Inv_others2070_c_e(cou,i+1) = sum(others_type(cou,:)-others_type0).*EP5_IX70(i,10)*10^3; % million USD
            end
        end
        p_others2070(i+1,10) = sum(sum(others(:,i+1)));
        Inv_others2070(i+1,10) = sum(Inv_others2070_c(:,i+1)); % 变化
        Inv_others2070_e(i+1,10) = sum(Inv_others2070_c_e(:,i+1));
    end
    others_type1a(:,i+1) = others_type_0123(:,1)-others_type(:,1);
    others_type2a(:,i+1) = others_type_0123(:,2)-others_type(:,2);
    others_type3a(:,i+1) = others_type_0123(:,3)-others_type(:,3);
    others_type4a(:,i+1) = others_type_0123(:,4)-others_type(:,4);
    others_type5a(:,i+1) = others_type_0123(:,5)-others_type(:,5);
    others_type_0123 = others_type;
end
clear Inv_others2070_c_c
clear Inv_others2070_c
clear Inv_others2070_c_e
save('H:\global-PV-wind\ANS\others_type1a_2070Cneutrality.mat','others_type1a','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_type2a_2070Cneutrality.mat','others_type2a','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_type3a_2070Cneutrality.mat','others_type3a','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_type4a_2070Cneutrality.mat','others_type4a','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_type5a_2070Cneutrality.mat','others_type5a','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_2070Cneutrality.mat','others','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\others_type_2070Cneutrality.mat','others_type','-v7.3')  % TWh/yr
% % 1.Biomass; 2.Geothermal; 3.Hydro; 4.Nuclear; 5.Ocean;
% 3.Hydro; 2.Geothermal; 4.Nuclear; 1.Biomass; 5.Ocean;

Aba_others_cou2 = others_type.*repmat(EF_mean,1,5)/1000; % Gt CO2/yr
Inv_others2070_end = sum(sum(repmat(MAC_others,192,1).*Aba_others_cou2*1000))+sum(sum(others_type,2).*EP_mean*10^3); % million USD
Inv_others2070_e_end = sum(sum(others_type,2).*EP_mean*10^3); % million USD
Inv_others2070(:,10) = flip(cumsum([Inv_others2070_end;-flip(Inv_others2070(2:end,10))]));
Inv_others2070_e(:,10) = flip(cumsum([Inv_others2070_e_end;-flip(Inv_others2070_e(2:end,10))]));
powerdemand = repmat(Ph_country2070_IIASA_2deg,1,size(others,2))-others;

% CCS
p_CCS_cou = powerdemand-Ph_cou_per;
save('H:\global-PV-wind\ANS\p_CCS_cou_2070Cneutrality.mat','p_CCS_cou','-v7.3')  % TWh/yr
p_CCS_cou(p_CCS_cou<0)=0;
Ph_CCS5(:,10) = (sum(p_CCS_cou))';
pccs = [0;diff(Ph_CCS5(:,10))];
Aba_CCS_cou_dif = sum((diff(p_CCS_cou')).*repmat(EF5_IX70(:,10),1,192)/1000)'; % Gt CO2/yr
Aba_CCS_cou_end = p_CCS_cou(:,end).*EF_mean/1000; % Gt CO2/yr
Aba_CCS_cou2(:,1) = Aba_CCS_cou_end-Aba_CCS_cou_dif;
Aba_CCS_cou2(Aba_CCS_cou2<0)=0;
Aba_CCS_cou2_2070 = Aba_CCS_cou2;

for i = 1:1:size(B_utilize_trans_storage_IX70,1)+1
    if i ==1
        Inv_IX2 = 0;
        Inv_noe_IX2 = 0;
        Inv_CCS2070(i,10) = 0; % million USD，变化值
        Inv_CCS2070_c(i,10) = 0; % million USD
        Inv_CCS2070_e(i,10) = 0; % million USD
    else
        Inv_IX2 = Inv_IX70_2(1:i-1);
        Inv_noe_IX2 = Inv_noe_IX70_2(1:i-1);
        Aba_CCS = pccs(i).*EF5_IX70(i-1,10)/1000;
        Inv_CCS2070(i,10) = sum(130.85*Aba_CCS*1000)+pccs(i).*EP5_IX70(i-1,10)*10^3; % million USD，变化值
        Inv_CCS2070_c(i,10) = sum(130.85*Aba_CCS*1000); % million USD
        Inv_CCS2070_e(i,10) = pccs(i).*EP5_IX70(i-1,10)*10^3; % million USD
    end
    Inv_PVwind2070(i,10) = sum(Inv_IX2);
    Inv_PVwind2070_noe(i,10) = sum(Inv_noe_IX2);
end

Aba_others_cou2 = p_CCS_cou(:,end).*EF_mean/1000; % Gt CO2/yr
% Inv_CCS2070_end = sum(sum(repmat(MAC_others,192,1).*Aba_others_cou2*1000))+sum(p_CCS_cou(:,end).*EP_mean*10^3); % million USD
% Inv_CCS2070_c_end = sum(sum(repmat(MAC_others,192,1).*Aba_others_cou2*1000)); % million USD
Inv_CCS2070_end = sum(130.85.*Aba_others_cou2*1000)+sum(p_CCS_cou(:,end).*EP_mean*10^3); % million USD
Inv_CCS2070_c_end = sum(130.85.*Aba_others_cou2*1000); % million USD
Inv_CCS2070_e_end = sum(p_CCS_cou(:,end).*EP_mean*10^3); % million USD
Inv_CCS2070(:,10) = flip(cumsum([Inv_CCS2070_end;-flip(Inv_CCS2070(2:end,10))]));
Inv_CCS2070_c(:,10) = flip(cumsum([Inv_CCS2070_c_end;-flip(Inv_CCS2070_c(2:end,10))]));
Inv_CCS2070_e(:,10) = flip(cumsum([Inv_CCS2070_e_end;-flip(Inv_CCS2070_e(2:end,10))]));


%% 规划前期建设
Inv_CCS2070_2 = Inv_CCS2070_c;
Ph_CCS5_2 = Ph_CCS5;
LR_CCS = 0.0317;
P_renewable2020 = 5171.93599; % TWh/yr
P_CCS2020 = 879.5886979;
P_CCS = Ph_CCS5(:,10);
Cost_CCS = Inv_CCS2070_c(:,10); % million USD
LLRR = (log(1- LR_CCS)/log(2));
for nnn = 1%:1:size(P_CCS,1)
    Cost_CCSnnn = Cost_CCS(nnn)/100;
    P_CCSnnn = P_CCS(nnn)/100;
    lrr = [1;(cumsum(ones(100,1)*P_CCSnnn)./P_CCS2020+1).^LLRR];
    Cost_CCSnnna = [0;cumsum(ones(100,1)*Cost_CCSnnn)];
    P_CCSnnna = [0;cumsum(ones(100,1)*P_CCSnnn)];
    costmin1=1e18;
    for p1=0:2:100
        Cost_xz(1,1) = Cost_CCSnnna(p1+1);
        for p2=0:2:(100-p1)
            Cost_xz(2,1) = Cost_CCSnnna(p2+1)*lrr(p1+1);
            for p3=0:2:(100-p1-p2)
                Cost_xz(3,1) = Cost_CCSnnna(p3+1)*lrr(p1+p2+1);
                for p4=0:2:(100-p1-p2-p3)
                    Cost_xz(4,1) = Cost_CCSnnna(p4+1)*lrr(p1+p2+p3+1);
                    for p5=0:2:(100-p1-p2-p3-p4)
                        Cost_xz(5,1) = Cost_CCSnnna(p5+1)*lrr(p1+p2+p3+p4+1);
                        for p6=0:2:(100-p1-p2-p3-p4-p5)
                            Cost_xz(6,1) = Cost_CCSnnna(p6+1)*lrr(p1+p2+p3+p4+p5+1);
                            for p7=0:2:(100-p1-p2-p3-p4-p5-p6)
                                Cost_xz(7,1) = Cost_CCSnnna(p7+1)*lrr(p1+p2+p3+p4+p5+p6+1);
                                for p8=0:2:(100-p1-p2-p3-p4-p5-p6-p7)
                                    Cost_xz(8,1) = Cost_CCSnnna(p8+1)*lrr(p1+p2+p3+p4+p5+p6+p7+1);
                                    for p9=0:2:(100-p1-p2-p3-p4-p5-p6-p7-p8)
                                        Cost_xz(9,1) = Cost_CCSnnna(p9+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+1);
                                        p10 = 100-p1-p2-p3-p4-p5-p6-p7-p8-p9;
                                        Cost_xz(10,1) = Cost_CCSnnna(p10+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+p9+1);
                                        costall = sum(Cost_xz);
                                        if costall<costmin1
                                            Inv_CCS2070_2(nnn,1:10) = (cumsum(Cost_xz))';
                                            mmn = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            Ph_CCS5_2(nnn,1:10) = P_CCSnnna(mmn+1)';
                                            prop_CCS = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            costmin1 = costall;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        p1
    end
end


for nnn = 2:1:size(P_CCS,1)
    Cost_CCSnnn = Cost_CCS(nnn)/100;
    P_CCSnnn = P_CCS(nnn)/100;
    lrr = [1;(cumsum(ones(100,1)*P_CCSnnn)./P_CCS2020+1).^LLRR];
    Cost_CCSnnna = [0;cumsum(ones(100,1)*Cost_CCSnnn)];
    P_CCSnnna = [0;cumsum(ones(100,1)*P_CCSnnn)];
    for p1=prop_CCS(1)
        Cost_xz(1,1) = Cost_CCSnnna(p1+1);
        for p2=prop_CCS(2)-prop_CCS(1)
            Cost_xz(2,1) = Cost_CCSnnna(p2+1)*lrr(p1+1);
            for p3=prop_CCS(3)-prop_CCS(2)
                Cost_xz(3,1) = Cost_CCSnnna(p3+1)*lrr(p1+p2+1);
                for p4=prop_CCS(4)-prop_CCS(3)
                    Cost_xz(4,1) = Cost_CCSnnna(p4+1)*lrr(p1+p2+p3+1);
                    for p5=prop_CCS(5)-prop_CCS(4)
                        Cost_xz(5,1) = Cost_CCSnnna(p5+1)*lrr(p1+p2+p3+p4+1);
                        for p6=prop_CCS(6)-prop_CCS(5)
                            Cost_xz(6,1) = Cost_CCSnnna(p6+1)*lrr(p1+p2+p3+p4+p5+1);
                            for p7=prop_CCS(7)-prop_CCS(6)
                                Cost_xz(7,1) = Cost_CCSnnna(p7+1)*lrr(p1+p2+p3+p4+p5+p6+1);
                                for p8=prop_CCS(8)-prop_CCS(7)
                                    Cost_xz(8,1) = Cost_CCSnnna(p8+1)*lrr(p1+p2+p3+p4+p5+p6+p7+1);
                                    for p9=prop_CCS(9)-prop_CCS(8)
                                        Cost_xz(9,1) = Cost_CCSnnna(p9+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+1);
                                        for p10=prop_CCS(10)-prop_CCS(9)
                                            Cost_xz(10,1) = Cost_CCSnnna(p10+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+p9+1);
                                            
                                            costall = sum(Cost_xz);
                                            lrr_CCS = [1;lrr(p1+1);lrr(p1+p2+1);lrr(p1+p2+p3+1);lrr(p1+p2+p3+p4+1);lrr(p1+p2+p3+p4+p5+1);lrr(p1+p2+p3+p4+p5+p6+1);lrr(p1+p2+p3+p4+p5+p6+p7+1);lrr(p1+p2+p3+p4+p5+p6+p7+p8+1);lrr(p1+p2+p3+p4+p5+p6+p7+p8+p9+1);lrr(end)];
                                            Inv_CCS2070_2(nnn,1:10) = (cumsum(Cost_xz))';
                                            mmn = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            Ph_CCS5_2(nnn,1:10) = P_CCSnnna(mmn+1)';
                                            
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
r_CCS_period = mmn;
save('H:\global-PV-wind\ANS\r_CCS_period_2070Cneutrality.mat','r_CCS_period','-v7.3')  % 
save('H:\global-PV-wind\ANS\lrr_CCS2070_8_2s_6_2a.mat','lrr_CCS','-v7.3')  % million $，+附带fossil fuel成本

Inv_others2070_2 = Inv_others2070;
p_others2070_2 = p_others2070;
LR_others = 0.05;
P_renewable2020 = 5171.93599; % TWh/yr
P_others = p_others2070(:,10);
Cost_others = Inv_others2070(:,10); % million USD
LLRR = (log(1- LR_others)/log(2));
for nnn = 1 % :1:size(P_others,1)
    Cost_CCSnnn = Cost_others(nnn)/100;
    P_CCSnnn = P_others(nnn)/100;
    lrr = [1;(cumsum(ones(100,1)*P_CCSnnn)./P_renewable2020+1).^LLRR];
    Cost_CCSnnna = [0;cumsum(ones(100,1)*Cost_CCSnnn)];
    P_CCSnnna = [0;cumsum(ones(100,1)*P_CCSnnn)];
    costmin1=1e18;
    for p1=0:2:100
        Cost_xz(1,1) = Cost_CCSnnna(p1+1);
        for p2=0:2:(100-p1)
            Cost_xz(2,1) = Cost_CCSnnna(p2+1)*lrr(p1+1);
            for p3=0:2:(100-p1-p2)
                Cost_xz(3,1) = Cost_CCSnnna(p3+1)*lrr(p1+p2+1);
                for p4=0:2:(100-p1-p2-p3)
                    Cost_xz(4,1) = Cost_CCSnnna(p4+1)*lrr(p1+p2+p3+1);
                    for p5=0:2:(100-p1-p2-p3-p4)
                        Cost_xz(5,1) = Cost_CCSnnna(p5+1)*lrr(p1+p2+p3+p4+1);
                        for p6=0:2:(100-p1-p2-p3-p4-p5)
                            Cost_xz(6,1) = Cost_CCSnnna(p6+1)*lrr(p1+p2+p3+p4+p5+1);
                            for p7=0:2:(100-p1-p2-p3-p4-p5-p6)
                                Cost_xz(7,1) = Cost_CCSnnna(p7+1)*lrr(p1+p2+p3+p4+p5+p6+1);
                                for p8=0:2:(100-p1-p2-p3-p4-p5-p6-p7)
                                    Cost_xz(8,1) = Cost_CCSnnna(p8+1)*lrr(p1+p2+p3+p4+p5+p6+p7+1);
                                    for p9=0:2:(100-p1-p2-p3-p4-p5-p6-p7-p8)
                                        Cost_xz(9,1) = Cost_CCSnnna(p9+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+1);
                                        p10 = 100-p1-p2-p3-p4-p5-p6-p7-p8-p9;
                                        Cost_xz(10,1) = Cost_CCSnnna(p10+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+p9+1);
                                        costall = sum(Cost_xz);
                                        
                                        if costall<costmin1
                                            Inv_others2070_2(nnn,1:10) = (cumsum(Cost_xz))';
                                            mmn = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            p_others2070_2(nnn,1:10) = P_CCSnnna(mmn+1)';
                                            prop_others = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            costmin1 = costall;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        p1
    end
end

for nnn = 2:1:size(P_others,1)
    Cost_CCSnnn = Cost_others(nnn)/100;
    P_CCSnnn = P_others(nnn)/100;
    lrr = [1;(cumsum(ones(100,1)*P_CCSnnn)./P_renewable2020+1).^LLRR];
    Cost_CCSnnna = [0;cumsum(ones(100,1)*Cost_CCSnnn)];
    P_CCSnnna = [0;cumsum(ones(100,1)*P_CCSnnn)];
    for p1=prop_others(1)
        Cost_xz(1,1) = Cost_CCSnnna(p1+1);
        for p2=prop_others(2)-prop_others(1)
            Cost_xz(2,1) = Cost_CCSnnna(p2+1)*lrr(p1+1);
            for p3=prop_others(3)-prop_others(2)
                Cost_xz(3,1) = Cost_CCSnnna(p3+1)*lrr(p1+p2+1);
                for p4=prop_others(4)-prop_others(3)
                    Cost_xz(4,1) = Cost_CCSnnna(p4+1)*lrr(p1+p2+p3+1);
                    for p5=prop_others(5)-prop_others(4)
                        Cost_xz(5,1) = Cost_CCSnnna(p5+1)*lrr(p1+p2+p3+p4+1);
                        for p6=prop_others(6)-prop_others(5)
                            Cost_xz(6,1) = Cost_CCSnnna(p6+1)*lrr(p1+p2+p3+p4+p5+1);
                            for p7=prop_others(7)-prop_others(6)
                                Cost_xz(7,1) = Cost_CCSnnna(p7+1)*lrr(p1+p2+p3+p4+p5+p6+1);
                                for p8=prop_others(8)-prop_others(7)
                                    Cost_xz(8,1) = Cost_CCSnnna(p8+1)*lrr(p1+p2+p3+p4+p5+p6+p7+1);
                                    for p9=prop_others(9)-prop_others(8)
                                        Cost_xz(9,1) = Cost_CCSnnna(p9+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+1);
                                        for p10=prop_others(10)-prop_others(9)
                                            Cost_xz(10,1) = Cost_CCSnnna(p10+1)*lrr(p1+p2+p3+p4+p5+p6+p7+p8+p9+1);
                                            
                                            costall = sum(Cost_xz);
                                            Inv_others2070_2(nnn,1:10) = (cumsum(Cost_xz))';
                                            mmn = [p1,p1+p2,p1+p2+p3,p1+p2+p3+p4,p1+p2+p3+p4+p5,p1+p2+p3+p4+p5+p6,p1+p2+p3+p4+p5+p6+p7,p1+p2+p3+p4+p5+p6+p7+p8,p1+p2+p3+p4+p5+p6+p7+p8+p9,100];
                                            p_others2070_2(nnn,1:10) = P_CCSnnna(mmn+1)';
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
r_others_period = mmn;
save('H:\global-PV-wind\ANS\r_others_period_2070Cneutrality.mat','r_others_period','-v7.3')  % 

Inv_CCS2070 = Inv_CCS2070_2+repmat(Inv_CCS2070_e(:,10),1,10).*(Ph_CCS5_2./repmat(Ph_CCS5_2(:,10),1,10));
bb = repmat(Inv_CCS2070_e(:,10),1,10).*(Ph_CCS5_2./repmat(Ph_CCS5_2(:,10),1,10));
rrc = Aba_CCS_cou2(:,1)./sum(Aba_CCS_cou2(:,1));
Inv_CCS2070_c_cou0(:,1:10) = repmat(Inv_CCS2070_2(1,1:10),192,1).*repmat(rrc,1,10); % +repmat(bb(1,:),192,1).*repmat(rrc,1,5);
save('H:\global-PV-wind\ANS\Inv_CCS2070_c_cou0_8_2s_6_2a.mat','Inv_CCS2070_c_cou0','-v7.3')  % million $，+附带fossil fuel成本
sum(Inv_CCS2070_c_cou0)-Inv_CCS2070_2(1,:)


Ph_CCS5 = Ph_CCS5_2;
Inv_others = Inv_others2070_2;
Inv_others_e = repmat(Inv_others2070_e(:,10),1,10).*(p_others2070_2./repmat(p_others2070_2(:,10),1,10));
Inv_others_c = Inv_others-Inv_others_e;

a = Inv_others2070_t0-Inv_others2070_e0;
rrc = a(:,1)./sum(a(:,1));
Inv_others2070_c_cou0(:,1:10) = repmat(Inv_others_c(1,1:10),192,1).*repmat(rrc,1,10); % +repmat(bb(1,:),192,1).*repmat(rrc,1,5);
Inv_others2070_cou0(:,1:10) = repmat(Inv_others(1,1:10),192,1).*repmat(rrc,1,10); % +repmat(bb(1,:),192,1).*repmat(rrc,1,5);
sum(Inv_others2070_c_cou0)-Inv_others_c(1,:)
save('H:\global-PV-wind\ANS\Inv_others2070_c_cou0_8_2s_6_2a.mat','Inv_others2070_c_cou0','-v7.3')  % million $，+附带fossil fuel成本
save('H:\global-PV-wind\ANS\Inv_others2070_cou0_8_2s_6_2a.mat','Inv_others2070_cou0','-v7.3')  % million $，+附带fossil fuel成本

Ph_others5 = p_others2070_2;
save('H:\global-PV-wind\ANS\Inv_others_e2070_8_2s_6_2a.mat','Inv_others_e','-v7.3')  % million $,抵消的fossil fuel的成本
save('H:\global-PV-wind\ANS\Inv_others2070_8_2s_6_2a.mat','Inv_others','-v7.3')  % million $
save('H:\global-PV-wind\ANS\Ph_others2070_8_2s_6_2a.mat','Ph_others5','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\Inv_CCS2070_8_2s_6_2a.mat','Inv_CCS2070','-v7.3')  % million $
save('H:\global-PV-wind\ANS\Ph_CCS2070_8_2s_6_2a.mat','Ph_CCS5','-v7.3')  % TWh/yr

Inv_CCS2070_c = Inv_CCS2070_2;
Inv_CCS2070_e = repmat(Inv_CCS2070_e(:,10),1,10).*(Ph_CCS5_2./repmat(Ph_CCS5_2(:,10),1,10))
save('H:\global-PV-wind\ANS\Inv_CCS2070_c_8_2s_6_2a.mat','Inv_CCS2070_c','-v7.3')  % million $,CCS投资
save('H:\global-PV-wind\ANS\Inv_CCS2070_e_8_2s_6_2a.mat','Inv_CCS2070_e','-v7.3')  % million $，CCS附带的Fossile fuel成本



%% 2065
[m,n]=find(unitmin_IX70>9);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2070_IIASA_C1toC8_median.mat');
load('H:\global-PV-wind\Data\Ph_country2060_IIASA_C1toC8_median.mat');
Ph_country2065_IIASA_2deg = Ph_country2060_IIASA_2deg/2+Ph_country2070_IIASA_2deg/2;
r = pcon./repmat(Ph_country2065_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2065_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,9) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,9) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,9).*r;
Inv_IX70_2=Inv_IX70(:,9).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,9).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,9) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,9) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,9) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,9) = [0;cumsum(Inv_noe_IX70_2)];


%% 2060
[m,n]=find(unitmin_IX70>8);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;
pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2060_IIASA_C1toC8_median.mat');
r = pcon./repmat(Ph_country2060_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2060_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,8) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,8) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,8).*r;
Inv_IX70_2=Inv_IX70(:,8).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,8).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,8) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,8) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,8) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,8) = [0;cumsum(Inv_noe_IX70_2)];

%% 2055
[m,n]=find(unitmin_IX70>7);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2060_IIASA_C1toC8_median.mat');
load('H:\global-PV-wind\Data\Ph_country2050_IIASA_C1toC8_median.mat');
Ph_country2055_IIASA_2deg = Ph_country2050_IIASA_2deg/2+Ph_country2060_IIASA_2deg/2;
r = pcon./repmat(Ph_country2055_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2055_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,7) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,7) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,7).*r;
Inv_IX70_2=Inv_IX70(:,7).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,7).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,7) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,7) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,7) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,7) = [0;cumsum(Inv_noe_IX70_2)];

%% 2050
[m,n]=find(unitmin_IX70>6);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2050_IIASA_C1toC8_median.mat');
r = pcon./repmat(Ph_country2050_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2050_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,6) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,6) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,6).*r;
Inv_IX70_2=Inv_IX70(:,6).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,6).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,6) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,6) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,6) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,6) = [0;cumsum(Inv_noe_IX70_2)];

%% 2045
[m,n]=find(unitmin_IX70>5);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2050_IIASA_C1toC8_median.mat');
load('H:\global-PV-wind\Data\Ph_country2040_IIASA_C1toC8_median.mat');
Ph_country2045_IIASA_2deg = Ph_country2040_IIASA_2deg/2+Ph_country2050_IIASA_2deg/2;
r = pcon./repmat(Ph_country2045_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2045_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,5) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,5) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,5).*r;
Inv_IX70_2=Inv_IX70(:,5).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,5).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,5) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,5) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,5) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,5) = [0;cumsum(Inv_noe_IX70_2)];

%% 2040
[m,n]=find(unitmin_IX70>4);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2040_IIASA_C1toC8_median.mat');
r = pcon./repmat(Ph_country2040_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2040_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,4) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,4) = p;

r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,4).*r;
Inv_IX70_2=Inv_IX70(:,4).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,4).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,4) = [0;phhall_all2_IX70_2];
Ph_cou_per = [zeros(192,1) (pcon)'];
CO2_cum_pvwind(:,4) = [0;-cumsum(CO2_year_utilize_trans_storage_IX70_2)/1000];
Inv_PVwind2070(:,4) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,4) = [0;cumsum(Inv_noe_IX70_2)];

%% 2035
[m,n]=find(unitmin_IX70>3);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2030_IIASA_C1toC8_median.mat'); 
load('H:\global-PV-wind\Data\Ph_country2040_IIASA_C1toC8_median.mat'); 
Ph_country2035_IIASA_2deg = Ph_country2030_IIASA_2deg/2+Ph_country2040_IIASA_2deg/2;
r = pcon./repmat(Ph_country2035_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2035_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,3) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,3) = p;
r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,3).*r;
Inv_IX70_2=Inv_IX70(:,3).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,3).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,3) = [0;phhall_all2_IX70_2];
Inv_PVwind2070(:,3) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,3) = [0;cumsum(Inv_noe_IX70_2)];

%% 2030
[m,n]=find(unitmin_IX70>2);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2030_IIASA_C1toC8_median.mat'); 
r = pcon./repmat(Ph_country2030_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2030_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,2) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,2) = p;
r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,2).*r;
Inv_IX70_2=Inv_IX70(:,2).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,2).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,2) = [0;phhall_all2_IX70_2];
Inv_PVwind2070(:,2) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,2) = [0;cumsum(Inv_noe_IX70_2)];


%% 2025
[m,n]=find(unitmin_IX70>1);
phhall_all2_IX70(m) = 0;
B_utilize_trans_storage_IX70(m,:) = 0;
unitmin_IX70(m) = 0;
EP5_IX70(m,:) = 0;
EF5_IX70(m,:) = 0;
Inv_IX70(m,:) = 0;
Inv_noe_IX70(m,:) = 0;
powerunit_country_IX70(m) = 0;
CO2_year_utilize_trans_storage_IX70(m,:) = 0;
etrans_t_IX70(m) = 0;
etrans_cou1_num_IX70(m,:,:) = 0;

pgen = cumsum(sum(etrans_cou1_num_IX70,3));
pcon = cumsum(reshape(sum(etrans_cou1_num_IX70,2),[size(B_utilize_trans_storage_IX70,1),192]));

load('H:\global-PV-wind\Data\Ph_country2025_IIASA_C1toC8_median.mat');
r = pcon./repmat(Ph_country2025_IIASA_2deg',size(pcon,1),1);
r(find(isnan(r)==1))=0;
r(r>1) =1;
pcon = r.*repmat(Ph_country2025_IIASA_2deg',size(pcon,1),1);
p = [pcon(1,:);diff(pcon)];
Ph_PVwind5_ef_cou(:,1) = sum(p);
Ph_PVwind5_ef_cou_p(:,:,1) = p;
r = sum(p,2)./sum(sum(etrans_cou1_num_IX70,3),2);
r(find(isnan(r)==1))=0;
r(r>1) =1;
CO2_year_utilize_trans_storage_IX70_2=CO2_year_utilize_trans_storage_IX70(:,1).*r;
Inv_IX70_2=Inv_IX70(:,1).*r;
Inv_noe_IX70_2 = Inv_noe_IX70(:,1).*r;
phhall_all2_IX70_2 = phhall_all2_IX70.*r;
Ph_PVwind5(:,1) = [0;phhall_all2_IX70_2];
Inv_PVwind2070(:,1) = [0;cumsum(Inv_IX70_2)];
Inv_PVwind2070_noe(:,1) = [0;cumsum(Inv_noe_IX70_2)];


%%
Inv_CCSPVWind = (Inv_CCS2070_c+Inv_PVwind2070)/10^6;
Inv_clean = (Inv_CCS2070_c+Inv_PVwind2070+Inv_others)/10^6;
save('H:\global-PV-wind\ANS\Inv_PVwind2070_noe_8_2s_6_2a.mat','Inv_PVwind2070_noe','-v7.3')  % million $
save('H:\global-PV-wind\ANS\Inv_PVwind2070_8_2s_6_2a.mat','Inv_PVwind2070','-v7.3')  % million $
save('H:\global-PV-wind\ANS\Inv_CCSPVWind2070_8_2s_6_2a.mat','Inv_CCSPVWind','-v7.3')  % trillion $
save('H:\global-PV-wind\ANS\Inv_clean2070_8_2s_6_2a.mat','Inv_clean','-v7.3')  % trillion $
save('H:\global-PV-wind\ANS\Ph_PVwind2070_8_2s_6_2a.mat','Ph_PVwind5','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\unitmin_IX70_2070_8_2s_6_2a.mat','unitmin_IX70','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\Ph_PVwind5_ef_cou_2070_8_2s_6_2a.mat','Ph_PVwind5_ef_cou','-v7.3')  % TWh/yr
save('H:\global-PV-wind\ANS\Ph_PVwind5_ef_cou_p_2070_8_2s_6_2a.mat','Ph_PVwind5_ef_cou_p','-v7.3')  % TWh/yr

figure
for i = 1:5
    plot(Inv_PVwind2070(:,i)/10^6,Inv_clean(:,i))
    hold on
end

figure
plot(sum(Inv_PVwind2070,2)/5/10^6,sum(Inv_clean,2)/5)
