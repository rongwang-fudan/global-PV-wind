tic
clear;
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
Inv_CCS_dif = diff(Inv_CCS2070_c);
Inv_others_dif = diff(Inv_others);
Inv_others_e_dif = diff(Inv_others_e);
Inv_CCS_dif(m_cho,:) = []; % 增加此电厂导致CCS净成本减少量
Inv_others_dif(m_cho,:) = [];
Inv_others_e_dif(m_cho,:) = [];
Inv_others_c_dif = Inv_others_dif+Inv_others_e_dif; % 增加此电厂导致other renewbales净成本减少量
% Cost_CCS = (Inv_CCS0+sum(Inv_CCS_dif))'/10^6; % trillion $
% Cost_others = (Inv_others0+sum(Inv_others_dif))'/10^6; % trillion $
% Rev_others = (Inv_others_e0+sum(Inv_others_e_dif))'/10^6; % trillion $

%
nn = 1;
for i = 1:192
    for j = 1:192
        if i ~=j
            index_ij_s(nn,1) = i;
            index_ij_s(nn,2) = j;
            nn = nn+1;
        end
    end
end
load('H:\Global PV and wind\ANS\powerunit_country_IX_IX_8_2040_2s_2020s_test6xz.mat'); % powerunit_country_IX_IX
load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat')
powerunit_country_IX_IX = powerunit_country_IX_IX(Plant_ID_IX,:);
Inv_CCS_dif_noconnect_cou = zeros(192,10);
Inv_others_dif_noconnect_cou = zeros(192,10);


%% 2025
load('H:\Global PV and wind\ANS\r_power_coop_2025.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2025.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2025.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2030,2);
[m,n]=find(a==0);
r_power_coop2030(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,1),[1,size(r_power_coop2030,2)]).*r_power_coop2030;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,1),[1,size(r_power_coop2030,2)]).*r_power_coop2030;

asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,1) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,1) = sum(asd(m));    
end

Inv_CCS2025 = (sum(Inv_CCS_dif_connect))';
Inv_others2025 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2025 = Inv_pvwind_mineral_powertrans_tech2030;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2030

load('H:\Global PV and wind\ANS\index_ij_IX_2025.mat')  % index_ij_IX
Inv_pvwind2025xz = Inv_pvwind2025*0;
Inv_others2025xz = Inv_others2025*0;
Inv_CCS2025xz = Inv_CCS2025*0;
Inv_pvwind2025xz(1) = Inv_pvwind2025(1);
Inv_others2025xz(1) = Inv_others2025(1);
Inv_CCS2025xz(1) = Inv_CCS2025(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2025xz(m+1) = Inv_CCS2025(i+1);
    Inv_others2025xz(m+1) = Inv_others2025(i+1);
    Inv_pvwind2025xz(m+1) = Inv_pvwind2025(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2025
clear Inv_others2025
clear Inv_pvwind2025

%% 2030
load('H:\Global PV and wind\ANS\r_power_coop_2030.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2030.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2030.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2030,2);
[m,n]=find(a==0);
r_power_coop2030(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,2),[1,size(r_power_coop2030,2)]).*r_power_coop2030;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,2),[1,size(r_power_coop2030,2)]).*r_power_coop2030;

asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,2) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,2) = sum(asd(m));    
end

Inv_CCS2030 = (sum(Inv_CCS_dif_connect))';
Inv_others2030 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2030 = Inv_pvwind_mineral_powertrans_tech2030;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2030

load('H:\Global PV and wind\ANS\index_ij_IX_2030xz.mat')  % index_ij_IX
% load('H:\Global PV and wind\ANS\index_ij_IX_2040xz.mat')  % index_ij_IX
Inv_pvwind2030xz = Inv_pvwind2030*0;
Inv_others2030xz = Inv_others2030*0;
Inv_CCS2030xz = Inv_CCS2030*0;
Inv_pvwind2030xz(1) = Inv_pvwind2030(1);
Inv_others2030xz(1) = Inv_others2030(1);
Inv_CCS2030xz(1) = Inv_CCS2030(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2030xz(m+1) = Inv_CCS2030(i+1);
    Inv_others2030xz(m+1) = Inv_others2030(i+1);
    Inv_pvwind2030xz(m+1) = Inv_pvwind2030(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2030
clear Inv_others2030
clear Inv_pvwind2030


%% 2030
load('H:\Global PV and wind\ANS\r_power_coop_2035.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2035.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2035.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2030,2);
[m,n]=find(a==0);
r_power_coop2030(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,3),[1,size(r_power_coop2030,2)]).*r_power_coop2030;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,3),[1,size(r_power_coop2030,2)]).*r_power_coop2030;

asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,3) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,3) = sum(asd(m));    
end

Inv_CCS2035 = (sum(Inv_CCS_dif_connect))';
Inv_others2035 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2035 = Inv_pvwind_mineral_powertrans_tech2030;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2030

load('H:\Global PV and wind\ANS\index_ij_IX_2035xz.mat')  % index_ij_IX
Inv_pvwind2035xz = Inv_pvwind2035*0;
Inv_others2035xz = Inv_others2035*0;
Inv_CCS2035xz = Inv_CCS2035*0;
Inv_pvwind2035xz(1) = Inv_pvwind2035(1);
Inv_others2035xz(1) = Inv_others2035(1);
Inv_CCS2035xz(1) = Inv_CCS2035(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2035xz(m+1) = Inv_CCS2035(i+1);
    Inv_others2035xz(m+1) = Inv_others2035(i+1);
    Inv_pvwind2035xz(m+1) = Inv_pvwind2035(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2035
clear Inv_others2035
clear Inv_pvwind2035

%% 2040
load('H:\Global PV and wind\ANS\r_power_coop_2040.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2040.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2040.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,4),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,4),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,4) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,4) = sum(asd(m));    
end

Inv_CCS2040 = (sum(Inv_CCS_dif_connect))';
Inv_others2040 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2040 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2040xz.mat')  % index_ij_IX
Inv_pvwind2040xz = Inv_pvwind2040*0;
Inv_others2040xz = Inv_others2040*0;
Inv_CCS2040xz = Inv_CCS2040*0;
Inv_pvwind2040xz(1) = Inv_pvwind2040(1);
Inv_others2040xz(1) = Inv_others2040(1);
Inv_CCS2040xz(1) = Inv_CCS2040(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2040xz(m+1) = Inv_CCS2040(i+1);
    Inv_others2040xz(m+1) = Inv_others2040(i+1);
    Inv_pvwind2040xz(m+1) = Inv_pvwind2040(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2040
clear Inv_others2040
clear Inv_pvwind2040


%% 2045
load('H:\Global PV and wind\ANS\r_power_coop_2045.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2045.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2045.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,5),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,5),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,5) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,5) = sum(asd(m));    
end

Inv_CCS2045 = (sum(Inv_CCS_dif_connect))';
Inv_others2045 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2045 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2045xz.mat')  % index_ij_IX
Inv_pvwind2045xz = Inv_pvwind2045*0;
Inv_others2045xz = Inv_others2045*0;
Inv_CCS2045xz = Inv_CCS2045*0;
Inv_pvwind2045xz(1) = Inv_pvwind2045(1);
Inv_others2045xz(1) = Inv_others2045(1);
Inv_CCS2045xz(1) = Inv_CCS2045(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2045xz(m+1) = Inv_CCS2045(i+1);
    Inv_others2045xz(m+1) = Inv_others2045(i+1);
    Inv_pvwind2045xz(m+1) = Inv_pvwind2045(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2045
clear Inv_others2045
clear Inv_pvwind2045

%% 2050
load('H:\Global PV and wind\ANS\r_power_coop_2050.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2050.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2050.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,6),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,6),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,6) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,6) = sum(asd(m));    
end

Inv_CCS2050 = (sum(Inv_CCS_dif_connect))';
Inv_others2050 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2050 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2050xz.mat')  % index_ij_IX
Inv_pvwind2050xz = Inv_pvwind2050*0;
Inv_others2050xz = Inv_others2050*0;
Inv_CCS2050xz = Inv_CCS2050*0;
Inv_pvwind2050xz(1) = Inv_pvwind2050(1);
Inv_others2050xz(1) = Inv_others2050(1);
Inv_CCS2050xz(1) = Inv_CCS2050(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2050xz(m+1) = Inv_CCS2050(i+1);
    Inv_others2050xz(m+1) = Inv_others2050(i+1);
    Inv_pvwind2050xz(m+1) = Inv_pvwind2050(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2050
clear Inv_others2050
clear Inv_pvwind2050

%% 2055
load('H:\Global PV and wind\ANS\r_power_coop_2055.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2055.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2055.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,7),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,7),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,7) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,7) = sum(asd(m));    
end

Inv_CCS2055 = (sum(Inv_CCS_dif_connect))';
Inv_others2055 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2055 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2055xz.mat')  % index_ij_IX
Inv_pvwind2055xz = Inv_pvwind2055*0;
Inv_others2055xz = Inv_others2055*0;
Inv_CCS2055xz = Inv_CCS2055*0;
Inv_pvwind2055xz(1) = Inv_pvwind2055(1);
Inv_others2055xz(1) = Inv_others2055(1);
Inv_CCS2055xz(1) = Inv_CCS2055(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2055xz(m+1) = Inv_CCS2055(i+1);
    Inv_others2055xz(m+1) = Inv_others2055(i+1);
    Inv_pvwind2055xz(m+1) = Inv_pvwind2055(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2055
clear Inv_others2055
clear Inv_pvwind2055

%% 2060
load('H:\Global PV and wind\ANS\r_power_coop_2060.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2060.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2060.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,8),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,8),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,8) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,8) = sum(asd(m));    
end

Inv_CCS2060 = (sum(Inv_CCS_dif_connect))';
Inv_others2060 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2060 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2060xz.mat')  % index_ij_IX
Inv_pvwind2060xz = Inv_pvwind2060*0;
Inv_others2060xz = Inv_others2060*0;
Inv_CCS2060xz = Inv_CCS2060*0;
Inv_pvwind2060xz(1) = Inv_pvwind2060(1);
Inv_others2060xz(1) = Inv_others2060(1);
Inv_CCS2060xz(1) = Inv_CCS2060(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2060xz(m+1) = Inv_CCS2060(i+1);
    Inv_others2060xz(m+1) = Inv_others2060(i+1);
    Inv_pvwind2060xz(m+1) = Inv_pvwind2060(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2060
clear Inv_others2060
clear Inv_pvwind2060

%% 2065
load('H:\Global PV and wind\ANS\r_power_coop_2065.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2065.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2065.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,9),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,9),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,9) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,9) = sum(asd(m));    
end

Inv_CCS2065 = (sum(Inv_CCS_dif_connect))';
Inv_others2065 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2065 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2065xz.mat')  % index_ij_IX
Inv_pvwind2065xz = Inv_pvwind2065*0;
Inv_others2065xz = Inv_others2065*0;
Inv_CCS2065xz = Inv_CCS2065*0;
Inv_pvwind2065xz(1) = Inv_pvwind2065(1);
Inv_others2065xz(1) = Inv_others2065(1);
Inv_CCS2065xz(1) = Inv_CCS2065(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2065xz(m+1) = Inv_CCS2065(i+1);
    Inv_others2065xz(m+1) = Inv_others2065(i+1);
    Inv_pvwind2065xz(m+1) = Inv_pvwind2065(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2065
clear Inv_others2065
clear Inv_pvwind2065

%% 2070
load('H:\Global PV and wind\ANS\r_power_coop_2070.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070.mat')  % Inv_pvwind_mineral_powertrans_tech2030
% Connection数目依次是：0,1,2,...,n
a=sum(r_power_coop2040,2);
[m,n]=find(a==0);
r_power_coop2040(m,1) = 1;

Inv_CCS_dif_connect = repmat(Inv_CCS_dif(:,10),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
Inv_others_c_dif_connect = repmat(Inv_others_c_dif(:,10),[1,size(r_power_coop2040,2)]).*r_power_coop2040;
asd = Inv_CCS_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_CCS_dif_noconnect_cou(i,10) = sum(asd(m));    
end
asd = Inv_others_c_dif_connect(:,1);
for i = 1:192
   [m,n]=find(powerunit_country_IX_IX==i); 
   Inv_others_dif_noconnect_cou(i,10) = sum(asd(m));    
end
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
Inv_CCS_dif_noconnect_cou2 = Inv_CCS_dif_noconnect_cou.*repmat(rr2',[192,1]);
Inv_others_dif_noconnect_cou2 = Inv_others_dif_noconnect_cou.*repmat(rr2',[192,1]);
Inv_CCS_dif_noconnect_cou2 = [zeros(size(Inv_CCS_dif_noconnect_cou2,1),1) Inv_CCS_dif_noconnect_cou2];
Inv_others_dif_noconnect_cou2 = [zeros(size(Inv_others_dif_noconnect_cou2,1),1) Inv_others_dif_noconnect_cou2];
for i = 1:10
    Inv_CCS_dif_noconnect_cou_dis2a(:,i) = sum(Inv_CCS_dif_noconnect_cou2(:,i:i+1),2)/2;
    Inv_others_dif_noconnect_cou_dis2a(:,i) = sum(Inv_others_dif_noconnect_cou2(:,i:i+1),2)/2;
end
save('H:\Global PV and wind\ANS\Inv_CCS_dif_noconnect_cou_dis2a_8_2040.mat','Inv_CCS_dif_noconnect_cou_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_others_dif_noconnect_cou_dis2a_8_2040.mat','Inv_others_dif_noconnect_cou_dis2a','-v7.3')  %

Inv_CCS2070 = (sum(Inv_CCS_dif_connect))';
Inv_others2070 = (sum(Inv_others_c_dif_connect))';
Inv_pvwind2070 = Inv_pvwind_mineral_powertrans_tech2040;
clear Inv_CCS_dif_connect
clear Inv_others_c_dif_connect
clear Inv_pvwind_mineral_powertrans_tech2040

load('H:\Global PV and wind\ANS\index_ij_IX_2070xz.mat')  % index_ij_IX
Inv_pvwind2070xz = Inv_pvwind2070*0;
Inv_others2070xz = Inv_others2070*0;
Inv_CCS2070xz = Inv_CCS2070*0;
Inv_pvwind2070xz(1) = Inv_pvwind2070(1);
Inv_others2070xz(1) = Inv_others2070(1);
Inv_CCS2070xz(1) = Inv_CCS2070(1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    Inv_CCS2070xz(m+1) = Inv_CCS2070(i+1);
    Inv_others2070xz(m+1) = Inv_others2070(i+1);
    Inv_pvwind2070xz(m+1) = Inv_pvwind2070(i+1);
    nn = nn+1;
    i
end
clear Inv_CCS2070
clear Inv_others2070
clear Inv_pvwind2070


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
Inv_CCSall_dis = [Inv_CCS2025xz Inv_CCS2030xz Inv_CCS2035xz Inv_CCS2040xz Inv_CCS2045xz Inv_CCS2050xz Inv_CCS2055xz Inv_CCS2060xz Inv_CCS2065xz Inv_CCS2070xz].*repmat(rr2',[size(Inv_CCS2070xz,1),1]);
Inv_othersall_dis = [Inv_others2025xz Inv_others2030xz Inv_others2035xz Inv_others2040xz Inv_others2045xz Inv_others2050xz Inv_others2055xz Inv_others2060xz Inv_others2065xz Inv_others2070xz].*repmat(rr2',[size(Inv_CCS2070xz,1),1]);
Inv_PWall_dis = [Inv_pvwind2025xz Inv_pvwind2030xz Inv_pvwind2035xz Inv_pvwind2040xz Inv_pvwind2045xz Inv_pvwind2050xz Inv_pvwind2055xz Inv_pvwind2060xz Inv_pvwind2065xz Inv_pvwind2070xz].*repmat(rr2',[size(Inv_CCS2070xz,1),1]);

Inv_CCS0_dis2 = [0 Inv_CCS0.*rr2'];
Inv_others_c0_dis2 = [0 Inv_others_c0.*rr2'];
Inv_CCSall_dis2 = [zeros(size(Inv_CCSall_dis,1),1) Inv_CCSall_dis];
Inv_othersall_dis2 = [zeros(size(Inv_othersall_dis,1),1) Inv_othersall_dis];
Inv_PWall_dis2 = [zeros(size(Inv_PWall_dis,1),1) Inv_PWall_dis];
for i = 1:10
    Inv_CCSall_dis2a(:,i) = sum(Inv_CCSall_dis2(:,i:i+1),2)/2;
    Inv_othersall_dis2a(:,i) = sum(Inv_othersall_dis2(:,i:i+1),2)/2;
    Inv_PWall_dis2a(:,i) = sum(Inv_PWall_dis2(:,i:i+1),2)/2;
    Inv_CCS0_dis2a(:,i) = sum(Inv_CCS0_dis2(i:i+1))/2;
    Inv_others_c0_dis2a(:,i) = sum(Inv_others_c0_dis2(i:i+1))/2;
end
sum(Inv_CCSall_dis2a)+sum(Inv_othersall_dis2a)+sum(Inv_PWall_dis2a)+Inv_CCS0_dis2a+Inv_others_c0_dis2a

save('H:\Global PV and wind\ANS\Inv_CCSall_dis2a_8_2040.mat','Inv_CCSall_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_othersall_dis2a_8_2040.mat','Inv_othersall_dis2a','-v7.3')  % 
save('H:\Global PV and wind\ANS\Inv_PWall_dis2a_8_2040.mat','Inv_PWall_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_CCS0_dis2a_8_2040.mat','Inv_CCS0_dis2a','-v7.3')  %
save('H:\Global PV and wind\ANS\Inv_others_c0_dis2a_8_2040.mat','Inv_others_c0_dis2a','-v7.3')  %
