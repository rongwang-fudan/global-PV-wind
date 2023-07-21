tic
clear;
load('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat')
load('H:\Global PV and wind\ANS\unitmin2040_8_2sxz.mat'); % unitmin
unitmin = unitmin(Plant_ID_IX,:);

MAC_CCS = 130.85;
% Baselinese
load('H:\Global PV and wind\ANS\B_utilize_trans_storage_nodiffuse_2040_2040baseline.mat')  % B_utilize_trans_storage_nodiffuse
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_2040_2040baseline.mat')  % CO2_year_utilize_trans_storage
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2040_2040baseline.mat')  % r_power_coop_notech_nopowertrans2070
[m,n]=find(unitmin>4);
r_power_coop_notech_nopowertrans2040(m,:) = [];
CO2_2070(:,1) = sum(r_power_coop_notech_nopowertrans2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,1) = B_utilize_trans_storage_nodiffuse;


% Case A
load('H:\Global PV and wind\ANS\B_utilize_trans_storage_nodiffuse_2040_2040CaseA.mat')  % B_utilize_trans_storage_nodiffuse
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_2040_2040CaseA.mat')  % CO2_year_utilize_trans_storage
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2040_2040CaseA.mat')  % r_power_coop_notech_nopowertrans2070
[m,n]=find(unitmin>4);
r_power_coop_notech_nopowertrans2040(m,:) = [];
CO2_2070(:,2) = sum(r_power_coop_notech_nopowertrans2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,2) = B_utilize_trans_storage_nodiffuse;

% Case B
load('H:\Global PV and wind\ANS\B_utilize_trans_storage_nodiffuse_2040_2040CaseB.mat')  % B_utilize_trans_storage_nodiffuse
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_2040_2040CaseB.mat')  % CO2_year_utilize_trans_storage
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2040_2040CaseB.mat')  % r_power_coop_notech_nopowertrans2070
[m,n]=find(unitmin>4);
r_power_coop_notech_nopowertrans2040(m,:) = [];
CO2_2070(:,3) = sum(r_power_coop_notech_nopowertrans2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,3) = B_utilize_trans_storage_nodiffuse;
clear r_power_coop_notech_nopowertrans2070
clear CO2_year_utilize_trans_storage
clear B_utilize_trans_storage_nodiffuse

%
load('H:\Global PV and wind\ANS\B_utilize_trans_storage_nodiffuse_2040_2040.mat')  % B_utilize_trans_storage_nodiffuse
load('H:\Global PV and wind\ANS\CO2_year_utilize_trans_storage_2040_2040.mat')  % CO2_year_utilize_trans_storage
% Case C: without national cooperation
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopower_nomineraltrans2040.mat')  % r_power_coop_notech_nopower_nomineraltrans2040
[m,n]=find(unitmin>4);
r_power_coop_notech_nopower_nomineraltrans2040(m,:) = [];
CO2_2070(:,4) = sum(r_power_coop_notech_nopower_nomineraltrans2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,4) = B_utilize_trans_storage_nodiffuse;
clear r_power_coop_notech_nopower_nomineraltrans2040

% Case D: with mineral trade
load('H:\Global PV and wind\ANS\r_power_coop_notech_nopowertrans_2040.mat')  % r_power_coop_notech_nopowertrans2070
[m,n]=find(unitmin>4);
r_power_coop_notech_nopowertrans2040(m,:) = [];
CO2_2070(:,5) = sum(r_power_coop_notech_nopowertrans2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,5) = B_utilize_trans_storage_nodiffuse;
clear r_power_coop_notech_nopowertrans2040

% Case E: with power transporatation
load('H:\Global PV and wind\ANS\r_power_coop_notech_2040.mat')  % r_power_coop_notech2070
[m,n]=find(unitmin>4);
r_power_coop_notech2040(m,:) = [];
CO2_2070(:,6) = sum(r_power_coop_notech2040,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,6) = B_utilize_trans_storage_nodiffuse;
clear r_power_coop_notech2040
clear B_utilize_trans_storage_nodiffuse

% Optimal path
load('H:\Global PV and wind\ANS\B_utilize_trans_storage_2040_2040.mat')  % B_utilize_trans_storage
load('H:\Global PV and wind\ANS\r_power_coop_2040xz.mat')  % r_power_coop2070xz
[m,n]=find(unitmin>4);
r_power_coop2040xz(m,:) = [];
CO2_2070(:,7) = sum(r_power_coop2040xz,2).*CO2_year_utilize_trans_storage;
MAC_2070(:,7) = B_utilize_trans_storage;
clear r_power_coop2040xz


%% Fig. 2A
figure
for i = 1:7
    CO2 = CO2_2070(:,i);
    MAC = MAC_2070(:,i);
    [m,n]=find(CO2>=0);
    CO2(m) = [];
    MAC(m) = [];
%     [m,n]=find(MAC>MAC_CCS);
%     CO2(m) = [];
%     MAC(m) = [];
    [MAC_IX,IX] = sort(MAC);
    plot(-cumsum(CO2(IX))/1000, MAC_IX)
    hold on
end
% axis([0 80 -200 150]);
axis([0 50 -200 130.85]);
% axis([0 50 -200 200]);
% xlabel('Abated emissions in 2070 (Gt CO2 y-1)')
% ylabel('MAC of PV, wind and other technologies (2020$ per tCO2)')
% legend('Baseline','Case A','Case B','Case C','Case D','Case E','Optimal path')
set(gca,'xticklabel',[])
set(gca,'xlabel',[])
set(gca,'ylabel',[])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'title',[])

min(MAC_IX)
max(MAC_IX)
-sum(CO2(IX))/1000

