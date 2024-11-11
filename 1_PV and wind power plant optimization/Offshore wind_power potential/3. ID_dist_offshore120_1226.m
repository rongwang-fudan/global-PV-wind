tic
clear;

load('H:\global-PV-wind\Data\SRTM_World2_xz.mat'); % m
load('H:\global-PV-wind\Data\dist_s120_1227_2_xz_county.mat'); % km
load('H:\global-PV-wind\Data\dist_id_pro120_1227_xz_county.mat'); % pro ID
% load('H:\world1\data\2dist_id_pro120_1227_2_xz.mat'); % pro ID
% 3:Fujian; 5:Guangdong; 6: Guangxi; 8:Hainan; 9:Hebei; 17:Jiangsu; 
% 18:Liaoning; 25: Shandong; 26:Shanghai; 28: Tianjin; 32: Yunnan; 33:Zhejiang; 
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\cost_offshorewind1227_xz2_xz_county.mat'); % $/W

dist_s120(dist_s120==0)=1;
dist_s120=ceil(dist_s120);


load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat');

ID_dist_offshore120 = zeros(180*120,360*120);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1 );% & dist_s120~=0
ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), m, n))= dist_s120(sub2ind(size(dist_s120), m, n));
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\ID_dist_offshore120_1227_xz2_xz_county.mat','ID_dist_offshore120','-v7.3'); %  km

