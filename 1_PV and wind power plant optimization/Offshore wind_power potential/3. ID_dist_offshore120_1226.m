tic
clear;

load('H:\Global PV and wind\Data\SRTM_World2_xz.mat'); % m
load('H:\Global PV and wind\Data\dist_s120_1227_2_xz_county.mat'); % km
load('H:\Global PV and wind\Data\dist_id_pro120_1227_xz_county.mat'); % pro ID
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\cost_offshorewind1227_xz2_xz_county.mat'); % $/W

dist_s120(dist_s120==0)=1;
dist_s120=ceil(dist_s120);


load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat');

ID_dist_offshore120 = zeros(180*120,360*120);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1 );% & dist_s120~=0
ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), m, n))= dist_s120(sub2ind(size(dist_s120), m, n));
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\ID_dist_offshore120_1227_xz2_xz_county.mat','ID_dist_offshore120','-v7.3'); %  km

