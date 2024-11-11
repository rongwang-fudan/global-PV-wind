tic
clear;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat'); 
load('H:\global-PV-wind\Data\SRTM_World2_xz.mat'); % m
SRTM30_sea1 = zeros(21600,43200);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1);
SRTM30_sea1(sub2ind(size(SRTM30_sea1), m, n))= SRTM30_sea(sub2ind(size(SRTM30_sea), m, n));

load('H:\global-PV-wind\Data\dist_s120_1227_2_xz_county.mat'); % km
C0 = 2000/1000; % $/W
C = C0.*(0.0084.*(-SRTM30_sea1)+0.8368).*(0.0057.*dist_s120+0.7714);
cost_offshorewind = zeros(21600,43200);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1); % & dist_s120~=0
cost_offshorewind(sub2ind(size(cost_offshorewind), m, n))= C(sub2ind(size(C), m, n));

save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\cost_offshorewind1227_xz2_xz_county.mat','cost_offshorewind','-v7.3'); % $/W

