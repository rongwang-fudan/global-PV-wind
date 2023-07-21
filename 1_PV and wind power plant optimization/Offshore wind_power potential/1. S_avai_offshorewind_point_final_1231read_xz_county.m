%%
tic
clear;
load('H:\Global PV and wind\Data\S_offshorewind_1227_xzz_county.mat');
load('H:\Global PV and wind\Data\SRTM_World2_xz.mat'); % m
[m,n]=find(SRTM30_sea>-1);
S_offshorewind(sub2ind(size(S_offshorewind), m, n))= 0;
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat','S_offshorewind','-v7.3');


%%
tic
clear;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat');
D = 164; %m
Pwp = 8*10^6;%W
CP_unit = Pwp/(7*D*7*D); %W/m2

Rearth =6371; % km
for i = 180*120:-1:(0*120+1)
    gridarea120(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
grid_area = gridarea120 *ones(1,360*120)*10^6; % m2
S_offshorewind=grid_area.*S_offshorewind;% m2

CP = CP_unit.*S_offshorewind/1000; % kW
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\CP1227_xzz_county.mat','CP','-v7.3'); % kW

%%
tic
clear;

UTI_coef=0.95;
ARR_coef=0.90;
Other_coef=0.98;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat');
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\CP1227_xzz_county.mat'); % kW
load('H:\Global PV and wind\Data\dist_s120_1227_2_xz_county.mat'); % km
load('H:\Global PV and wind\Data\SRTM_World2_xz.mat'); % m
clear m
clear n
SRTM30_sea2 = zeros(21600,43200);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1);
SRTM30_sea2(sub2ind(size(SRTM30_sea2), m, n))= SRTM30_sea(sub2ind(size(SRTM30_sea), m, n));

Depth = SRTM30_sea2.*(-1);
Ele_loss =(2.07+(0.073*dist_s120)+(-0.0016*dist_s120.^2 )+(0.000017*dist_s120.^3 )+(-0.000000086*dist_s120.^4 )+(0.000000000157*dist_s120.^5 )+0.0015*Depth+(-0.0000047*Depth.^2 )+(0.0000000082*Depth.^3 )+(-0.0000000000041*Depth.^4 ))/100;
Ele_loss(Ele_loss>1)=1;
Ele_loss(Ele_loss<0)=0;
Ele_coef = ones(21600,43200)-Ele_loss;
clear Ele_loss
clear dist_s120
clear Depth
Ele_coef120 = zeros(21600,43200);
[m,n]=find(S_offshorewind~=0 & SRTM30_sea<=-1);
Ele_coef120(sub2ind(size(Ele_coef120), m, n))= Ele_coef(sub2ind(size(Ele_coef), m, n));
clear S_offshorewind
clear Ele_coef
clear m
clear n
clear SRTM30_sea

load('H:\Global PV and wind\Data\CF_offshore_mean1224.mat'); % CF_offshore_mean1224
Ph = CP.*CF.*UTI_coef.*ARR_coef.*Ele_coef120.*Other_coef.*8760; % kwh/year
P_noCF = CP.*UTI_coef.*ARR_coef.*Ele_coef120.*Other_coef; % kwh/h
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\offshorewind_power_all1227_xzz_county.mat','Ph','-v7.3'); % kwh/year
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\P_noCF_offshorewind_county.mat','P_noCF','-v7.3'); % kwh/h
