tic
clear;
load('H:\Global PV and wind\Data\GADM_country120_xz2.mat')
load('H:\Global PV and wind\Data\GADM_pro120_xz2.mat')
pro_id = unique(GADM_pro120);
pro_id(1)=[];

Rearth    =  6371.3;      % km average radium of the earth
for i = 180*120:-1:(0*120+1)
    gridarea1200(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea = gridarea1200 * 4 *ones(1,360*30);% 1/120*1/30, km2
clear gridarea1200

area_pro = zeros(size(pro_id,1),1);
for i = 1:size(pro_id,1)
    [m,n]=find(GADM_pro120==pro_id(i));
    area_pro(pro_id(i)+1,1) = sum(sum(gridarea(sub2ind(size(gridarea), m, n))));
    i
end

load('H:\Global PV and wind\Data\mineral_production2021.mat');  %  thousand metric tons/year
load('H:\Global PV and wind\Data\mineral_reserve2021.mat');  %  thousand metric tons
% 1 Copper, 2 Zinc, 3 Nickel, 4 Silicon, 
% 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% Molybdenum is assumed to be enough（col 7）
mineral_production_pro = zeros(size(pro_id,1),8);
mineral_reserve_pro = zeros(size(pro_id,1),8);
load('H:\Global PV and wind\Data\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
r_area_pro = zeros(size(pro_id,1),1);
for i = 1:1:192
    [m,n]=find(GADM_country120==i);
    area_cou(i,1) = sum(sum(gridarea(sub2ind(size(gridarea), m, n))));
    [m,n]=find(ID_pro(:,3)==i);
    r_area_pro(ID_pro(m,1)+1) = area_pro(ID_pro(m,1)+1)./ area_cou(i,1);    
    mineral_production_pro(ID_pro(m,1)+1,:) = r_area_pro(ID_pro(m,1)+1,1).*(mineral_production2021(i,:));    
    mineral_reserve_pro(ID_pro(m,1)+1,:) = r_area_pro(ID_pro(m,1)+1,1).*(mineral_reserve2021(i,:));   
    i
end
mineral_production_pro(:,7) = 10^10;
mineral_reserve_pro(:,7) = 10^10;
mineral_reserve_pro(:,4) = 10^10;

save('H:\Global PV and wind\ANS\mineral_production_pro.mat','mineral_production_pro') 
save('H:\Global PV and wind\ANS\mineral_reserve_pro.mat','mineral_reserve_pro') 