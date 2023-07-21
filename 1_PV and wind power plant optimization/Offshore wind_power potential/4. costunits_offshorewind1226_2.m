%%
tic
clear;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\cost_offshorewind1227_xz2_xz_county.mat'); % $/W
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\S_offshorewind_1227_xzz_county2.mat'); 

D = 164; %m
Pwp = 8*10^6;%W
CP_unit = Pwp/(7*D*7*D); %W/m2
Rearth =6371; % km
for i = 180*120:-1:(0*120+1)
    gridarea120(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
grid_area = gridarea120 *ones(1,360*120)*10^6; % m2
clear gridarea120

S_offshorewind=grid_area.*S_offshorewind;% m2
clear grid_area
CP = CP_unit.*S_offshorewind/1000; % kW
Cost = cost_offshorewind.*CP*1000/10^6; % million dollar
clear cost_offshorewind

load('H:\Global PV and wind\Data\carbonsink12030_xz.mat');  % g C m-2 yr-1
land_sink = zeros(21600,43200);
for i = 1:21600
    for j = 1:10800
        land_sink(i,j*4-3:j*4)=carbonsink12030(i,j);        
    end
    i
end
clear carbonsink12030

Rearth    =  6371.3;      % km average radium of the earth
for i = 180*120:-1:(0*120+1)
    gridarea1200(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea = gridarea1200 *ones(1,360*120);% 1/120*1/120 单位：km2
clear gridarea1200
land_sink_2=land_sink.*gridarea.*S_offshorewind; % land carbon sink gC/m2/yr -> ton C/yr
clear land_sink
clear gridarea
clear suitablity

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\ID_dist_offshore120_1227_xz2_xz_county.mat'); %
load('H:\Global PV and wind\Data\dist_id_pro120_1227_xz_county.mat'); % ID

 dist_id_pro120_1 = ones(21600,43200)*(-1);
 [m,n]=find(S_offshorewind~=0);
 dist_id_pro120_1(sub2ind(size(dist_id_pro120_1), m, n))= dist_id_pro120(sub2ind(size(dist_id_pro120), m, n));
clear dist_id_pro120
clear m
clear n

load('H:\Global PV and wind\Data\dist_id_country120_1227_xz_county.mat');
bb=unique(dist_id_country120);
bb1 = bb(2:end,1);
clear bb

load('H:\Global PV and wind\Data\Wind_speed_100m_mean120.mat'); % m/s
Wind_speed_100m_mean120(find(isnan(Wind_speed_100m_mean120)==1))=0;


load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\offshorewind_power_all1227_xzz_county.mat'); % ,'Ph','-v7.3' kwh/year

costunits_offshorewind = zeros(10000000,10);
tt=1;
nnp = 1;
for coun = 1:size(bb1,1)
    [mmm1,nnn1]=find(dist_id_country120==bb1(coun,1));
    mm1m = min(min(mmm1));
    nnn1m = min(min(nnn1));
    dist_id_pro120_2 = ones(max(mmm1)-min(mmm1)+1,max(nnn1)-min(nnn1)+1)*(-1);
    dist_id_pro120_2(sub2ind(size(dist_id_pro120_2), mmm1-mm1m+1, nnn1-nnn1m+1)) = dist_id_pro120_1(sub2ind(size(dist_id_pro120_1), mmm1, nnn1));
    nnnnnn=unique(dist_id_pro120_2);    
for i=2:(size(nnnnnn,1))
    proo = nnnnnn(i,1);
    [m,n]=find(dist_id_pro120_2==proo);
    ID_dist_offshore120_pro = zeros(max(max(m))-min(min(m))+1,max(max(n))-min(min(n))+1);
    mm1 = min(min(m));
    nn1 = min(min(n));
    ID_dist_offshore120_pro(sub2ind(size(ID_dist_offshore120_pro), m-mm1+1, n-nn1+1))= ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), m+mm1m-1, n+nnn1m-1));
    for j=1:max(max(ID_dist_offshore120_pro))
        [m,n]=find(ID_dist_offshore120_pro==j);
        if ~isempty (m)
        costunits_offshorewind(tt,1)= nnp; % power plant ID
        costunits_offshorewind(tt,8)= proo; % county ID for other countries, pro ID for China
        costunits_offshorewind(tt,9)= bb1(coun,1); % country ID
        costunits_offshorewind(tt,2)= j; % 距离
        costunits_offshorewind(tt,3)= sum(sum(CP(sub2ind(size(CP), m+mm1+mm1m-2, n+nn1+nnn1m-2)))); % CP,kW
        costunits_offshorewind(tt,4)= sum(sum(Cost(sub2ind(size(Cost), m+mm1+mm1m-2, n+nn1+nnn1m-2))));% Cost,million dollar
        costunits_offshorewind(tt,5)= sum(sum(Ph(sub2ind(size(Ph), m+mm1+mm1m-2, n+nn1+nnn1m-2)))); % kwh/year
        costunits_offshorewind(tt,6)= sum(sum(S_offshorewind(sub2ind(size(S_offshorewind), m+mm1+mm1m-2, n+nn1+nnn1m-2))))/10^6; % km2
        costunits_offshorewind(tt,7)= sum(sum(Wind_speed_100m_mean120(sub2ind(size(Wind_speed_100m_mean120), m+mm1+mm1m-2, n+nn1+nnn1m-2))))/size(m,1); % mean wind speed, m/s
        costunits_offshorewind(tt,10)= sum(sum(land_sink_2(sub2ind(size(land_sink_2), m+mm1+mm1m-2, n+nn1+nnn1m-2))))*0.02; % ton C/yr
        tt=tt+1;
        end
%         ja
    end
%     i
        if max(max(ID_dist_offshore120_pro))~=0
    nnp = nnp+1;
        end
end
coun
end
clear land_sink_2
clear Wind_speed_100m_mean120
clear S_offshorewind
clear Ph
clear Cost
clear CP

asd = sum(costunits_offshorewind,2);
[m,n]=find(asd==0);
costunits_offshorewind(m,:)=[];

save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\costunits_offshorewind1227_xz2_xz_county.dat','costunits_offshorewind');
