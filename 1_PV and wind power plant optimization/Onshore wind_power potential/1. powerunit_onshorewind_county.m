% Author: Rong Wang
% Date: 2021.5.16

tic
clear;
clear global;
global np0 npfraction

load('H:\global-PV-wind\Data\GADM_county120_xz2.mat')
clear GADM_county120_dev
load('H:\global-PV-wind\Data\GADM_country120_xz2.mat')
Rearth =6371; % km
for i = 180*120:-1:(0*120+1)
    gridarea120(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea=gridarea120 * 4 *ones(1,360*30);% 1/120*1/30 单位：Km2
clear gridarea120

% ENERGY 1
load('H:\global-PV-wind\Data\windpower_100m_12to18_day_2.mat'); % wind_kineticenergy(4800x1950) (GWh/grid/day) for 100-m height
% P_h=S×ρ×CF×UTI_coef×ARR_coef÷1000
% S: Available land area（m2）
% ρ=P_Wp/(7D*5D): power density for the current turbine model, W/m2；P_Wp=2500kW；D=103m: Rotor diameter；
% CF: the capacity factor for the curr ent wind turbine model and hub height in the cell；
% UTI_coef= 0.95: the utilization efficiency, due to technical failures and so on (literature values ranged from 0.94 to 0.98)；
% ARR_coef= 0.90: the array efficiency factor, due to wake effects in the wind turbine arrays. (literature values range from 0.7 to 0.925).

% GEOGRAPHIC DATA 2
load('H:\global-PV-wind\Data\LC_CN_xz.mat'); %
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 closed Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
load('H:\global-PV-wind\Data\Slope_world.mat');
% Ground slope (%)
load('H:\global-PV-wind\Data\CF_mean.mat') %单位 km2
CF_mean = CF;
load('H:\global-PV-wind\Data\CF_mean_max.mat') %单位 km2

% suitable land
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 closed Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
load('H:\global-PV-wind\Data\nature_reserve_wildlife120_xz.mat'); %
load('H:\global-PV-wind\Data\DEM_world_xz.mat'); %

suitableland=zeros(21600,10800);
idx1=find(LC_CN==6 | LC_CN==7 | LC_CN==8 | LC_CN==9 | LC_CN==10 | LC_CN==12 | LC_CN==14 | LC_CN==16);
s1=zeros(21600,10800); s1(idx1)=1;
idx2=find(s1==1 & Slope_world<20 & CF>0.2 & nature_reserve_wildlife120~=100 & DEM_world<3000);
suitableland(idx2)=1;
clear CF
clear s1
clear idx1
clear idx2

clear LC_CN
clear DEM_world
clear CF_mean
clear nature_reserve_wildlife120
clear s1
clear idx1
clear idx2
clear Slope_world
clear gridarea
clear windpower_100m_12to18_day

suitableland_onshorewind = suitableland;
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\suitableland_onshorewind.mat','suitableland_onshorewind','-V7.3');
clear suitableland_onshorewind

GADM_country120_2(:,1:5400) = GADM_country120(:,5401:10800);
GADM_country120_2(:,5401:10800) = GADM_country120(:,1:5400);
suitableland_2(:,1:5400) = suitableland(:,5401:10800);
suitableland_2(:,5401:10800) = suitableland(:,1:5400);
GADM_county120_2(:,1:5400) = GADM_county120(:,5401:10800);
GADM_county120_2(:,5401:10800) = GADM_county120(:,1:5400);

a1 = [1:1:10000]';
xnet11 = repmat(a1,1,10000);
a1 = [1:1:10000];
ynet11 = repmat(a1,10000,1);
clear a1
idx=find(xnet11<10000 & ynet11<10000);
a=xnet11(idx);
b=ynet11(idx);
np0 = npoint( a , b );
clear xnet11


% Spatial Analysis
for coun = 1:34 %max(max(GADM_country120)) %1:43% 45:59 %61:
    county_CN1=zeros(21600,10800);
    [m,n] = find(GADM_country120==coun);
    if max(n)-min(n)<4000
        if (~isempty(m))
            county_CN1(sub2ind(size(county_CN1), m, n))=GADM_county120(sub2ind(size(GADM_county120), m, n));
            county_CN2 = county_CN1(min(m):max(m),min(n):max(n));
            position_county(coun,1)=min(m);
            position_county(coun,2)=max(m);
            position_county(coun,3)=min(n);
            position_county(coun,4)=max(n);
            suitableland11 = suitableland(min(m):max(m),min(n):max(n));
            
            powerunit=zeros(200000,6);
            id=0;
            unitid=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_cy=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_all=zeros(max(m)-min(m)+1,max(n)-min(n)+1);% all the grid id of power unit in space
            npfraction=0.6;
            
            aa=unique(county_CN2);
            if max(max(county_CN2))~=0
                % for cy=1:max(max(county_CN))
                for cyy123=2:(size(aa,1))
                    cy = aa(cyy123,1);
                    % 5 types of topograph considered by R. Wang
                    numgrid=zeros(5,13);
                    numgrid_all=zeros(5,13);
                    [mcy,ncy] = find(county_CN2==cy);
                    a1 = [1:1:max(mcy)-min(mcy)+1]';
                    xnet = repmat(a1,1,max(ncy)-min(ncy)+1);
                    a1 = [1:1:max(ncy)-min(ncy)+1];
                    ynet = repmat(a1,max(mcy)-min(mcy)+1,1);
                    clear a1
                    %     idx=find(xnet<4000 & ynet<4000);
                    %     a=xnet(idx);
                    %     b=ynet(idx);
                    %     np0 = npoint( a , b );
                    
                    county_CN= county_CN2(min(mcy):max(mcy),min(ncy):max(ncy));
                    suitableland1= suitableland11(min(mcy):max(mcy),min(ncy):max(ncy));
                    
                    % TYPE 1
                    idx=find(county_CN==cy & suitableland1==1); %
                    idx2=find(county_CN==cy); %
                    if size(idx,1)<=1
                        nnn(floor(cy)) = 0;
                        continue; % no grids
                    end
                    [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                    numgrid(1,1)=num; % suitable grids number of the PV plant
                    numgrid(1,2)=x1+min(mcy)-1;
                    numgrid(1,3)=y1+min(ncy)-1;
                    numgrid_all(1,1)=size(xy12,1);
                    % Grid Information
                    xy_suitabe=zeros(num,13*2);
                    xy_suitabe(1:num,1)=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                    xy_suitabe(1:num,2)=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                    xy_all=zeros(num,13*2);
                    xy_all(1:size(xy12,1),1)=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                    xy_all(1:size(xy12,1),2)=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    
                    
                    % TYPE 2
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & xnet>x1); %  Right
                            idx2=find(county_CN==cy & xnet>x1); %  Right
                        else
                            idx=find(county_CN==cy & suitableland1==1 & xnet<=x1); %  Left
                            idx2=find(county_CN==cy & xnet<=x1); %  Left
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(2,1)=numgrid(2,1)+num;
                        numgrid(2,9+dom)=num;
                        numgrid(2,dom*2)=x+min(mcy)-1;
                        numgrid(2,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
                        numgrid_all(2,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+1))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+2))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+1))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+2))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 3
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1); %  Upper
                            idx2=find(county_CN==cy & ynet>y1); %  Upper
                        else
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1); %  Bottom
                            idx2=find(county_CN==cy & ynet<=y1); %  Bottom
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(3,1)=numgrid(3,1)+num;
                        numgrid(3,9+dom)=num;
                        numgrid(3,dom*2)=x+min(mcy)-1;
                        numgrid(3,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
                        numgrid_all(3,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+5))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+6))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+5))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+6))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 4
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(4,1)=numgrid(4,1)+num;
                        numgrid(4,9+dom)=num;
                        numgrid(4,dom*2)=x+min(mcy)-1;
                        numgrid(4,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
                        numgrid_all(4,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+9))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+10))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+9))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+10))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 5
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet<=x1); %  Upleft
                            idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet>x1); %  Upright
                            idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet>x1); %  Downright
                            idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet<=x1); %  Downleft
                            idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(5,1)=numgrid(5,1)+num;
                        numgrid(5,9+dom)=num;
                        numgrid(5,dom*2)=x+min(mcy)-1;
                        numgrid(5,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
                        numgrid_all(5,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+17))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+18))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+17))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+18))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % Find the best topography for power plants
                    zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
                    if idx1(1)==1
                        if numgrid(1,1)>2
                            id=id+1; % 1 unit
                            powerunit(id,1)=numgrid(1,2);
                            powerunit(id,2)=numgrid(1,3);
                            powerunit(id,3)=numgrid(1,1);
                            powerunit(id,4)=1;
                            powerunit(id,5)=coun;
                            powerunit(id,6)=cy;
                            for i=1:numgrid(1,1)
                                unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
                                unitid_cy(xy_suitabe(i,1),xy_suitabe(i,2))=cy;
                            end
                            for i=1:numgrid_all(1,1)
                                unitid_all(xy_all(i,1),xy_all(i,2))=id;
                            end
                        end
                    elseif idx1(1)==2
                        for dom=1:2
                            if numgrid(2,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(2,dom*2);
                                powerunit(id,2)=numgrid(2,dom*2+1);
                                powerunit(id,3)=numgrid(2,9+dom);
                                powerunit(id,4)=2;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(2,9+dom)
                                    unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=cy;
                                end
                                for i=1:numgrid_all(2,9+dom)
                                    unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                                end
                            end
                        end
                    elseif idx1(1)==3
                        for dom=1:2
                            if numgrid(3,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(3,dom*2);
                                powerunit(id,2)=numgrid(3,dom*2+1);
                                powerunit(id,3)=numgrid(3,9+dom);
                                powerunit(id,4)=3;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(3,9+dom)
                                    unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=cy;
                                end
                                for i=1:numgrid_all(3,9+dom)
                                    unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                                end
                            end
                        end
                    elseif idx1(1)==4
                        for dom=1:4
                            if numgrid(4,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(4,dom*2);
                                powerunit(id,2)=numgrid(4,dom*2+1);
                                powerunit(id,3)=numgrid(4,9+dom);
                                powerunit(id,4)=4;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(4,9+dom)
                                    unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=cy;
                                end
                                for i=1:numgrid_all(4,9+dom)
                                    unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                                end
                            end
                        end
                    elseif idx1(1)==5
                        for dom=1:4
                            if numgrid(5,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(5,dom*2);
                                powerunit(id,2)=numgrid(5,dom*2+1);
                                powerunit(id,3)=numgrid(5,9+dom);
                                powerunit(id,4)=5;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(5,9+dom)
                                    unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=cy;
                                end
                                for i=1:numgrid_all(5,9+dom)
                                    unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                                end
                            end
                        end
                    end
                    %     cy
                end
            end
            t=powerunit;
            powerunit=t(1:id,:); % delete unused data
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat');
            save(filename, 'powerunit','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_all','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_cy','-v7.3');
        end
    end
    if max(n)-min(n)>4000
        [m,n] = find(GADM_country120_2==coun);
        if (~isempty(m))
            county_CN1(sub2ind(size(county_CN1), m, n))=GADM_county120_2(sub2ind(size(GADM_county120_2), m, n));
            county_CN2 = county_CN1(min(m):max(m),min(n):max(n));
            position_county(coun,1)=min(m);
            position_county(coun,2)=max(m);
            position_county(coun,3)=min(n);
            position_county(coun,4)=max(n);
            position_county(coun,5)=2; %是左右位置调换之后的结果
            suitableland11 = suitableland_2(min(m):max(m),min(n):max(n));
            
            powerunit=zeros(200000,6);
            id=0;
            unitid=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_cy=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_all=zeros(max(m)-min(m)+1,max(n)-min(n)+1);% all the grid id of power unit in space
            npfraction=0.6;
            
            aa=unique(county_CN2);
            if max(max(county_CN2))~=0
                % for cy=1:max(max(county_CN))
                for cyy123=2:(size(aa,1))
                    cy = aa(cyy123,1);
                    % 5 types of topograph considered by R. Wang
                    numgrid=zeros(5,13);
                    numgrid_all=zeros(5,13);
                    [mcy,ncy] = find(county_CN2==cy);
                    a1 = [1:1:max(mcy)-min(mcy)+1]';
                    xnet = repmat(a1,1,max(ncy)-min(ncy)+1);
                    a1 = [1:1:max(ncy)-min(ncy)+1];
                    ynet = repmat(a1,max(mcy)-min(mcy)+1,1);
                    clear a1
                    %     idx=find(xnet<4000 & ynet<4000);
                    %     a=xnet(idx);
                    %     b=ynet(idx);
                    %     np0 = npoint( a , b );
                    
                    county_CN= county_CN2(min(mcy):max(mcy),min(ncy):max(ncy));
                    suitableland1= suitableland11(min(mcy):max(mcy),min(ncy):max(ncy));
                    % TYPE 1
                    idx=find(county_CN==cy & suitableland1==1); %
                    idx2=find(county_CN==cy); %
                    if size(idx,1)<=1
                        nnn(floor(cy)) = 0;
                        continue; % no grids
                    end
                    [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                    numgrid(1,1)=num; % suitable grids number of the PV plant
                    numgrid(1,2)=x1+min(mcy)-1;
                    numgrid(1,3)=y1+min(ncy)-1;
                    numgrid_all(1,1)=size(xy12,1);
                    % Grid Information
                    xy_suitabe=zeros(num,13*2);
                    xy_suitabe(1:num,1)=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                    xy_suitabe(1:num,2)=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                    xy_all=zeros(num,13*2);
                    xy_all(1:size(xy12,1),1)=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                    xy_all(1:size(xy12,1),2)=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    
                    % TYPE 2
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & xnet>x1); %  Right
                            idx2=find(county_CN==cy & xnet>x1); %  Right
                        else
                            idx=find(county_CN==cy & suitableland1==1 & xnet<=x1); %  Left
                            idx2=find(county_CN==cy & xnet<=x1); %  Left
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(2,1)=numgrid(2,1)+num;
                        numgrid(2,9+dom)=num;
                        numgrid(2,dom*2)=x+min(mcy)-1;
                        numgrid(2,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
                        numgrid_all(2,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+1))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+2))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+1))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+2))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 3
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1); %  Upper
                            idx2=find(county_CN==cy & ynet>y1); %  Upper
                        else
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1); %  Bottom
                            idx2=find(county_CN==cy & ynet<=y1); %  Bottom
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(3,1)=numgrid(3,1)+num;
                        numgrid(3,9+dom)=num;
                        numgrid(3,dom*2)=x+min(mcy)-1;
                        numgrid(3,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
                        numgrid_all(3,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+5))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+6))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+5))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+6))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 4
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(4,1)=numgrid(4,1)+num;
                        numgrid(4,9+dom)=num;
                        numgrid(4,dom*2)=x+min(mcy)-1;
                        numgrid(4,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
                        numgrid_all(4,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+9))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+10))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+9))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+10))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 5
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet<=x1); %  Upleft
                            idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet>x1); %  Upright
                            idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet>x1); %  Downright
                            idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet<=x1); %  Downleft
                            idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(5,1)=numgrid(5,1)+num;
                        numgrid(5,9+dom)=num;
                        numgrid(5,dom*2)=x+min(mcy)-1;
                        numgrid(5,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
                        numgrid_all(5,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+17))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+18))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+17))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+18))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % Find the best topography for power plants
                    zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
                    if idx1(1)==1
                        if numgrid(1,1)>2
                            id=id+1; % 1 unit
                            powerunit(id,1)=numgrid(1,2);
                            powerunit(id,2)=numgrid(1,3);
                            powerunit(id,3)=numgrid(1,1);
                            powerunit(id,4)=1;
                            powerunit(id,5)=coun;
                            powerunit(id,6)=cy;
                            for i=1:numgrid(1,1)
                                unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
                                unitid_cy(xy_suitabe(i,1),xy_suitabe(i,2))=cy;
                            end
                            for i=1:numgrid_all(1,1)
                                unitid_all(xy_all(i,1),xy_all(i,2))=id;
                            end
                        end
                    elseif idx1(1)==2
                        for dom=1:2
                            if numgrid(2,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(2,dom*2);
                                powerunit(id,2)=numgrid(2,dom*2+1);
                                powerunit(id,3)=numgrid(2,9+dom);
                                powerunit(id,4)=2;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(2,9+dom)
                                    unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=cy;
                                end
                                for i=1:numgrid_all(2,9+dom)
                                    unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                                end
                            end
                        end
                    elseif idx1(1)==3
                        for dom=1:2
                            if numgrid(3,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(3,dom*2);
                                powerunit(id,2)=numgrid(3,dom*2+1);
                                powerunit(id,3)=numgrid(3,9+dom);
                                powerunit(id,4)=3;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(3,9+dom)
                                    unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=cy;
                                end
                                for i=1:numgrid_all(3,9+dom)
                                    unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                                end
                            end
                        end
                    elseif idx1(1)==4
                        for dom=1:4
                            if numgrid(4,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(4,dom*2);
                                powerunit(id,2)=numgrid(4,dom*2+1);
                                powerunit(id,3)=numgrid(4,9+dom);
                                powerunit(id,4)=4;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(4,9+dom)
                                    unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=cy;
                                end
                                for i=1:numgrid_all(4,9+dom)
                                    unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                                end
                            end
                        end
                    elseif idx1(1)==5
                        for dom=1:4
                            if numgrid(5,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(5,dom*2);
                                powerunit(id,2)=numgrid(5,dom*2+1);
                                powerunit(id,3)=numgrid(5,9+dom);
                                powerunit(id,4)=5;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(5,9+dom)
                                    unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=cy;
                                end
                                for i=1:numgrid_all(5,9+dom)
                                    unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                                end
                            end
                        end
                    end
                    %     cy
                end
            end
            t=powerunit;
            powerunit=t(1:id,:); % delete unused data
            

            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat');
            save(filename, 'powerunit','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_all','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_cy','-v7.3');
        end
    end
    
    
    coun
    clear powerunit
end

for coun = 35
    position_county(coun,1)  = 35*120+1;
    position_county(coun,2) = 75*120;
    position_county(coun,3) =(180+73)*30+1;
    position_county(coun,4) =(180+138)*30;
    xnet=zeros(4800,1950);
    ynet=zeros(4800,1950);
    for i=1:4800
        for j=1:1950
            xnet(i,j)=i;
            ynet(i,j)=j;
        end
    end
    idx=find(xnet<1800 & ynet<1800);
    a=xnet(idx);
    b=ynet(idx);
    np0  =  npoint( a , b );
    
    
    county_CN1=zeros(21600,10800);
    suitableland11=zeros(21600,10800);
    [m,n] = find(GADM_country120==35);
    county_CN1(sub2ind(size(county_CN1), m, n))=GADM_county120(sub2ind(size(GADM_county120), m, n));
    suitableland11(sub2ind(size(suitableland11), m, n))=suitableland(sub2ind(size(GADM_county120), m, n));
    county_CN = county_CN1(35*120+1:75*120,(180+73)*30+1:(180+138)*30);
    suitableland22 = suitableland11(35*120+1:75*120,(180+73)*30+1:(180+138)*30);
    clear county_CN1
    clear suitableland11
    
    powerunit=zeros(10000,4);
    id=0;
    unitid=zeros(4800,1950); % suitable grid id of power unit in space
    unitid_all=zeros(4800,1950); % all the grid id of power unit in space
    unitid_cy=zeros(4800,1950); % all the grid id of power unit in space
    npfraction=0.6;
    for cy=1:2373
        % 5 types of topograph considered by R. Wang
        numgrid=zeros(5,13);
        numgrid_all=zeros(5,13);
        % TYPE 1
        idx=find(county_CN==cy & suitableland22==1); %
        idx2=find(county_CN==cy); %
        if size(idx,1)<1
            nnn(cy) = 0;
            continue; % no grids
        end
        [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
        numgrid(1,1)=num; % suitable grids number of the PV plant
        numgrid(1,2)=x1;
        numgrid(1,3)=y1;
        numgrid_all(1,1)=size(xy12,1);
        % Grid Information
        xy_suitabe=zeros(num,13*2);
        xy_suitabe(1:num,1:2)=xy; % suitable grids id of the PV plant
        xy_all=zeros(num,13*2);
        xy_all(1:size(xy12,1),1:2)=xy12; % All the grids id of the PV plant
        
        % TYPE 2
        for dom=1:2
            if dom==1
                idx=find(county_CN==cy & suitableland22==1 & xnet>x1); %  Right
                idx2=find(county_CN==cy & xnet>x1); %  Right
            else
                idx=find(county_CN==cy & suitableland22==1 & xnet<=x1); %  Left
                idx2=find(county_CN==cy & xnet<=x1); %  Left
            end
            if size(idx,1)<1
                continue; % no grids
            end
            [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
            numgrid(2,1)=numgrid(2,1)+num;
            numgrid(2,9+dom)=num;
            numgrid(2,dom*2)=x;
            numgrid(2,dom*2+1)=y;
            numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
            numgrid_all(2,9+dom)=size(xy12,1);
            xy_suitabe(1:num,(dom*2+1):(dom*2+2))=xy;
            xy_all(1:size(xy12,1),(dom*2+1):(dom*2+2))=xy12;
        end
        
        % TYPE 3
        for dom=1:2
            if dom==1
                idx=find(county_CN==cy & suitableland22==1 & ynet>y1); %  Upper
                idx2=find(county_CN==cy & ynet>y1); %  Upper
            else
                idx=find(county_CN==cy & suitableland22==1 & ynet<=y1); %  Bottom
                idx2=find(county_CN==cy & ynet<=y1); %  Bottom
            end
            if size(idx,1)<1
                continue; % no grids
            end
            [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
            numgrid(3,1)=numgrid(3,1)+num;
            numgrid(3,9+dom)=num;
            numgrid(3,dom*2)=x;
            numgrid(3,dom*2+1)=y;
            numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
            numgrid_all(3,9+dom)=size(xy12,1);
            xy_suitabe(1:num,(dom*2+5):(dom*2+6))=xy;
            xy_all(1:size(xy12,1),(dom*2+5):(dom*2+6))=xy12;
        end
        
        % TYPE 4
        for dom=1:4
            if dom==1
                idx=find(county_CN==cy & suitableland22==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
            elseif dom==2
                idx=find(county_CN==cy & suitableland22==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
            elseif dom==3
                idx=find(county_CN==cy & suitableland22==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
            elseif dom==4
                idx=find(county_CN==cy & suitableland22==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
            end
            if size(idx,1)<1
                continue; % no grids
            end
            [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
            numgrid(4,1)=numgrid(4,1)+num;
            numgrid(4,9+dom)=num;
            numgrid(4,dom*2)=x;
            numgrid(4,dom*2+1)=y;
            numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
            numgrid_all(4,9+dom)=size(xy12,1);
            xy_suitabe(1:num,(dom*2+9):(dom*2+10))=xy;
            xy_all(1:size(xy12,1),(dom*2+9):(dom*2+10))=xy12;
        end
        
        % TYPE 5
        for dom=1:4
            if dom==1
                idx=find(county_CN==cy & suitableland22==1 & ynet>y1 & xnet<=x1); %  Upleft
                idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
            elseif dom==2
                idx=find(county_CN==cy & suitableland22==1 & ynet>y1 & xnet>x1); %  Upright
                idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
            elseif dom==3
                idx=find(county_CN==cy & suitableland22==1 & ynet<=y1 & xnet>x1); %  Downright
                idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
            elseif dom==4
                idx=find(county_CN==cy & suitableland22==1 & ynet<=y1 & xnet<=x1); %  Downleft
                idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
            end
            if size(idx,1)<1
                continue; % no grids
            end
            [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
            numgrid(5,1)=numgrid(5,1)+num;
            numgrid(5,9+dom)=num;
            numgrid(5,dom*2)=x;
            numgrid(5,dom*2+1)=y;
            numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
            numgrid_all(5,9+dom)=size(xy12,1);
            xy_suitabe(1:num,(dom*2+17):(dom*2+18))=xy;
            xy_all(1:size(xy12,1),(dom*2+17):(dom*2+18))=xy12;
        end
        
        % Find the best topography for power plants
        zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
        if idx1(1)==1
            if numgrid(1,1)>2
                id=id+1; % 1 unit
                powerunit(id,1)=numgrid(1,2);
                powerunit(id,2)=numgrid(1,3);
                powerunit(id,3)=numgrid(1,1);
                powerunit(id,4)=1;
                powerunit(id,5)=coun;
                powerunit(id,6)=cy;
                for i=1:numgrid(1,1)
                    unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
                    unitid_cy(xy_suitabe(i,1),xy_suitabe(i,2))=cy;
                end
                for i=1:numgrid_all(1,1)
                    unitid_all(xy_all(i,1),xy_all(i,2))=id;
                end
            end
        elseif idx1(1)==2
            for dom=1:2
                if numgrid(2,9+dom)>2
                    id=id+1; % 2 unit
                    powerunit(id,1)=numgrid(2,dom*2);
                    powerunit(id,2)=numgrid(2,dom*2+1);
                    powerunit(id,3)=numgrid(2,9+dom);
                    powerunit(id,4)=2;
                    powerunit(id,5)=coun;
                    powerunit(id,6)=cy;
                    for i=1:numgrid(2,9+dom)
                        unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                        unitid_cy(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=cy;
                    end
                    for i=1:numgrid_all(2,9+dom)
                        unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                    end
                end
            end
        elseif idx1(1)==3
            for dom=1:2
                if numgrid(3,9+dom)>2
                    id=id+1; % 2 unit
                    powerunit(id,1)=numgrid(3,dom*2);
                    powerunit(id,2)=numgrid(3,dom*2+1);
                    powerunit(id,3)=numgrid(3,9+dom);
                    powerunit(id,4)=3;
                    powerunit(id,5)=coun;
                    powerunit(id,6)=cy;
                    for i=1:numgrid(3,9+dom)
                        unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                        unitid_cy(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=cy;
                    end
                    for i=1:numgrid_all(3,9+dom)
                        unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                    end
                end
            end
        elseif idx1(1)==4
            for dom=1:4
                if numgrid(4,9+dom)>2
                    id=id+1; % 4 unit
                    powerunit(id,1)=numgrid(4,dom*2);
                    powerunit(id,2)=numgrid(4,dom*2+1);
                    powerunit(id,3)=numgrid(4,9+dom);
                    powerunit(id,4)=4;
                    powerunit(id,5)=coun;
                    powerunit(id,6)=cy;
                    for i=1:numgrid(4,9+dom)
                        unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                        unitid_cy(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=cy;
                    end
                    for i=1:numgrid_all(4,9+dom)
                        unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                    end
                end
            end
        elseif idx1(1)==5
            for dom=1:4
                if numgrid(5,9+dom)>2
                    id=id+1; % 4 unit
                    powerunit(id,1)=numgrid(5,dom*2);
                    powerunit(id,2)=numgrid(5,dom*2+1);
                    powerunit(id,3)=numgrid(5,9+dom);
                    powerunit(id,4)=5;
                    powerunit(id,5)=coun;
                    powerunit(id,6)=cy;
                    for i=1:numgrid(5,9+dom)
                        unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                        unitid_cy(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=cy;
                    end
                    for i=1:numgrid_all(5,9+dom)
                        unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                    end
                end
            end
        end
        cy
    end
    t=powerunit;
    powerunit=t(1:id,:); % delete unused data
    
    filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat');
    save(filename, 'powerunit','-v7.3');
    
    filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat');
    save(filename, 'unitid','-v7.3');
    
    filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat');
    save(filename, 'unitid_all','-v7.3');
    
    filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat');
    save(filename, 'unitid_cy','-v7.3');
end

for coun = 36:1:max(max(GADM_country120)) %1:43% 45:59 %61:
    county_CN1=zeros(21600,10800);
    [m,n] = find(GADM_country120==coun);
    if max(n)-min(n)<4000
        if (~isempty(m))
            county_CN1(sub2ind(size(county_CN1), m, n))=GADM_county120(sub2ind(size(GADM_county120), m, n));
            county_CN2 = county_CN1(min(m):max(m),min(n):max(n));
            position_county(coun,1)=min(m);
            position_county(coun,2)=max(m);
            position_county(coun,3)=min(n);
            position_county(coun,4)=max(n);
            suitableland11 = suitableland(min(m):max(m),min(n):max(n));
            
            powerunit=zeros(200000,6);
            id=0;
            unitid=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_cy=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_all=zeros(max(m)-min(m)+1,max(n)-min(n)+1);% all the grid id of power unit in space
            npfraction=0.6;
            
            aa=unique(county_CN2);
            if max(max(county_CN2))~=0
                % for cy=1:max(max(county_CN))
                for cyy123=2:(size(aa,1))
                    cy = aa(cyy123,1);
                    % 5 types of topograph considered by R. Wang
                    numgrid=zeros(5,13);
                    numgrid_all=zeros(5,13);
                    [mcy,ncy] = find(county_CN2==cy);
                    a1 = [1:1:max(mcy)-min(mcy)+1]';
                    xnet = repmat(a1,1,max(ncy)-min(ncy)+1);
                    a1 = [1:1:max(ncy)-min(ncy)+1];
                    ynet = repmat(a1,max(mcy)-min(mcy)+1,1);
                    clear a1
                    %     idx=find(xnet<4000 & ynet<4000);
                    %     a=xnet(idx);
                    %     b=ynet(idx);
                    %     np0 = npoint( a , b );
                    
                    county_CN= county_CN2(min(mcy):max(mcy),min(ncy):max(ncy));
                    suitableland1= suitableland11(min(mcy):max(mcy),min(ncy):max(ncy));
                    
                    % TYPE 1
                    idx=find(county_CN==cy & suitableland1==1); %
                    idx2=find(county_CN==cy); %
                    if size(idx,1)<=1
                        nnn(floor(cy)) = 0;
                        continue; % no grids
                    end
                    [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                    numgrid(1,1)=num; % suitable grids number of the PV plant
                    numgrid(1,2)=x1+min(mcy)-1;
                    numgrid(1,3)=y1+min(ncy)-1;
                    numgrid_all(1,1)=size(xy12,1);
                    % Grid Information
                    xy_suitabe=zeros(num,13*2);
                    xy_suitabe(1:num,1)=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                    xy_suitabe(1:num,2)=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                    xy_all=zeros(num,13*2);
                    xy_all(1:size(xy12,1),1)=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                    xy_all(1:size(xy12,1),2)=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    
                    
                    % TYPE 2
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & xnet>x1); %  Right
                            idx2=find(county_CN==cy & xnet>x1); %  Right
                        else
                            idx=find(county_CN==cy & suitableland1==1 & xnet<=x1); %  Left
                            idx2=find(county_CN==cy & xnet<=x1); %  Left
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(2,1)=numgrid(2,1)+num;
                        numgrid(2,9+dom)=num;
                        numgrid(2,dom*2)=x+min(mcy)-1;
                        numgrid(2,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
                        numgrid_all(2,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+1))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+2))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+1))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+2))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 3
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1); %  Upper
                            idx2=find(county_CN==cy & ynet>y1); %  Upper
                        else
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1); %  Bottom
                            idx2=find(county_CN==cy & ynet<=y1); %  Bottom
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(3,1)=numgrid(3,1)+num;
                        numgrid(3,9+dom)=num;
                        numgrid(3,dom*2)=x+min(mcy)-1;
                        numgrid(3,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
                        numgrid_all(3,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+5))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+6))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+5))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+6))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 4
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(4,1)=numgrid(4,1)+num;
                        numgrid(4,9+dom)=num;
                        numgrid(4,dom*2)=x+min(mcy)-1;
                        numgrid(4,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
                        numgrid_all(4,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+9))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+10))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+9))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+10))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 5
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet<=x1); %  Upleft
                            idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet>x1); %  Upright
                            idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet>x1); %  Downright
                            idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet<=x1); %  Downleft
                            idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(5,1)=numgrid(5,1)+num;
                        numgrid(5,9+dom)=num;
                        numgrid(5,dom*2)=x+min(mcy)-1;
                        numgrid(5,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
                        numgrid_all(5,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+17))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+18))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+17))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+18))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % Find the best topography for power plants
                    zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
                    if idx1(1)==1
                        if numgrid(1,1)>2
                            id=id+1; % 1 unit
                            powerunit(id,1)=numgrid(1,2);
                            powerunit(id,2)=numgrid(1,3);
                            powerunit(id,3)=numgrid(1,1);
                            powerunit(id,4)=1;
                            powerunit(id,5)=coun;
                            powerunit(id,6)=cy;
                            for i=1:numgrid(1,1)
                                unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
                                unitid_cy(xy_suitabe(i,1),xy_suitabe(i,2))=cy;
                            end
                            for i=1:numgrid_all(1,1)
                                unitid_all(xy_all(i,1),xy_all(i,2))=id;
                            end
                        end
                    elseif idx1(1)==2
                        for dom=1:2
                            if numgrid(2,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(2,dom*2);
                                powerunit(id,2)=numgrid(2,dom*2+1);
                                powerunit(id,3)=numgrid(2,9+dom);
                                powerunit(id,4)=2;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(2,9+dom)
                                    unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=cy;
                                end
                                for i=1:numgrid_all(2,9+dom)
                                    unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                                end
                            end
                        end
                    elseif idx1(1)==3
                        for dom=1:2
                            if numgrid(3,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(3,dom*2);
                                powerunit(id,2)=numgrid(3,dom*2+1);
                                powerunit(id,3)=numgrid(3,9+dom);
                                powerunit(id,4)=3;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(3,9+dom)
                                    unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=cy;
                                end
                                for i=1:numgrid_all(3,9+dom)
                                    unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                                end
                            end
                        end
                    elseif idx1(1)==4
                        for dom=1:4
                            if numgrid(4,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(4,dom*2);
                                powerunit(id,2)=numgrid(4,dom*2+1);
                                powerunit(id,3)=numgrid(4,9+dom);
                                powerunit(id,4)=4;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(4,9+dom)
                                    unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=cy;
                                end
                                for i=1:numgrid_all(4,9+dom)
                                    unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                                end
                            end
                        end
                    elseif idx1(1)==5
                        for dom=1:4
                            if numgrid(5,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(5,dom*2);
                                powerunit(id,2)=numgrid(5,dom*2+1);
                                powerunit(id,3)=numgrid(5,9+dom);
                                powerunit(id,4)=5;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(5,9+dom)
                                    unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=cy;
                                end
                                for i=1:numgrid_all(5,9+dom)
                                    unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                                end
                            end
                        end
                    end
                    %     cy
                end
            end
            t=powerunit;
            powerunit=t(1:id,:); % delete unused data

            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat');
            save(filename, 'powerunit','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_all','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_cy','-v7.3');
        end
    end
    if max(n)-min(n)>4000
        [m,n] = find(GADM_country120_2==coun);
        if (~isempty(m))
            county_CN1(sub2ind(size(county_CN1), m, n))=GADM_county120_2(sub2ind(size(GADM_county120_2), m, n));
            county_CN2 = county_CN1(min(m):max(m),min(n):max(n));
            position_county(coun,1)=min(m);
            position_county(coun,2)=max(m);
            position_county(coun,3)=min(n);
            position_county(coun,4)=max(n);
            position_county(coun,5)=2; %是左右位置调换之后的结果
            suitableland11 = suitableland_2(min(m):max(m),min(n):max(n));
            
            powerunit=zeros(200000,6);
            id=0;
            unitid=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_cy=zeros(max(m)-min(m)+1,max(n)-min(n)+1); % suitable grid id of power unit in space
            unitid_all=zeros(max(m)-min(m)+1,max(n)-min(n)+1);% all the grid id of power unit in space
            npfraction=0.6;
            
            aa=unique(county_CN2);
            if max(max(county_CN2))~=0
                % for cy=1:max(max(county_CN))
                for cyy123=2:(size(aa,1))
                    cy = aa(cyy123,1);
                    % 5 types of topograph considered by R. Wang
                    numgrid=zeros(5,13);
                    numgrid_all=zeros(5,13);
                    [mcy,ncy] = find(county_CN2==cy);
                    a1 = [1:1:max(mcy)-min(mcy)+1]';
                    xnet = repmat(a1,1,max(ncy)-min(ncy)+1);
                    a1 = [1:1:max(ncy)-min(ncy)+1];
                    ynet = repmat(a1,max(mcy)-min(mcy)+1,1);
                    clear a1
                    %     idx=find(xnet<4000 & ynet<4000);
                    %     a=xnet(idx);
                    %     b=ynet(idx);
                    %     np0 = npoint( a , b );
                    
                    county_CN= county_CN2(min(mcy):max(mcy),min(ncy):max(ncy));
                    suitableland1= suitableland11(min(mcy):max(mcy),min(ncy):max(ncy));
                    % TYPE 1
                    idx=find(county_CN==cy & suitableland1==1); %
                    idx2=find(county_CN==cy); %
                    if size(idx,1)<=1
                        nnn(floor(cy)) = 0;
                        continue; % no grids
                    end
                    [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                    numgrid(1,1)=num; % suitable grids number of the PV plant
                    numgrid(1,2)=x1+min(mcy)-1;
                    numgrid(1,3)=y1+min(ncy)-1;
                    numgrid_all(1,1)=size(xy12,1);
                    % Grid Information
                    xy_suitabe=zeros(num,13*2);
                    xy_suitabe(1:num,1)=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                    xy_suitabe(1:num,2)=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                    xy_all=zeros(num,13*2);
                    xy_all(1:size(xy12,1),1)=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                    xy_all(1:size(xy12,1),2)=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    
                    % TYPE 2
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & xnet>x1); %  Right
                            idx2=find(county_CN==cy & xnet>x1); %  Right
                        else
                            idx=find(county_CN==cy & suitableland1==1 & xnet<=x1); %  Left
                            idx2=find(county_CN==cy & xnet<=x1); %  Left
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(2,1)=numgrid(2,1)+num;
                        numgrid(2,9+dom)=num;
                        numgrid(2,dom*2)=x+min(mcy)-1;
                        numgrid(2,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
                        numgrid_all(2,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+1))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+2))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+1))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+2))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 3
                    for dom=1:2
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1); %  Upper
                            idx2=find(county_CN==cy & ynet>y1); %  Upper
                        else
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1); %  Bottom
                            idx2=find(county_CN==cy & ynet<=y1); %  Bottom
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(3,1)=numgrid(3,1)+num;
                        numgrid(3,9+dom)=num;
                        numgrid(3,dom*2)=x+min(mcy)-1;
                        numgrid(3,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
                        numgrid_all(3,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+5))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+6))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+5))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+6))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 4
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(4,1)=numgrid(4,1)+num;
                        numgrid(4,9+dom)=num;
                        numgrid(4,dom*2)=x+min(mcy)-1;
                        numgrid(4,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
                        numgrid_all(4,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+9))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+10))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+9))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+10))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % TYPE 5
                    for dom=1:4
                        if dom==1
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet<=x1); %  Upleft
                            idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
                        elseif dom==2
                            idx=find(county_CN==cy & suitableland1==1 & ynet>y1 & xnet>x1); %  Upright
                            idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
                        elseif dom==3
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet>x1); %  Downright
                            idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
                        elseif dom==4
                            idx=find(county_CN==cy & suitableland1==1 & ynet<=y1 & xnet<=x1); %  Downleft
                            idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
                        end
                        if size(idx,1)<=1
                            continue; % no grids
                        end
                        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
                        numgrid(5,1)=numgrid(5,1)+num;
                        numgrid(5,9+dom)=num;
                        numgrid(5,dom*2)=x+min(mcy)-1;
                        numgrid(5,dom*2+1)=y+min(ncy)-1;
                        numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
                        numgrid_all(5,9+dom)=size(xy12,1);
                        xy_suitabe(1:num,(dom*2+17))=xy(:,1)+min(mcy)-1; % suitable grids id of the PV plant
                        xy_suitabe(1:num,(dom*2+18))=xy(:,2)+min(ncy)-1; % suitable grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+17))=xy12(:,1)+min(mcy)-1;  % All the grids id of the PV plant
                        xy_all(1:size(xy12,1),(dom*2+18))=xy12(:,2)+min(ncy)-1; % All the grids id of the PV plant
                    end
                    
                    % Find the best topography for power plants
                    zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
                    if idx1(1)==1
                        if numgrid(1,1)>2
                            id=id+1; % 1 unit
                            powerunit(id,1)=numgrid(1,2);
                            powerunit(id,2)=numgrid(1,3);
                            powerunit(id,3)=numgrid(1,1);
                            powerunit(id,4)=1;
                            powerunit(id,5)=coun;
                            powerunit(id,6)=cy;
                            for i=1:numgrid(1,1)
                                unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
                                unitid_cy(xy_suitabe(i,1),xy_suitabe(i,2))=cy;
                            end
                            for i=1:numgrid_all(1,1)
                                unitid_all(xy_all(i,1),xy_all(i,2))=id;
                            end
                        end
                    elseif idx1(1)==2
                        for dom=1:2
                            if numgrid(2,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(2,dom*2);
                                powerunit(id,2)=numgrid(2,dom*2+1);
                                powerunit(id,3)=numgrid(2,9+dom);
                                powerunit(id,4)=2;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(2,9+dom)
                                    unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=cy;
                                end
                                for i=1:numgrid_all(2,9+dom)
                                    unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                                end
                            end
                        end
                    elseif idx1(1)==3
                        for dom=1:2
                            if numgrid(3,9+dom)>2
                                id=id+1; % 2 unit
                                powerunit(id,1)=numgrid(3,dom*2);
                                powerunit(id,2)=numgrid(3,dom*2+1);
                                powerunit(id,3)=numgrid(3,9+dom);
                                powerunit(id,4)=3;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(3,9+dom)
                                    unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=cy;
                                end
                                for i=1:numgrid_all(3,9+dom)
                                    unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                                end
                            end
                        end
                    elseif idx1(1)==4
                        for dom=1:4
                            if numgrid(4,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(4,dom*2);
                                powerunit(id,2)=numgrid(4,dom*2+1);
                                powerunit(id,3)=numgrid(4,9+dom);
                                powerunit(id,4)=4;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(4,9+dom)
                                    unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=cy;
                                end
                                for i=1:numgrid_all(4,9+dom)
                                    unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                                end
                            end
                        end
                    elseif idx1(1)==5
                        for dom=1:4
                            if numgrid(5,9+dom)>2
                                id=id+1; % 4 unit
                                powerunit(id,1)=numgrid(5,dom*2);
                                powerunit(id,2)=numgrid(5,dom*2+1);
                                powerunit(id,3)=numgrid(5,9+dom);
                                powerunit(id,4)=5;
                                powerunit(id,5)=coun;
                                powerunit(id,6)=cy;
                                for i=1:numgrid(5,9+dom)
                                    unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                                    unitid_cy(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=cy;
                                end
                                for i=1:numgrid_all(5,9+dom)
                                    unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                                end
                            end
                        end
                    end
                    %     cy
                end
            end
            t=powerunit;
            powerunit=t(1:id,:); % delete unused data
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat');
            save(filename, 'powerunit','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_all','-v7.3');
            
            filename = strcat('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat');
            save(filename, 'unitid_cy','-v7.3');
        end
    end
    
    
    coun
    clear powerunit
end

save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\position_county_onshorewind_500MW_sameC.dat','position_county');


