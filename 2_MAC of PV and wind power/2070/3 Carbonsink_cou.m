tic
clear;
load('H:\global-PV-wind\Data\carbonsink12030_xz.mat');  % g C m-2 yr-1
% Hengmao Wang-2022-Global Terrestrial Ecosystem Carbon Flux Inferred from TanSat XCO2 Retrievals

Rearth    =  6371.3;      % km average radium of the earth
for i = 180*120:-1:(0*120+1)
    gridarea1200(180*120+1-i,1)=abs(Rearth^2*(sin(((i/120-90)+1/120)*pi/180)-sin((i/120-90)*pi/180))*1/120*pi/180); %km2
end
gridarea = gridarea1200 * 4 *ones(1,360*30);% 1/120*1/30 单位：km2
clear gridarea1200

land_sink = carbonsink12030.*gridarea;
clear carbonsink12030
% land carbon sink g C/m2/yr -> ton C/yr

load('H:\global-PV-wind\Data\GADM_country120_xz.mat')
for i = 1:192
    [m,n]=find(GADM_country120==i);
    land_sink_country(i,1) = sum(sum(land_sink(sub2ind(size(land_sink), m, n)))); % ton C/yr
end


load('H:\global-PV-wind\Data\fossilfuel_emissionfactor.mat')  % kg CO2/kWh
% fossilfuel_emissionfactor(35) = 0.783;;
pd_landsink=-land_sink_country./0.2727./fossilfuel_emissionfactor/10^6; % TWh/year
pd_landsink(pd_landsink<0)=0;
save('H:\global-PV-wind\ANS\pd_landsink.mat', 'pd_landsink', '-v7.3')  % TWh/year
