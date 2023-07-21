%%
tic
clear;
load('H:\Global PV and wind\ANS\UHV_Station_country_all.mat')  
position1(:,1)=UHV_Station_country(:,2); % row
position1(:,2)=UHV_Station_country(:,3); % col
position1(:,3:16)=0;
position1(:,10)=UHV_Station_country(:,6);% pro ID(0-3638)
position1(:,11)=UHV_Station_country(:,4);% country ID
position1(:,17)=UHV_Station_country(:,5);% region ID
position1(:,18)=UHV_Station_country(:,7);% power demand (TWh/year)
position_xy = position1;
save('H:\Global PV and wind\ANS\position_substation_2_CN0811_all.mat','position_xy','-v7.3'); % lat;lon
% 1 row； 2 col； 10 pro ID；11 country ID；17 region ID；19 power demand (TWh/year)
