tic
clear;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\position_county_onshorewind_500MW_sameC.dat','-mat');
load('H:\Global PV and wind\Data\GADM_country120_xz2.mat') 

powerunit_w = [];
unitid_w = zeros(21600,10800);
unitid_all_w = zeros(21600,10800);
unitid_cy_w = zeros(21600,10800);
num=0;
for iii = 1:max(max(GADM_country120))
    coun = iii;
    load(strcat('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_',num2str(coun),'_sameC.mat'));
    load(strcat('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_',num2str(coun),'_sameC.mat'));
    load(strcat('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_',num2str(coun),'_sameC.mat'));
    load(strcat('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_',num2str(coun),'_sameC.mat'));
    if (~isempty(powerunit))
    unitid_w11 = zeros(21600,10800);
    unitid_all_w11 = zeros(21600,10800);
    unitid_cy_w11 = zeros(21600,10800);
    if position_county(coun,5)~=2
        powerunit_w1(:,1) = powerunit(:,1)+position_county(coun,1)-1;
        powerunit_w1(:,2) = powerunit(:,2)+position_county(coun,3)-1;
        powerunit_w1(:,3:6) = powerunit(:,3:6);
        unitid_w11(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid;
        unitid_all_w11(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid_all;
        unitid_cy_w11(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid_cy;
    end
    
    if position_county(coun,5)==2
        unitid_w111 = zeros(21600,10800);
        unitid_all_w111 = zeros(21600,10800);
        unitid_pro_w111 = zeros(21600,10800);
        powerunit_w1(:,1) = powerunit(:,1)+position_county(coun,1)-1;
        powerunit_w1(:,2) = powerunit(:,2)+position_county(coun,3)-1;
        powerunit_w1(:,3:6) = powerunit(:,3:6);
        for i = 1:size(powerunit_w1,1) 
            if powerunit_w1(i,2)<=5400
                powerunit_w1(i,2) = powerunit_w1(i,2)+5400;
            else if powerunit_w1(i,2)>5400
                powerunit_w1(i,2) = powerunit_w1(i,2)-5400;
                end
            end
        end
        unitid_w111(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid;
        unitid_all_w111(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid_all;
        unitid_pro_w111(position_county(coun,1):position_county(coun,2),position_county(coun,3):position_county(coun,4))=unitid_cy;

        unitid_w11(:,1:5400) = unitid_w111(:,5401:10800);
        unitid_w11(:,5401:10800) = unitid_w111(:,1:5400);
        unitid_all_w11(:,1:5400) = unitid_all_w111(:,5401:10800);
        unitid_all_w11(:,5401:10800) = unitid_all_w111(:,1:5400);
        unitid_cy_w11(:,1:5400) = unitid_pro_w111(:,5401:10800);
        unitid_cy_w11(:,5401:10800) = unitid_pro_w111(:,1:5400);
    end
    
    powerunit_w = [powerunit_w;powerunit_w1];
    clear powerunit_w1
    [m,n]=find(unitid_w11~=0);
    unitid_w(sub2ind(size(unitid_w), m, n))= unitid_w11(sub2ind(size(unitid_w11), m, n))+num;
    [m,n]=find(unitid_all_w11~=0);
    unitid_all_w(sub2ind(size(unitid_all_w), m, n))= unitid_all_w11(sub2ind(size(unitid_all_w11), m, n))+num;
    [m,n]=find(unitid_cy_w11~=0);
    unitid_cy_w(sub2ind(size(unitid_cy_w), m, n))= unitid_cy_w11(sub2ind(size(unitid_cy_w11), m, n));
    num = num+max(max(unitid_w11));
    end
    coun
end
unitid_cy_w = floor(unitid_cy_w);

save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_w_onshorewind_county_500MW.dat','unitid_w');
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_all_w_onshorewind_county_500MW.dat','unitid_all_w');
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\unitid_cy_w_onshorewind_county_500MW.dat','unitid_cy_w');
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_w_onshorewind_county_500MW.dat','powerunit_w');
