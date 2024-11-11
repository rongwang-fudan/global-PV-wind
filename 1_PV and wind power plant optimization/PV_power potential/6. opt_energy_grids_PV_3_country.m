%%
tic
clear;
load('H:\global-PV-wind\Data\discount_country.mat')
discount = discount/100;
discount(35)=0.05; % per year

load('H:\global-PV-wind\Data\fossilfuel_emissionfactor.mat')  % kg CO2/kWh
% load('H:\world\code\数据处理\initialcost_ratio_country_0111low.mat')  % 1.PV 2.Onshorewind的成本与China的对比
load('H:\global-PV-wind\Data\Line_trans_land_ratio_country.mat')  % Line_trans_land_ratio_country,line, trans, land成本与中国比值，和建造何种电厂无关
load('H:\global-PV-wind\Data\Fixedcost_PV_country.mat')  % 2020USD/kW
% 1.Fixed cost except inverter;	2.Inverter;

load('H:\global-PV-wind\Data\GADM_country120_xz2.mat')
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_w_pv_county.dat','-mat');
powerunit = powerunit_w;
clear powerunit_w
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\costunits01_2_all2.dat','-mat');
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_w_pv_county.dat','-mat');
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_all_w_pv_county.dat','-mat'); % 电厂id，电厂范围内所有面积
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat'); % lines
unitid = unitid_w;
unitid_all = unitid_all_w;
clear unitid_w
clear unitid_all_w

numpowerunit=size(powerunit,1);
rmb2us=1/6.8967; % RMB to USD2019
% moduleprice21=moduleprice+mountprice+equi2thprice+planprice+installprice+PaCprice+accessnetprice;
% inverterprice1

clear moduleprice
clear mountprice
clear equi2thprice
clear installprice
clear PaCprice
clear planprice
clear accessnetprice

% moduleprice2=moduleprice+mountprice+equi2thprice+planprice+installprice+PaCprice+equi1thprice+accessnetprice;
% ratio_module = moduleprice/moduleprice2;
lifetime_inverter=10; % renewed per 10 years
lifetime_power=25;

discount1yr=zeros(192,1);
for i = 1:192
    for t=1:lifetime_power
        discount1yr(i,1)=discount1yr(i,1)+1/(1+discount(i))^(t-1);
    end
end

degrat1yr=zeros(192,1);
degration = 0.0;
for i = 1:192
    for t=1:lifetime_power
        degrat1yr(i,1)=degrat1yr(i,1)+(1-degration)^(t-1)/(1+discount(i,1))^(t-1);
    end
end

discountinverter=zeros(192,1);
for i = 1:192
    for t=1:floor(lifetime_power/lifetime_inverter)
        discountinverter(i,1)=discountinverter(i,1)+1/(1+discount(i,1))^((t-1)*lifetime_inverter);
    end
end

OMratio_majorline=0.01;
OMratio_substation=0.01;
CO2_C=0.2727;
carbonprice=0; % USD2019/tonCO2

powers=zeros(numpowerunit,30);
CPPP111=zeros(numpowerunit,1);
unitid_lcoe = zeros(21600,10800);
unitid_lcoe_all = zeros(21600,10800);
% figure
% i = 0;
load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
region_ID(region_ID==8)=100;
region_ID(region_ID==9)=100;
region_ID(region_ID==10)=100;
region_ID(region_ID~=100)=0;
for coun = 1:192
    moduleprice2 = Fixedcost_PV_country(coun,1)/1000; % 2020USD/W
    inverterprice = Fixedcost_PV_country(coun,2)/1000;% 2020USD/W
    [mp,np]=find(powerunit(:,5)==coun);
    powerunit1 = powerunit(mp,:);
    [mm1,nn1]=find(GADM_country120==coun);
    unitid2=zeros(21600,10800);
    unitid2(sub2ind(size(unitid2), mm1, nn1))=unitid(sub2ind(size(unitid), mm1, nn1));
    unitid22 = unitid2(min(mm1):max(mm1),min(nn1):max(nn1));
    clear unitid2
    
    unitid_all2=zeros(21600,10800);
    unitid_all2(sub2ind(size(unitid_all2), mm1, nn1))=unitid_all(sub2ind(size(unitid_all), mm1, nn1));
    unitid_all22 = unitid_all2(min(mm1):max(mm1),min(nn1):max(nn1));
    clear unitid_all2
    
    for ii=1:size(powerunit1,1)
        %     i=i+1;
        [m,n] = find(unitid_all22==mp(ii));
        mpp = min(m)-1+min(mm1)-1;
        npp = min(n)-1+min(nn1)-1;
        xnet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
        ynet=zeros(max(m)-min(m)+1,max(n)-min(n)+1);
        for ii11=1:max(m)-min(m)+1
            for jj11=1:max(n)-min(n)+1
                xnet(ii11,jj11)=ii11;
                ynet(ii11,jj11)=jj11;
            end
        end
        
        idx_unitid=find(unitid22(min(m):max(m),min(n):max(n))==mp(ii));
        [m_unitid,n_unitid]=find(unitid22(min(m):max(m),min(n):max(n))==mp(ii));
        m_unitid2=m_unitid(:);
        n_unitid2=n_unitid(:);
        
        idx_unitid_all=find(unitid_all22(min(m):max(m),min(n):max(n))==mp(ii));
        [m_unitid_all,n_unitid_all]=find(unitid_all22(min(m):max(m),min(n):max(n))==mp(ii));
        m_unitid2_all=m_unitid_all(:);
        n_unitid2_all=n_unitid_all(:);
        
        centerx=powerunit(mp(ii),1)-min(mm1)+1-min(m)+1; % lat 1/120
        centery=powerunit(mp(ii),2)-min(nn1)+1-min(n)+1; % lon 1/30
        z=floor(abs(xnet(idx_unitid)-centerx)+abs(ynet(idx_unitid)-centery)*4)+1;
        z_all=floor(abs(xnet(idx_unitid_all)-centerx)+abs(ynet(idx_unitid_all)-centery)*4)+1;
        idx=find(costs(:,1)==mp(ii));
        %     idx=find(costs1(:,1)==i);
        lcoeunit=zeros(size(idx,1),20);
        lcoeunit1=zeros(size(idx,1),1);
        CPPP=zeros(size(idx,1),1);
        jopt=0;
        lcoemin=1000;
        for j=1:size(idx,1)
            electricity=costs(idx(j),3); % electricity TWh / year
            capacity=costs(idx(j),11); % capacity potential MW
            if j>1
                %             SRR(j,1)=SRR(j,1)+SR2(idx(j),1);
                %             SRR(j,2)=SRR(j,2)+SR2(idx(j),2);
                CPPP(j,1)=CPPP(j-1,1)+capacity; % capacity potential MW
                lcoeunit(j,1)=lcoeunit(j-1,1)+electricity; % electricity TWh / year
                lcoeunit(j,2)=lcoeunit(j-1,2)+costs(idx(j),9)*rmb2us*Line_trans_land_ratio_country(coun,1); % cost of connection to national grid
                lcoeunit(j,3)=lcoeunit(j-1,3)+costs(idx(j),6)*rmb2us*Line_trans_land_ratio_country(coun,1);% cost of substation
                lcoeunit(j,4)=lcoeunit(j-1,4)+costs(idx(j),7)*rmb2us*Line_trans_land_ratio_country(coun,1); % cost of expanding power unit
                %             lcoeunit(j,4)=lcoeunit(j-1,4)+capacity*lineprice; % cost of expanding power unit
                lcoeunit(j,5)=lcoeunit(j-1,5)+capacity*moduleprice2;   % cost of PV module
                lcoeunit(j,6)=lcoeunit(j-1,6)+capacity*inverterprice; % cost of PV inverter
                lcoeunit(j,7)=lcoeunit(j-1,7)+costs(idx(j),12)*rmb2us*Line_trans_land_ratio_country(coun,1);   % cost of land
                lcoeunit(j,8)=lcoeunit(j-1,8)-electricity*fossilfuel_emissionfactor(coun,1)*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
                lcoeunit(j,9)=lcoeunit(j-1,9)-costs(idx(j),5)*1e-6;  % annual land carbon sink Mton C / year
                lcoeunit(j,10)=lcoeunit(j-1,10)+costs(idx(j),4)*1e-6;   % total land use change emission Mton C
            else
                CPPP(j,1)=capacity; % capacity potential MW
                %             SRR(j,1)=SR2(idx(j),1);
                %             SRR(j,2)=SR2(idx(j),2);
                lcoeunit(j,1)=electricity; % electricity TWh / year
                lcoeunit(j,2)=costs(idx(j),9)*rmb2us*Line_trans_land_ratio_country(coun,1); % cost of connection to national grid
                lcoeunit(j,3)=costs(idx(j),6)*rmb2us*Line_trans_land_ratio_country(coun,1); % cost of substation
                lcoeunit(j,4)=costs(idx(j),7)*rmb2us*Line_trans_land_ratio_country(coun,1); % cost of expanding power unit
                %             lcoeunit(j,4)=capacity*lineprice; % cost of expanding power unit
                lcoeunit(j,5)=capacity*moduleprice2;   % cost of PV module
                lcoeunit(j,6)=capacity*inverterprice; % cost of PV inverter
                lcoeunit(j,7)=costs(idx(j),12)*rmb2us*Line_trans_land_ratio_country(coun,1);   % cost of land
                lcoeunit(j,8)=-electricity*fossilfuel_emissionfactor(coun,1)*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
                lcoeunit(j,9)=-costs(idx(j),5)*1e-6;  % emission from reduced annual land carbon sink Mton C / year
                lcoeunit(j,10)=costs(idx(j),4)*1e-6;   % total land use change emission Mton C
            end
            lcoeunit(j,11)=lcoeunit(j,2)*(1+OMratio_majorline*discount1yr(coun,1)); % cost of connection to national grid
            lcoeunit(j,12)=lcoeunit(j,3)*(1+OMratio_substation*discount1yr(coun,1)); % cost of substation
            lcoeunit(j,13)=lcoeunit(j,4)*(1+OMratio_substation*discount1yr(coun,1)); % cost of distance in power unit
            lcoeunit(j,14)=lcoeunit(j,5)*(1+OMratio_substation*discount1yr(coun,1)); % cost of PV module
            lcoeunit(j,15)=lcoeunit(j,6)*(discountinverter(coun,1)+OMratio_substation*discount1yr(coun,1)); % cost of PV inverter
            lcoeunit(j,16)=lcoeunit(j,7)*(1+OMratio_substation*discount1yr(coun,1)); % cost of land
            %         lcoeunit(j,17)=lcoeunit(j,8)/CO2_C*carbonprice*discount1yr; % saved tax for abated annual CO2 emission from fossil fuel
            %         lcoeunit(j,18)=lcoeunit(j,9)/CO2_C*carbonprice*discount1yr; % tax for emission from reduced annual land carbon sink
            %         lcoeunit(j,19)=lcoeunit(j,10)/CO2_C*carbonprice; % tax for land use change emission
            % LCoE
            lcoeunit(j,20)=sum(lcoeunit(j,11:19),2)/lcoeunit(j,1)/degrat1yr(coun,1)/1000; % LCoE million USD2019/TWh->USD2019/kWh
            %         lcoeunit1(j,1)=sum(lcoeunit(j,12:19),2)/(lcoeunit(j,1))/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh
            if lcoeunit(j,20)<lcoemin
                jopt=j;
                lcoemin=lcoeunit(j,20);
            end
            %         if lcoeunit1(j,1)<lcoemin1 & lcoeunit1(j,1)>0
            %             jopt1=j;
            %             lcoemin1=lcoeunit1(j,1);
            %         end
        end
        % plot LCoE cost
        % list the cost at the lowest LCoE
        if jopt~=0
            %             if region_ID(coun,1)==100
                if sum(costs(idx(1):idx(jopt),11),1)<=100000 % power capacity MW
                    j=jopt;
                else if sum(costs(idx(1):idx(jopt),11),1)>100000
                        for iiiiii = 1:jopt
                            cp_1230(iiiiii,1) = sum(costs(idx(1):idx(iiiiii),11),1);
                        end
%                         [~,Index123] = min(abs(cp_1230-100000));
%                         j=Index123;
                        [Index123,~] = find(cp_1230<=100000);
                        [mm,nn]=find(lcoeunit(Index123,20)==min(min(lcoeunit(Index123,20))));
                        j=mm;
                        clear cp_1230
                    end
                end
            [mmm,nnn] = find(z<=costs(idx(j),8));
            [mmm_all,nnn_all] = find(z_all<=costs(idx(j),8));
            if ~isempty (mmm)
                m_unitid_plant=m_unitid2(sub2ind(size(m_unitid2), mmm, nnn));
                n_unitid_plant=n_unitid2(sub2ind(size(n_unitid2), mmm, nnn));
                unitid_lcoe(sub2ind(size(unitid_lcoe), m_unitid_plant+mpp, n_unitid_plant+npp))=mp(ii);
                m_unitid_plant_all=m_unitid2_all(sub2ind(size(m_unitid2_all), mmm_all, nnn_all));
                n_unitid_plant_all=n_unitid2_all(sub2ind(size(n_unitid2_all), mmm_all, nnn_all));
                unitid_lcoe_all(sub2ind(size(unitid_lcoe_all), m_unitid_plant_all+mpp, n_unitid_plant_all+npp))=mp(ii);
            end
            
            
            powers(mp(ii),1:20)=lcoeunit(j,1:20); % 20 for LCoE USD2019/kWh
            powers(mp(ii),21:29)=lcoeunit(j,11:19)/lcoeunit(j,1)/degrat1yr(coun,1)/1000; % breakdown of LCoE USD2019/kWh
            powers(mp(ii),30)=sum(costs(idx(1):idx(j),11),1); % power capacity MW
            CPPP111(mp(ii),1)=CPPP(j,1); % 20 for LCoE USD2019/kWh
        end
        
        %    ii
    end
    coun
end

save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powers_PV_100GW_3_2_all2_5%_inilow.mat','powers','-v7.3'); % USD2019/kWh
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\CPPP111_PV_100GW_3_2_all2_5%_inilow.mat','CPPP111','-v7.3'); % GW
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_lcoe_100GW_3_2_all2_5%_inilow.dat','unitid_lcoe');
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\unitid_lcoe_all_100GW_3_2_all2_5%_inilow.dat','unitid_lcoe_all');

%%
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_w_pv_county.dat','-mat');
powerunit = powerunit_w;
clear powerunit_w
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','-mat'); % lines
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powers_PV_100GW_3_2_all2_5%_inilow.mat')
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\CPPP111_PV_100GW_3_2_all2_5%_inilow.mat')
numpowerunit=size(powers,1);
[B,IX]=sort(powers(:,20),1);
optpowerunit=zeros(numpowerunit,30);
CPPP111_IX=zeros(numpowerunit,1);
powerunit_IX=zeros(numpowerunit,1);
lines_IX=zeros(size(lines,1),size(lines,2));
powerunit_num_IX=zeros(size(powers,1),6);
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
powerunit_IX=IX;
CPPP111_IX = CPPP111(IX,:);
powerunit_num_IX = powerunit(IX,:);
lines_IX = lines(numlines+IX,:);
optpowerunit(:,1:30) = powers(IX,1:30); % lat lon
optpowerunit(:,31) = IX; % lat lon
optpowerunit(:,32)=optpowerunit(:,8)+optpowerunit(:,9);  % net CO2 emissions Mt C / year = (Abatement-Landsink)
optpowerunit(:,33)=cumsum(optpowerunit(:,32));
optpowerunit(:,34)=cumsum(optpowerunit(:,1));
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\CPPP111_PV_IX_100GW_3_2_all2_5%_inilow.mat','CPPP111_IX','-v7.3'); % MW
optpowerunit_PV = optpowerunit(:,1:34);
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat','optpowerunit_PV'); % power ID
powerunit_IX_PV = powerunit_IX;
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat','powerunit_IX_PV'); % power ID
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_IX_100GW_3_2_all2_5%_inilow.mat','lines_IX'); %
powerunit_num_IX_PV = powerunit_num_IX;
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat','powerunit_num_IX_PV'); %

[m,n] = find(powerunit_num_IX_PV(:,5)==35);
sum(optpowerunit_PV(m,30))/10^6 % TW
sum(optpowerunit_PV(m,1)) % TWh/year

