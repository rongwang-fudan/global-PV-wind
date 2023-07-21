tic
clear;
r_module_hardwarecost_ons = 0.704960836;
zerotracker_year_country = ones(192,1)*2040;
load('H:\Global PV and wind\Data\Country_classify.mat')
% 1:developed country; 2: developing country; 3: Least Developed Countries
Country_classify(Country_classify>=2)=2;
load('H:\Global PV and wind\Data\CP_PV_ons_off2020.mat'); % MW 2020
% 1 Solar photovoltaic, 2 Onshore wind energy, 3 Offshore wind energy
CP_PV_ons_off2020_2 = CP_PV_ons_off2020;
CP_PV_ons_off2020 (CP_PV_ons_off2020==0)=NaN;
CP_mean = mean (CP_PV_ons_off2020,'omitnan');  % MW 2020
clear CP_PV_ons_off2020
CP_PV_ons_off2020 = CP_PV_ons_off2020_2;
load('H:\Global PV and wind\Data\powerdemand_monhour2070_ele288_1112_2_max.mat')  % 288*34 % TWh/h 电气化后
p_all2060=sum(powerdemand_monhour2070_ele288)'/288*8760; % 	TWh/year
clear powerdemand_monhour2060_ele

load('H:\Global PV and wind\Data\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
for i = 1:1:192
    [m,n]=find(ID_pro(:,3)==i);
    FID_pro = unique(ID_pro(m,1));
    p_all2060_country(i,1) =sum(p_all2060(FID_pro+1,1));
    clear FID_pro
    i
end
clear p_all2060
clear ID_pro

%%
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_IX_100GW_3_2_all2_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % 电厂编号
clear powerunit_IX_PV

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\optpowerunit_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\tranmission_lines_IX_100GW_3_2_all_5%_inilow.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_num_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind;

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\optpowerunit_offshorewind_100GW_county_5%.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_IX_offshorewind_100GW_county_5%.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_IX_100GW_county_5%.mat');  % lines_IX
lines_IX(size(optpowerunit_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind;

% %
% load('H:\Global PV and wind\ANS\ID_pro3.mat') % FID 	FIRST_ID_0  ID_country120_0214	FIRST_ID_1
% for i = 1:size(off_pro_IX,1)
%     a = off_pro_IX(i,1);
%     [m1,n1]=find(ID_pro(:,1)==a);
%     off_pro_IX(i,2) = unique(ID_pro(m1,3));
%     i
% end
% clear ID_pro2

% powerunit_country_IX = [powerunit_num_IX_PV(:,5);powerunit_num_IX_onshorewind(:,5);off_pro_IX(:,2)];
powerunit_num_IX_PV_ori = powerunit_num_IX_PV(:,5);
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind(:,5);
off_pro_IX_ori = off_pro_IX(:,2);

load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
a = zeros(size(powerunit_num_IX_PV,1),1);
b = zeros(size(powerunit_num_IX_onshorewind,1),1);
c = zeros(size(off_pro_IX,1),1);
for country =  1:1:192
    [m,n] = find(powerunit_num_IX_PV(:,5)==country);
    a(m,1) = region_ID(country,1);
    [m,n] = find(powerunit_num_IX_onshorewind(:,5)==country);
    b(m,1) = region_ID(country,1);
    [m,n] = find(off_pro_IX(:,2)==country);
    c(m,1) = region_ID(country,1);
end
powerunit_num_IX_PV(:,5)=a;
powerunit_num_IX_onshorewind(:,5)=b;
off_pro_IX(:,3)=c;
clear a
clear b
clear c
clear m

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的onshorewind电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat') % 按照成本排序&考虑建厂时间后保留的offshorewind电厂原始序号

[is,pos]=ismember(index_mineral_pv_time2,optpowerunit_PV(:,40));
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
optpowerunit_PV2 = optpowerunit_PV(pos,:);
powerunit_num_IX_PV2 = powerunit_num_IX_PV(pos,:);
powerunit_num_IX_PV_ori2 = powerunit_num_IX_PV_ori(pos,:);
lines_IX_PV2 = lines_IX_PV(pos,:);
clear optpowerunit_PV
clear powerunit_num_IX_PV
clear lines_IX_PV
clear powerunit_num_IX_PV_ori
optpowerunit_PV = optpowerunit_PV2;
powerunit_num_IX_PV = powerunit_num_IX_PV2;
powerunit_num_IX_PV_ori = powerunit_num_IX_PV_ori2;
lines_IX_PV = lines_IX_PV2;
clear optpowerunit_PV2
clear powerunit_num_IX_PV2
clear index_mineral_pv
clear lines_IX_PV2
clear powerunit_num_IX_PV_ori2

[is,pos]=ismember(index_mineral_ons_time2,optpowerunit_onshorewind(:,40));
optpowerunit_onshorewind2 = optpowerunit_onshorewind(pos,:);
powerunit_num_IX_onshorewind2 = powerunit_num_IX_onshorewind(pos,:);
lines_IX_onshorewind2 = lines_IX_onshorewind(pos,:);
powerunit_num_IX_onshorewind_ori2 = powerunit_num_IX_onshorewind_ori(pos,:);
clear optpowerunit_onshorewind
clear powerunit_num_IX_onshorewind
clear powerunit_num_IX_onshorewind_ori
clear lines_IX_onshorewind
optpowerunit_onshorewind = optpowerunit_onshorewind2;
powerunit_num_IX_onshorewind = powerunit_num_IX_onshorewind2;
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind_ori2;
lines_IX_onshorewind = lines_IX_onshorewind2;
clear optpowerunit_onshorewind2
clear powerunit_num_IX_onshorewind2
clear index_mineral_ons
clear lines_IX_onshorewind2
clear powerunit_num_IX_onshorewind_ori2

[is,pos]=ismember(index_mineral_off_time2,optpowerunit_offshorewind(:,40));
optpowerunit_offshorewind2 = optpowerunit_offshorewind(pos,:);
off_pro_IX2 = off_pro_IX(pos,:);
off_pro_IX_ori2 = off_pro_IX_ori(pos,:);
% powerunit_num_IX_offshorewind2 = powerunit_num_IX_offshorewind(pos,:);
lines_IX_offshorewind2 = lines_IX_offshorewind(pos,:);
clear optpowerunit_offshorewind
% clear powerunit_num_IX_offshorewind
clear off_pro_IX
clear off_pro_IX_ori
clear lines_IX_offshorewind
optpowerunit_offshorewind = optpowerunit_offshorewind2;
off_pro_IX = off_pro_IX2;
off_pro_IX_ori = off_pro_IX_ori2;
% powerunit_num_IX_offshorewind = powerunit_num_IX_offshorewind2;
lines_IX_offshorewind = lines_IX_offshorewind2;
clear optpowerunit_offshorewind2
clear off_pro_IX2
clear off_pro_IX_ori2
% clear powerunit_num_IX_offshorewind2
clear lines_IX_offshorewind2

%%
optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
clear optpowerunit_PV
clear optpowerunit_onshorewind
clear optpowerunit_offshorewind
lines_IX_offshorewind(:,12:18)=0;
lines_IX = [lines_IX_PV;lines_IX_onshorewind;lines_IX_offshorewind];
clear lines_IX_PV
clear lines_IX_onshorewind
clear lines_IX_offshorewind
powerunit_num_IX = [powerunit_num_IX_PV(:,5);powerunit_num_IX_onshorewind(:,5);off_pro_IX(:,3)];
powerunit_country_IX = [powerunit_num_IX_PV_ori;powerunit_num_IX_onshorewind_ori;off_pro_IX_ori];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
powerunit_IX_all(:,1)=IX;
optpowerunit_IX_all(:,1:40)=optpowerunit(IX,1:40);
lines_IX_IX(:,1:15)=lines_IX(IX,1:15);
powerunit_num_IX_IX(:,:)=powerunit_num_IX(IX,:);
powerunit_country_IX_IX(:,:)=powerunit_country_IX(IX,:);
clear optpowerunit
clear lines_IX
clear powerunit_num_IX
clear powerunit_country_IX
clear IX

%%
load('H:\Global PV and wind\ANS\country_needadjust2040_8_2s.mat'); % country_needadjust
load('H:\Global PV and wind\ANS\unitmin2040_8_2s_1.mat'); % unitmin
unitmin_all = unitmin;
learnprice_region = zeros(192,6);
cost_all2020_cou = zeros(192,1);
cost_ave_cou = zeros(192,1);
costmin_cou = zeros(192,1);
for iiaa = 1:size(country_needadjust,1)
    country = country_needadjust(iiaa);
    [m22,n]=find(powerunit_country_IX_IX==country & optpowerunit_IX_all(:,35)==2 & unitmin_all<=4);
    if ~isempty(m22)
        optpowerunit = [optpowerunit_IX_all(m22,:)];
        numpowerunit_pv(country,1) = size(optpowerunit,1);
        if ~isempty(optpowerunit)
            [B,IX]=sort(optpowerunit(:,20),1);
            numpowerunit = size(optpowerunit,1);
            for i=1:numpowerunit
                i2=IX(i);
                powerunit_IX(i,1)=i2;
                optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
            end
            
            power_country(country,1) = sum(optpowerunit_IX(:,1));
            
            % onshorewind
            if Country_classify(country)==1 % developed country 发达国家
                LR_onshorewind1 = 0.073; %0.18; % 0.37; %  0.37 PV module
                LR_onshorewind2= 0.1643;
            else
                LR_onshorewind1 = 0.122; %0.18; % 0.37; %  0.37 PV module
                LR_onshorewind2= 0.25;
            end
            CP0_onshorewind = CP_PV_ons_off2020(country,2); % 278324; % MW
            CP_onshorewind = zeros(numpowerunit,1);
            CP_onshorewind(1,1) = CP0_onshorewind;
            for i = 1: numpowerunit
                if optpowerunit_IX(i,35)==2
                    CP_onshorewind(i+1,1) =  CP_onshorewind(i,1) + optpowerunit_IX(i,30);
                else
                    CP_onshorewind(i+1,1) =  CP_onshorewind(i,1);
                end
            end
            CP0_onshorewind_cum = CP_onshorewind(2:end,1)-CP0_onshorewind;
            if CP0_onshorewind<CP_mean(2)
                CP0_onshorewind = CP_mean(2);
            end
            r_onshorewind = CP_onshorewind./CP0_onshorewind;
            r_onshorewind(r_onshorewind<1)=1;
            
            for i = 1:size(r_onshorewind,1)
                LR_onshorewind(i,1) = LR_onshorewind1;
                learnprice_onshorewind(i,1) = (r_onshorewind(i,1)).^(log(1- LR_onshorewind(i,1))/log(2));
                LR_onshorewind_2(i,1) = LR_onshorewind2;
                learnprice_onshorewind_2(i,1) = (r_onshorewind(i,1)).^(log(1- LR_onshorewind_2(i,1))/log(2));
            end
            clear LR_onshorewind
            clear LR_onshorewind_2
            learnprice_onshorewind(find(isnan(learnprice_onshorewind)==1))=1;
            learnprice_onshorewind_2(find(isnan(learnprice_onshorewind_2)==1))=1;
            learnprice_region(country,3)=learnprice_onshorewind(end);
            learnprice_region(country,4)=learnprice_onshorewind_2(end);
            for i = 1:size(CP_onshorewind,1)-1
                CP_onshorewind22(i,1) = CP_onshorewind(i+1,1)-CP_onshorewind(1,1);
                CP_onshorewind22(i,2) = (CP_onshorewind(i+1,1)-CP_onshorewind(1,1))/(CP_onshorewind(end,1)-CP_onshorewind(1,1));
            end
            cost_onshorewind = zeros(numpowerunit,19);
            if optpowerunit_IX(1,35)==2
                cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
            end
            for i = 2: numpowerunit
                if optpowerunit_IX(i,35)==2
                    cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
                else
                    cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
                end
            end
            clear CP_onshorewind
            clear CP0_onshorewind
            %%
            P_PV111= zeros(numpowerunit,1);
            P_onshorewind111= zeros(numpowerunit,1);
            P_offshorewind111= zeros(numpowerunit,1);
            for i = 1: numpowerunit
                if optpowerunit_IX(i,35)==1
                    P_PV111(i,1) =  optpowerunit_IX(i,1);
                else if optpowerunit_IX(i,35)==2
                        P_onshorewind111(i,1) =  optpowerunit_IX(i,1);
                    else if optpowerunit_IX(i,35)==3
                            P_offshorewind111(i,1) =  optpowerunit_IX(i,1);
                        end
                    end
                end
            end
            P_ratio_3type(:,1) = cumsum(P_PV111)./p_all2060_country(country,1);
            P_ratio_3type(:,2) = cumsum(P_onshorewind111)./p_all2060_country(country,1);
            P_ratio_3type(:,3) = cumsum(P_offshorewind111)./p_all2060_country(country,1);
            %
            optpowerunit1= zeros(numpowerunit,30);
            optpowerunit1(1,30)=optpowerunit_IX(1,30);
            optpowerunit1(1,1)=optpowerunit_IX(1,1);
            
            for i=2:numpowerunit
                optpowerunit1(i,30)=optpowerunit1(i-1,30)+optpowerunit_IX(i,30); % cumulative capacity
                optpowerunit1(i,1)=optpowerunit1(i-1,1)+optpowerunit_IX(i,1); % electricity generation, TWh/year
            end
            optpowerunit_IX(:,36)=optpowerunit1(:,30)/optpowerunit1(end,30);
            pe_ratio = optpowerunit1(:,1)/p_all2060_country(country,1);
            
            costmin=1e18;
            % costmax = 0;
            module_price_offshore = 1.008; % USD/W
            
            if size(optpowerunit_IX,1)>1
                if zerotracker_year_country(country)==2070
                    for p1=0:1:100
                        costmin1=1e18;
                        costmax1 = 0;
                        %     display(p1);
                        pene=zeros(5,2);
                        pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                        for p2=0:1:(100-p1)
                            pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040
                            for p3=0:1:(100-p1-p2)
                                pene(3,1)=(p1+p2+p3)/100; % peneratration ratio in 2040-2050
                                for p4=0:1:(100-p1-p2-p3)
                                    pene(4,1)=(p1+p2+p3+p4)/100; % peneratration ratio in 2050-2060
                                    pene(5,1)=1; % peneratration ratio in 2050-2060
                                    costall_PV=0;
                                    costall_onshorewind=0;
                                    costall_offshorewind=0;
                                    cumulativecapacity=0;
                                    unit1=1;
                                    for yy=1:5
                                        idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
                                        unit2(yy,1)=idx(1);
                                        if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                            if unit1==1
                                                % onshorewind
                                                costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14)+cost_onshorewind(unit2(yy,1),16)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD
                                            else
                                                % onshorewind
                                                costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*r_module_hardwarecost_ons*learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*(1-r_module_hardwarecost_ons)*learnprice_onshorewind_2(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD
                                            end
                                            pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                            pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                            unit1=unit2(yy,1)+1;
                                        end
                                    end
                                    
                                    for aaa= 1:1:4
                                        asdf = pene(aaa+1,2); % peneratration ratio CP
                                        if asdf==0
                                            pene(aaa+1,2)= pene(aaa,2);
                                        end
                                    end
                                    
                                    costall = costall_PV + costall_onshorewind + costall_offshorewind;
                                    if costall<costmin
                                        costmin=costall;
                                        penemin=pene;
                                        penemin_real=pene_real;
                                        unitmin=unit2;
                                        asdf=[penemin;costmin*ones(1,2)];
                                        asdf_real=[penemin_real;costmin*ones(1,2)];
                                    end
                                    if unique(pene(:,1))==1
                                        cost_all2020 = costall;
                                    end
                                    if pene(1,1)==0.2 && pene(2,1)==0.4 && pene(3,1)==0.6 && pene(4,1)==0.8 && pene(5,1)==1
                                        cost_ave = costall;
                                    end
                                end
                            end
                        end
                    end
                end
                if zerotracker_year_country(country)==2060
                    for p1=0:1:100
                        costmin1=1e18;
                        costmax1 = 0;
                        %     display(p1);
                        pene=zeros(4,2);
                        pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                        for p2=0:1:(100-p1)
                            pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040
                            for p3=0:1:(100-p1-p2)
                                pene(3,1)=(p1+p2+p3)/100; % peneratration ratio in 2040-2050
                                pene(4,1)=1; % peneratration ratio in 2050-2060
                                costall_PV=0;
                                costall_onshorewind=0;
                                costall_offshorewind=0;
                                cumulativecapacity=0;
                                unit1=1;
                                for yy=1:4
                                    idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
                                    unit2(yy,1)=idx(1);
                                    if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                        if unit1==1
                                            % onshorewind
                                            costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14)+cost_onshorewind(unit2(yy,1),16)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD
                                        else
                                            % onshorewind
                                            costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*r_module_hardwarecost_ons*learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*(1-r_module_hardwarecost_ons)*learnprice_onshorewind_2(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD
                                        end
                                        pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                        pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                        unit1=unit2(yy,1)+1;
                                    end
                                end
                                
                                for aaa= 1:1:3
                                    asdf = pene(aaa+1,2); % peneratration ratio CP
                                    if asdf==0
                                        pene(aaa+1,2)= pene(aaa,2);
                                    end
                                end
                                
                                costall = costall_PV + costall_onshorewind + costall_offshorewind;
                                if costall<costmin
                                    costmin=costall;
                                    penemin=pene;
                                    penemin_real=pene_real;
                                    unitmin=unit2;
                                    asdf=[penemin;costmin*ones(1,2)];
                                    asdf_real=[penemin_real;costmin*ones(1,2)];
                                end
                                if unique(pene(:,1))==1
                                    cost_all2020 = costall;
                                end
                                if pene(1,1)==0.25 && pene(2,1)==0.5 && pene(3,1)==0.75 && pene(4,1)==1
                                    cost_ave = costall;
                                end
                            end
                        end
                    end
                end
                if zerotracker_year_country(country)==2050
                    for p1=0:1:100
                        costmin1=1e18;
                        costmax1 = 0;
                        pene=zeros(3,2);
                        pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                        for p2=0:1:(100-p1)
                            pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040
                            pene(3,1)=1; % peneratration ratio in 2040-2050
                            costall_PV=0;
                            costall_onshorewind=0;
                            costall_offshorewind=0;
                            cumulativecapacity=0;
                            unit1=1;
                            for yy=1:3
                                idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
                                unit2(yy,1)=idx(1);
                                if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                    if unit1==1
                                        % onshorewind
                                        costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14)+cost_onshorewind(unit2(yy,1),16)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD
                                    else
                                        % onshorewind
                                        costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*r_module_hardwarecost_ons*learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*(1-r_module_hardwarecost_ons)*learnprice_onshorewind_2(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD
                                    end
                                    pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                    pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                    unit1=unit2(yy,1)+1;
                                end
                            end
                            
                            for aaa= 1:1:2
                                asdf = pene(aaa+1,2); % peneratration ratio CP
                                if asdf==0
                                    pene(aaa+1,2)= pene(aaa,2);
                                end
                            end
                            
                            costall = costall_PV + costall_onshorewind + costall_offshorewind;
                            if costall<costmin
                                costmin=costall;
                                penemin=pene;
                                penemin_real=pene_real;
                                unitmin=unit2;
                                asdf=[penemin;costmin*ones(1,2)];
                                asdf_real=[penemin_real;costmin*ones(1,2)];
                            end
                            if unique(pene(:,1))==1
                                cost_all2020 = costall;
                            end
                            if pene(1,1)==0.33 && pene(2,1)==0.67 && pene(3,1)==1
                                cost_ave = costall;
                            end
                        end
                    end
                end
                if zerotracker_year_country(country)==2040
                    for p1=0:1:100
                        costmin1=1e18;
                        costmax1 = 0;
                        pene=zeros(2,2);
                        pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                        pene(2,1)=1; % peneratration ratio in 2030-2040
                        costall_PV=0;
                        costall_onshorewind=0;
                        costall_offshorewind=0;
                        cumulativecapacity=0;
                        unit1=1;
                        for yy=1:2
                            idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
                            unit2(yy,1)=idx(1);
                            if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                if unit1==1
                                    % onshorewind
                                    costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14)+cost_onshorewind(unit2(yy,1),16)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD
                                else
                                    % onshorewind
                                    costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*r_module_hardwarecost_ons*learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*(1-r_module_hardwarecost_ons)*learnprice_onshorewind_2(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD
                                end
                                pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                unit1=unit2(yy,1)+1;
                            end
                        end
                        
                        for aaa= 1
                            asdf = pene(aaa+1,2); % peneratration ratio CP
                            if asdf==0
                                pene(aaa+1,2)= pene(aaa,2);
                            end
                        end
                        
                        costall = costall_PV + costall_onshorewind + costall_offshorewind;
                        if costall<costmin
                            costmin=costall;
                            penemin=pene;
                            penemin_real=pene_real;
                            unitmin=unit2;
                            asdf=[penemin;costmin*ones(1,2)];
                            asdf_real=[penemin_real;costmin*ones(1,2)];
                        end
                        if unique(pene(:,1))==1
                            cost_all2020 = costall;
                        end
                        if pene(1,1)==0.5 && pene(2,1)==1
                            cost_ave = costall;
                        end
                    end
                end
                
                
                peneminA = round(penemin(:,1)*100);
                costmin=1e18;
                for p1=0:0.5:peneminA(1,1)
                    costmin1=1e18;
                    costmax1 = 0;
                    pene=zeros(4,2);
                    pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                    idx=find(optpowerunit_IX(:,36)>=pene(1,1)); % percentage of energy
                    unit2(1,1)=idx(1);
                    for p2=peneminA(1,1)-p1
                        pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040
                        idx=find(optpowerunit_IX(:,36)>=pene(2,1)); % percentage of energy
                        unit2(2,1)=idx(1);
                        for p3=peneminA(1,1)-p1-p2:1:peneminA(2,1)-p1-p2
                            pene(3,1)=(p1+p2+p3)/100; % peneratration ratio in 2040-2050
                            idx=find(optpowerunit_IX(:,36)>=pene(3,1)); % percentage of energy
                            unit2(3,1)=idx(1);
                            pene(4,1)= 1; % peneratration ratio in 2050-2060
                            unit2(4,1)=size(optpowerunit_IX,1);
                            costall_PV=0;
                            costall_onshorewind=0;
                            costall_offshorewind=0;
                            cumulativecapacity=0;
                            unit1=1;
                            for yy=1:4
                                if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                    if unit1==1
                                        % onshorewind
                                        costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14)+cost_onshorewind(unit2(yy,1),16)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD
                                    else
                                        % onshorewind
                                        costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*r_module_hardwarecost_ons*learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14))*(1-r_module_hardwarecost_ons)*learnprice_onshorewind_2(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD
                                    end
                                    pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                    pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                    unit1=unit2(yy,1)+1;
                                end
                            end
                            
                            for aaa= 1:1:3
                                asdf = pene(aaa+1,2); % peneratration ratio CP
                                if asdf==0
                                    pene(aaa+1,2)= pene(aaa,2);
                                end
                            end
                            
                            costall = costall_PV + costall_onshorewind + costall_offshorewind;
                            if costall<costmin
                                costmin=costall;
                                penemin=pene;
                                penemin_real=pene_real;
                                unitmin=unit2;
                                asdf=[penemin;costmin*ones(1,2)];
                            end
                            
                        end
                    end
                end
                
                cost_all2020_cou(country,1) = cost_all2020;
                cost_ave_cou(country,1) = cost_ave;
                costmin_cou(country,1) = costmin;
            end
            
            if size(optpowerunit_IX,1)==1
                unitmin=ones((zerotracker_year_country(country)-2020)/5,1);
            end
            
            %%
            unit_PV=[];
            unit_onshorewind=[];
            unit_offshorewind=[];
            for i = 1:size(powerunit_IX,1)
                %     aa=powerunit_IX(i);
                aa=i;
                bb=optpowerunit_IX(aa,40);
                if optpowerunit_IX(aa,35) ==1
                    unit_PV(bb,1)=i;% 排序之前的
                end
                if optpowerunit_IX(aa,35) ==2
                    unit_onshorewind(bb,1)=i;% 排序之前的
                end
                if optpowerunit_IX(aa,35) ==3
                    unit_offshorewind(bb,1)=i;% 排序之前的
                end
                %     i
            end
            
            if zerotracker_year_country(country)==2040
                if ~isempty(unit_PV)
                    unit_PV(unit_PV<=unitmin(1,1)& unit_PV>0)=1;
                    unit_PV(unit_PV<=unitmin(2,1)& unit_PV>unitmin(1,1))=2;
                    unit_PV(unit_PV<=unitmin(3,1)& unit_PV>unitmin(2,1))=3;
                    unit_PV(unit_PV<=unitmin(4,1)& unit_PV>unitmin(3,1))=4;
                end
                if ~isempty(unit_onshorewind)
                    unit_onshorewind(unit_onshorewind<=unitmin(1,1)& unit_onshorewind>0)=1;
                    unit_onshorewind(unit_onshorewind<=unitmin(2,1)& unit_onshorewind>unitmin(1,1))=2;
                    unit_onshorewind(unit_onshorewind<=unitmin(3,1)& unit_onshorewind>unitmin(2,1))=3;
                    unit_onshorewind(unit_onshorewind<=unitmin(4,1)& unit_onshorewind>unitmin(3,1))=4;
                end
                if ~isempty(unit_offshorewind)
                    unit_offshorewind(unit_offshorewind<=unitmin(1,1)& unit_offshorewind>0)=1;
                    unit_offshorewind(unit_offshorewind<=unitmin(2,1)& unit_offshorewind>unitmin(1,1))=2;
                    unit_offshorewind(unit_offshorewind<=unitmin(3,1)& unit_offshorewind>unitmin(2,1))=3;
                    unit_offshorewind(unit_offshorewind<=unitmin(4,1)& unit_offshorewind>unitmin(3,1))=4;
                end
            end
            
            
            [m,n]=find(unit_PV~=0);
            if ~isempty(m)
                unit_PV_global(m,1)=unit_PV(m,1);
            end
            [m,n]=find(unit_onshorewind~=0);
            if ~isempty(m)
                unit_onshorewind_global(m,1)=unit_onshorewind(m,1);
            end
            [m,n]=find(unit_offshorewind~=0);
            if ~isempty(m)
                unit_offshorewind_global(m,1)=unit_offshorewind(m,1);
            end
            
            clear powerunit_IX;
            clear optpowerunit_IX;
            clear learnprice_PV;
            clear learnprice_PV_2;
            clear CP_PV22;
            clear learnprice_onshorewind;
            clear learnprice_onshorewind_2;
            clear CP_onshorewind22;
            clear CP_offshorewind22;
            clear learnprice_offshorewind;
            clear learnprice_offshorewind_2;
            clear P_ratio_3type
            clear unit_PV;
            clear unit_onshorewind;
            clear unit_offshorewind;
            clear penex;
            clear pene_real;
            clear pene;
            clear penex_real;
            clear penemin_real;
            clear P_ratio_3type;clear pe_ratio
            
            clear pene_reg;
            clear costmax1;
            clear costmin1
            clear penemax1
            clear penemin1
            
        end
    end
end
save('H:\Global PV and wind\ANS\unit_onshorewind_global_2040.mat','unit_onshorewind_global');


