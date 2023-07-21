% 实际各国learning rate，当为nan和inf时和全球统一
tic
clear;
load('H:\Global PV and wind\Data\r_module_hardwarecost_pv.mat');

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
p_all2060=sum(powerdemand_monhour2070_ele288)'; % 	TWh/year
clear powerdemand_monhour2060_ele

load('H:\Global PV and wind\Data\ID_pro3.mat') % FID  FIRST_ID_0  ID_country120_0214	FIRST_ID_1
for i = 1:1:192
    [m,n]=find(ID_pro(:,3)==i);
    FID_pro = unique(ID_pro(m,1));
    p_all2060_country(i,1) =sum(p_all2060(FID_pro+1,1));
    clear FID_pro
    i
end
clear p_all2060
clear ID_pro

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\optpowerunit_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_num_IX_PV_100GW_3_2_all2_5%_inilow.mat'); %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % 电厂编号

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\optpowerunit_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Onshore wind_power potential\ANS\powerunit_num_IX_onshorewind_100GW_3_2_all_5%_inilow.mat'); %
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % 电厂编号


load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\optpowerunit_offshorewind_100GW_county_5%.mat'); %
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_IX_offshorewind_100GW_county_5%.mat'); %
% load('H:\world\code\数据处理\unitid_lcoe_offshorewind1227_2_xz_100GW_county_5%.mat');
% unitid_lcoe_offshorewind = unitid_lcoe;
% clear unitid_lcoe
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat'); %
% 1.power plant ID;2 country;3 pro;4 county

% clear off_pro_IX
optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % 电厂编号

r_pv_world = sum(optpowerunit_PV(:,30))/sum(CP_PV_ons_off2020(:,1))
r_ons_world = sum(optpowerunit_onshorewind(:,30))/sum(CP_PV_ons_off2020(:,2))
r_off_world = sum(optpowerunit_offshorewind(:,30))/sum(CP_PV_ons_off2020(:,3))

for country =  1:1:192
    clear powerunit_IX;
    [id_m_pv,n]=find(powerunit_num_IX_PV(:,5)==country);
    [id_m_ons,n]=find(powerunit_num_IX_onshorewind(:,5)==country);
    [id_m_off,n]=find(off_pro_IX(:,2)==country);
    a(country,1) = sum(optpowerunit_PV(id_m_pv,1));
    a(country,2) = sum(optpowerunit_onshorewind(id_m_ons,1));
    a(country,3) = sum(optpowerunit_offshorewind(id_m_off,1));
    country
end
Ph_country = sum(a,2);
% save('H:\world\code\optimized_Time\Ph_country_8_2070.mat','Ph_country'); % TWh/year

powerunit_num_IX_PV_ori = powerunit_num_IX_PV(:,5);
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind(:,5);
off_pro_IX_ori = off_pro_IX(:,2);

%
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
off_pro_IX(:,2)=c;
clear a
clear b
clear c
clear m

%%
unit_PV_global = zeros(size(optpowerunit_PV,1),1);
unit_onshorewind_global = zeros(size(optpowerunit_onshorewind,1),1);
unit_offshorewind_global = zeros(size(optpowerunit_offshorewind,1),1);

%% 根据mineral计算的太阳能和风能在40年内最多建厂数目
load('H:\Global PV and wind\ANS\index_mineral_pv_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat') % 按照成本排序后保留的PV电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_ons_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat') % 按照成本排序后保留的onshorewind电厂原始序号
load('H:\Global PV and wind\ANS\index_mineral_off_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat') % 按照成本排序后保留的offshorewind电厂原始序号

[is,pos]=ismember(index_mineral_pv,optpowerunit_PV(:,40));
% is是与B大小一致的向量，如果在A中为1，不在为0
% pos是B中元素如果在A中出现，出现的位置。
optpowerunit_PV2 = optpowerunit_PV(pos,:);
powerunit_num_IX_PV2 = powerunit_num_IX_PV(pos,:);
powerunit_num_IX_PV_ori2 = powerunit_num_IX_PV_ori(pos,:);
clear optpowerunit_PV
clear powerunit_num_IX_PV
clear powerunit_num_IX_PV_ori
optpowerunit_PV = optpowerunit_PV2;
powerunit_num_IX_PV = powerunit_num_IX_PV2;
powerunit_num_IX_PV_ori = powerunit_num_IX_PV_ori2;
clear powerunit_num_IX_PV_ori2
clear optpowerunit_PV2
clear powerunit_num_IX_PV2
clear index_mineral_pv

[is,pos]=ismember(index_mineral_ons,optpowerunit_onshorewind(:,40));
optpowerunit_onshorewind2 = optpowerunit_onshorewind(pos,:);
powerunit_num_IX_onshorewind2 = powerunit_num_IX_onshorewind(pos,:);
powerunit_num_IX_onshorewind_ori2 = powerunit_num_IX_onshorewind_ori(pos,:);
clear powerunit_num_IX_onshorewind_ori
clear optpowerunit_onshorewind
clear powerunit_num_IX_onshorewind
optpowerunit_onshorewind = optpowerunit_onshorewind2;
powerunit_num_IX_onshorewind = powerunit_num_IX_onshorewind2;
powerunit_num_IX_onshorewind_ori = powerunit_num_IX_onshorewind_ori2;
clear optpowerunit_onshorewind2
clear powerunit_num_IX_onshorewind2
clear index_mineral_ons
clear powerunit_num_IX_onshorewind_ori2

[is,pos]=ismember(index_mineral_off,optpowerunit_offshorewind(:,40));
optpowerunit_offshorewind2 = optpowerunit_offshorewind(pos,:);
off_pro_IX2 = off_pro_IX(pos,:);
off_pro_IX_ori2 = off_pro_IX_ori(pos,:);
clear off_pro_IX_ori
clear optpowerunit_offshorewind
clear off_pro_IX
optpowerunit_offshorewind = optpowerunit_offshorewind2;
off_pro_IX = off_pro_IX2;
off_pro_IX_ori = off_pro_IX_ori2;
clear off_pro_IX_ori2
clear optpowerunit_offshorewind2
clear powerunit_num_IX_onshorewind2
clear off_pro_IX2

sum(optpowerunit_PV(:,1))+sum(optpowerunit_onshorewind(:,1))+sum(optpowerunit_offshorewind(:,1))
% powerunit_country_IX = [powerunit_num_IX_PV_ori;powerunit_num_IX_onshorewind_ori;off_pro_IX_ori];
%%
% learning rate :2020-B Sweerts-Joule-Evaluating the role of unit size in learning-by-doing of energy technologies
learnprice_region = zeros(192,6);
cost_all2020_cou = zeros(192,1);
cost_ave_cou = zeros(192,1);
costmin_cou = zeros(192,1);
% unitid_s1_PV=zeros(21600,43200);
% unitid_s1_onshorewind=zeros(21600,43200);
% unitid_s1_offshorewind=zeros(21600,43200);
for country = 1:1:192
    [id_m_pv,n]=find(powerunit_num_IX_PV_ori(:,1)==country);
    
    optpowerunit = [optpowerunit_PV(id_m_pv,:)];
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
        
        % PV
        if Country_classify(country)==1 % developed country 发达国家
            LR_PV1 = 0.18; % PV module
            LR_PV2 = 0.1643; % BOS
        else
            LR_PV1 = 0.26; %0.18; % 0.37; %  0.37 PV module
            LR_PV2 = 0.25; % BOS
        end
        CP0_PV = CP_PV_ons_off2020(country,1); %  cumulative capacity potential in China in 2020, MW
        CP_PV = zeros(numpowerunit,1);
        CP_PV(1,1) = CP0_PV;
        for i = 1: numpowerunit
            if optpowerunit_IX(i,35)==1
                CP_PV(i+1,1) =  CP_PV(i,1) + optpowerunit_IX(i,30);
            else
                CP_PV(i+1,1) =  CP_PV(i,1);
            end
        end
        CP0_PV_cum = CP_PV(2:end,1)-CP0_PV;
        if CP0_PV<CP_mean(1)
            CP0_PV = CP_mean(1);
        end
        r_pv = CP_PV./CP0_PV;
        r_pv(r_pv<1)=1;
        
        for i = 1:size(r_pv,1)
            LR_PV(i,1) = LR_PV1;
            learnprice_PV(i,1) = (r_pv(i)).^(log(1- LR_PV(i,1))/log(2));
            LR_PV_2(i,1) = LR_PV2;
            learnprice_PV_2(i,1) = (r_pv(i)).^(log(1- LR_PV_2(i,1))/log(2));
        end
        clear LR_PV
        clear LR_PV_2
        learnprice_PV(find(isnan(learnprice_PV)==1))=1;
        learnprice_PV_2(find(isnan(learnprice_PV_2)==1))=1;
        learnprice_region(country,1)=learnprice_PV(end);
        learnprice_region(country,2)=learnprice_PV_2(end);
        for i = 1:size(CP_PV,1)-1
            CP_PV22(i,1) = CP_PV(i+1,1)-CP_PV(1,1);
            CP_PV22(i,2) = (CP_PV(i+1,1)-CP_PV(1,1))/(CP_PV(end,1)-CP_PV(1,1));
        end
        cost_PV = zeros(numpowerunit,19);
        cost_PV2 = zeros(numpowerunit,19);
        if optpowerunit_IX(1,35)==1
            cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
            cost_PV2(1,11:19) =  optpowerunit_IX(1,11:19);
        end
        for i = 2: numpowerunit
            if optpowerunit_IX(i,35)==1
                cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
                cost_PV2(i,11:19) =  optpowerunit_IX(i,11:19);
            else
                cost_PV(i,11:19) =  cost_PV(i-1,11:19);
                cost_PV2(i,11:19) =  0;
            end
        end
        clear CP_PV
        clear CP0_PV
        
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
            for p1=0:1:100
                costmin1=1e18;
                costmax1 = 0;
                pene=zeros(5,2);
                pene(1,1)=p1/100; % peneratration ratio in 2020-2030
                idx=find(optpowerunit_IX(:,36)>=pene(1,1)); % percentage of energy
                unit2(1,1)=idx(1);
                for p2=0:1:(100-p1)
                    pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040
                    idx=find(optpowerunit_IX(:,36)>=pene(2,1)); % percentage of energy
                    unit2(2,1)=idx(1);
                    for p3=0:1:(100-p1-p2)
                        pene(3,1)=(p1+p2+p3)/100; % peneratration ratio in 2040-2050
                        idx=find(optpowerunit_IX(:,36)>=pene(3,1)); % percentage of energy
                        unit2(3,1)=idx(1);
                        for p4=0:1:(100-p1-p2-p3)
                            pene(4,1)=(p1+p2+p3+p4)/100; % peneratration ratio in 2050-2060
                            pene(5,1)= 1; % peneratration ratio in 2050-2060
                            idx=find(optpowerunit_IX(:,36)>=pene(4,1)); % percentage of energy
                            unit2(4,1)=idx(1);
                            unit2(5,1)=size(optpowerunit_IX,1);                                     costall_PV=0;
                                costall_onshorewind=0;
                                costall_offshorewind=0;
                                cumulativecapacity=0;
                                unit1=1;
                                for yy=1:5
%                                     idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
%                                     unit2(yy,1)=idx(1);
                                    if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                        if unit1==1
                                            % PV
                                            costall_PV=costall_PV+cost_PV(unit2(yy,1),14)+cost_PV(unit2(yy,1),16) +sum(cost_PV(unit2(yy,1),11:13),2)+ cost_PV(unit2(yy,1),15); % cost million USD
                                        else
                                            % PV
                                            costall_PV=costall_PV+(cost_PV(unit2(yy,1),14)-cost_PV(unit1-1,14))*r_module_hardwarecost_pv(country,1)*learnprice_PV(unit1)+(cost_PV(unit2(yy,1),14)-cost_PV(unit1-1,14))*(1-r_module_hardwarecost_pv(country,1))*learnprice_PV_2(unit1)+(cost_PV(unit2(yy,1),16)-cost_PV(unit1-1,16))+(sum(cost_PV(unit2(yy,1),11:13),2)-sum(cost_PV(unit1-1,11:13),2))*learnprice_PV_2(unit1)+(cost_PV(unit2(yy,1),15)-cost_PV(unit1-1,15))*learnprice_PV_2(unit1); % cost million USD
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
                                
                                costall = costall_PV;
                                if costall<costmin
                                    costmin=costall;
                                    penemin=pene;
                                    penemin_real=pene_real;
                                    unitmin=unit2;
                                    asdf=[penemin;costmin*ones(1,2)];
%                                     asdf_real=[penemin_real;costmin*ones(1,2)];
                                end
%                                 if unique(pene(:,1))==1
%                                     cost_all2020 = costall;
%                                 end
%                                 if pene(1,1)==0.2 && pene(2,1)==0.4 && pene(3,1)==0.6 && pene(4,1)==0.8 && pene(5,1)==1
%                                     cost_ave = costall;
%                                 end
                                
                            end
                        end
                    end
            end
            peneminA = round(penemin(:,1)*100);
            costmin=1e18;
            for p1=0:0.5:peneminA(1,1)
                costmin1=1e18;
                costmax1 = 0;
                pene=zeros(10,2);
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
                        for p4=peneminA(2,1)-p1-p2-p3
                            pene(4,1)=(p1+p2+p3+p4)/100; % peneratration ratio in 2050-2060
                            idx=find(optpowerunit_IX(:,36)>=pene(4,1)); % percentage of energy
                            unit2(4,1)=idx(1);
                            for p5=peneminA(2,1)-p1-p2-p3-p4:1:peneminA(3,1)-p1-p2-p3-p4
                                pene(5,1)=(p1+p2+p3+p4+p5)/100; % peneratration ratio in 2050-2060
                                idx=find(optpowerunit_IX(:,36)>=pene(5,1)); % percentage of energy
                                unit2(5,1)=idx(1);
                                for p6=peneminA(3,1)-p1-p2-p3-p4-p5
                                    pene(6,1)=(p1+p2+p3+p4+p5+p6)/100; % peneratration ratio in 2050-2060
                                    idx=find(optpowerunit_IX(:,36)>=pene(6,1)); % percentage of energy
                                    unit2(6,1)=idx(1);
                                    for p7=peneminA(3,1)-p1-p2-p3-p4-p5-p6:1:peneminA(4,1)-p1-p2-p3-p4-p5-p6
                                        pene(7,1)=(p1+p2+p3+p4+p5+p6+p7)/100; % peneratration ratio in 2050-2060
                                        idx=find(optpowerunit_IX(:,36)>=pene(7,1)); % percentage of energy
                                        unit2(7,1)=idx(1);
                                        for p8=peneminA(4,1)-p1-p2-p3-p4-p5-p6-p7
                                            pene(8,1)=(p1+p2+p3+p4+p5+p6+p7+p8)/100; % peneratration ratio in 2050-2060
                                            idx=find(optpowerunit_IX(:,36)>=pene(8,1)); % percentage of energy
                                            unit2(8,1)=idx(1);
                                            for p9=peneminA(4,1)-p1-p2-p3-p4-p5-p6-p7-p8:1:peneminA(5,1)-p1-p2-p3-p4-p5-p6-p7-p8
                                                pene(9,1)=(p1+p2+p3+p4+p5+p6+p7+p8+p9)/100; % peneratration ratio in 2050-2060
                                                pene(10,1)= 1; % peneratration ratio in 2050-2060
                                                idx=find(optpowerunit_IX(:,36)>=pene(9,1)); % percentage of energy
                                                unit2(9,1)=idx(1);
                                                unit2(10,1)=size(optpowerunit_IX,1);
                                                costall_PV=0;
                                                costall_onshorewind=0;
                                                costall_offshorewind=0;
                                                cumulativecapacity=0;
                                                unit1=1;
                                                for yy=1:10
                                                    if unit2(yy,1)>=unit1 && unit2(yy,1)~=1
                                                        if unit1==1
                                                            % PV
                                                            costall_PV=costall_PV+cost_PV(unit2(yy,1),14)+cost_PV(unit2(yy,1),16) +sum(cost_PV(unit2(yy,1),11:13),2)+ cost_PV(unit2(yy,1),15); % cost million USD
                                                        else
                                                            % PV
                                                            costall_PV=costall_PV+(cost_PV(unit2(yy,1),14)-cost_PV(unit1-1,14))*r_module_hardwarecost_pv(country,1)*learnprice_PV(unit1)+(cost_PV(unit2(yy,1),14)-cost_PV(unit1-1,14))*(1-r_module_hardwarecost_pv(country,1))*learnprice_PV_2(unit1)+(cost_PV(unit2(yy,1),16)-cost_PV(unit1-1,16))+(sum(cost_PV(unit2(yy,1),11:13),2)-sum(cost_PV(unit1-1,11:13),2))*learnprice_PV_2(unit1)+(cost_PV(unit2(yy,1),15)-cost_PV(unit1-1,15))*learnprice_PV_2(unit1); % cost million USD
                                                        end
                                                        pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio CP
                                                        pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % 发电量占2060年总需电量的比值
                                                        unit1=unit2(yy,1)+1;
                                                    end
                                                end
                                                
                                                for aaa= 1:1:9
                                                    asdf = pene(aaa+1,2); % peneratration ratio CP
                                                    if asdf==0
                                                        pene(aaa+1,2)= pene(aaa,2);
                                                    end
                                                end
                                                
                                                costall = costall_PV;
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
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if size(optpowerunit_IX,1)==1
            unitmin=ones(10,1);
%             filename = strcat('H:\Global PV and wind\ANS\3unitmin_region',num2str(country),'_county0811_all_pv_pro.mat');
%             save(filename, 'unitmin','-v7.3');
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
        
        if ~isempty(unit_PV)
            unit_PV(unit_PV<=unitmin(1,1)& unit_PV>0)=1;
            unit_PV(unit_PV<=unitmin(2,1)& unit_PV>unitmin(1,1))=2;
            unit_PV(unit_PV<=unitmin(3,1)& unit_PV>unitmin(2,1))=3;
            unit_PV(unit_PV<=unitmin(4,1)& unit_PV>unitmin(3,1))=4;
            unit_PV(unit_PV<=unitmin(5,1)& unit_PV>unitmin(4,1))=5;
            unit_PV(unit_PV<=unitmin(6,1)& unit_PV>unitmin(5,1))=6;
            unit_PV(unit_PV<=unitmin(7,1)& unit_PV>unitmin(6,1))=7;
            unit_PV(unit_PV<=unitmin(8,1)& unit_PV>unitmin(7,1))=8;
            unit_PV(unit_PV<=unitmin(9,1)& unit_PV>unitmin(8,1))=9;
            unit_PV(unit_PV<=unitmin(10,1)& unit_PV>unitmin(9,1))=10;
        end
        if ~isempty(unit_onshorewind)
            unit_onshorewind(unit_onshorewind<=unitmin(1,1)& unit_onshorewind>0)=1;
            unit_onshorewind(unit_onshorewind<=unitmin(2,1)& unit_onshorewind>unitmin(1,1))=2;
            unit_onshorewind(unit_onshorewind<=unitmin(3,1)& unit_onshorewind>unitmin(2,1))=3;
            unit_onshorewind(unit_onshorewind<=unitmin(4,1)& unit_onshorewind>unitmin(3,1))=4;
            unit_onshorewind(unit_onshorewind<=unitmin(5,1)& unit_onshorewind>unitmin(4,1))=5;
            unit_onshorewind(unit_onshorewind<=unitmin(6,1)& unit_onshorewind>unitmin(5,1))=6;
            unit_onshorewind(unit_onshorewind<=unitmin(7,1)& unit_onshorewind>unitmin(6,1))=7;
            unit_onshorewind(unit_onshorewind<=unitmin(8,1)& unit_onshorewind>unitmin(7,1))=8;
            unit_onshorewind(unit_onshorewind<=unitmin(9,1)& unit_onshorewind>unitmin(8,1))=9;
            unit_onshorewind(unit_onshorewind<=unitmin(10,1)& unit_onshorewind>unitmin(9,1))=10;
        end
        if ~isempty(unit_offshorewind)
            unit_offshorewind(unit_offshorewind<=unitmin(1,1)& unit_offshorewind>0)=1;
            unit_offshorewind(unit_offshorewind<=unitmin(2,1)& unit_offshorewind>unitmin(1,1))=2;
            unit_offshorewind(unit_offshorewind<=unitmin(3,1)& unit_offshorewind>unitmin(2,1))=3;
            unit_offshorewind(unit_offshorewind<=unitmin(4,1)& unit_offshorewind>unitmin(3,1))=4;
            unit_offshorewind(unit_offshorewind<=unitmin(5,1)& unit_offshorewind>unitmin(4,1))=5;
            unit_offshorewind(unit_offshorewind<=unitmin(6,1)& unit_offshorewind>unitmin(5,1))=6;
            unit_offshorewind(unit_offshorewind<=unitmin(7,1)& unit_offshorewind>unitmin(6,1))=7;
            unit_offshorewind(unit_offshorewind<=unitmin(8,1)& unit_offshorewind>unitmin(7,1))=8;
            unit_offshorewind(unit_offshorewind<=unitmin(9,1)& unit_offshorewind>unitmin(8,1))=9;
            unit_offshorewind(unit_offshorewind<=unitmin(10,1)& unit_offshorewind>unitmin(9,1))=10;
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
        country
    end
    
end
save('H:\Global PV and wind\ANS\3power_country_region_county0811_all_country_5%_LRglobal_inilow_pv_pro20702_8_2070.mat','power_country'); % TWh/year
save('H:\Global PV and wind\ANS\3unit_PV_region_county0811_all_country_5%_LRglobal_inilow_pv_pro20702_8_2070.mat','unit_PV_global');
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_country_5%_LRglobal_inilow_pv_pro20702_8_2070.mat','learnprice_region'); % 10 regions * 6
save('H:\Global PV and wind\ANS\numpowerunit_pv_pro20702_8_2070.mat','numpowerunit_pv'); % numpowerunit_pv,各国PV电厂个数
save('H:\Global PV and wind\ANS\cost_all2020_cou_pv_pro20702_8_2070.mat','cost_all2020_cou');
save('H:\Global PV and wind\ANS\cost_ave_cou_pv_pro20702_8_2070.mat','cost_ave_cou'); %
save('H:\Global PV and wind\ANS\costmin_cou_pv_pro20702_8_2070.mat','costmin_cou'); %
