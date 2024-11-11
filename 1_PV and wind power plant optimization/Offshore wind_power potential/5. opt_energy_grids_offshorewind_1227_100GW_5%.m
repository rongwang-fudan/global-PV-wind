tic
clear;
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\ID_dist_offshore120_1227_xz2_xz_county.mat'); % $/W
load('H:\global-PV-wind\Data\S_offshorewind_1227_xzz_county.mat');
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\costunits_offshorewind1227_xz2_xz_county.dat','-mat');
% 1 pro; 2 Powerunit_ID; 3 CP,kW; 4 Cost,million dollar; 5 kwh/year;
% 6 S km2; 7 mean wind speeed,m/s
load('H:\global-PV-wind\Data\dist_id_pro120_1227_xz_county.mat'); % pro ID
load('H:\global-PV-wind\Data\initialcost_ratio_country_0111low_off.mat')  % offshorewind的成本与China的对比

dist_id_pro120_1 = ones(21600,43200)*(-1);
[m,n]=find(S_offshorewind~=0);
dist_id_pro120_1(sub2ind(size(dist_id_pro120_1), m, n))= dist_id_pro120(sub2ind(size(dist_id_pro120), m, n));
clear dist_id_pro120
clear S_offshorewind

module_price = 1.008; % USD/W
rmb2us=1/6.8967; % RMB to USD2019
lifetime_power=25;
load('H:\global-PV-wind\Data\discount_country.mat')
discount = discount/100;
discount(35) = 0.05;

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
OMratio_substation=0.03;
load('H:\global-PV-wind\Data\fossilfuel_emissionfactor.mat')  %kg CO2/kWh
CO2_C=0.2727;

load('H:\global-PV-wind\Data\dist_id_country120_1227_xz_county.mat');
bb=unique(dist_id_country120);
bb1 = bb(2:end,1);
clear bb
load('H:\global-PV-wind\Data\pro_county.mat'); % 第一列是pro的ID，第二列是county

% figure
cmap=jet(12);
offshorewind_power_lcoemin = zeros(21600,43200);
unitid_lcoe = ones(21600,43200)*(-1);
numpowerunit=size(unique(costunits_offshorewind(:,1)),1);
r_unitid_lcoe = zeros(numpowerunit,1);
powers = zeros(numpowerunit,10);
off_ID = zeros(numpowerunit,4);
nnnnp = 1;
% aa1_CN=[1003, 1120, 1121, 1166, 1209, 1366, 1723, 2862, 2863, 3134, 3498]';
load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
region_ID(region_ID==8)=100;
region_ID(region_ID==9)=100;
region_ID(region_ID==10)=100;
region_ID(region_ID~=100)=0;
for coun = 1:size(bb1,1)
    idx=find(costunits_offshorewind(:,9)==bb1(coun,1));
    costunits_offshorewind1 = costunits_offshorewind(idx,:);
    [mmm1,nnn1]=find(dist_id_country120==bb1(coun,1));
    dist_id_pro120_2 = ones(max(mmm1)-min(mmm1)+1,max(nnn1)-min(nnn1)+1)*(-1);
    ID_dist_offshore120_2 = zeros(max(max(mmm1))-min(min(mmm1))+1,max(max(nnn1))-min(min(nnn1))+1);
    mm1m = min(min(mmm1));
    nnn1m = min(min(nnn1));
    dist_id_pro120_2(sub2ind(size(dist_id_pro120_2), mmm1-mm1m+1, nnn1-nnn1m+1)) = dist_id_pro120_1(sub2ind(size(dist_id_pro120_1), mmm1, nnn1));
    ID_dist_offshore120_2(sub2ind(size(ID_dist_offshore120_2), mmm1-mm1m+1, nnn1-nnn1m+1)) = ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), mmm1, nnn1));
    nnnnnn = unique(costunits_offshorewind(idx,1)); % power plant ID
    mmmmmm = unique(costunits_offshorewind(idx,8)); % county ID for other countries, pro ID for China
%     if bb1(coun,1)~=35
        [is,pos]=ismember(mmmmmm,pro_county(:,2));
        % is是与B大小一致的向量，如果在A中为1，不在为0
        % pos是B中元素如果在A中出现，出现的位置。
        pppppp = pro_county(pos,1);% pro
        off_ID(nnnnp:nnnnp+size(pppppp,1)-1,1)=nnnnnn; % power plant ID
        off_ID(nnnnp:nnnnp+size(pppppp,1)-1,2)=ones(size(pppppp,1),1)*bb1(coun,1); % country
        off_ID(nnnnp:nnnnp+size(pppppp,1)-1,3)=pppppp; % pro
        off_ID(nnnnp:nnnnp+size(pppppp,1)-1,4)=mmmmmm; % county
%     end
%     if bb1(coun,1)==35
%         pppppp = mmmmmm;% pro
%         off_ID(nnnnp:nnnnp+size(pppppp,1)-1,1)=nnnnnn; % power plant ID
%         off_ID(nnnnp:nnnnp+size(pppppp,1)-1,2)=ones(size(pppppp,1),1)*bb1(coun,1); % country
%         off_ID(nnnnp:nnnnp+size(pppppp,1)-1,3)=pppppp; % pro
%         off_ID(nnnnp:nnnnp+size(pppppp,1)-1,4)=zeros(size(pppppp,1),1); % county
%     end
    nnnnp = nnnnp+size(pppppp,1);
    
    for i=1:(size(nnnnnn,1))
        %     figure
        lcoemin=1000;
        [m,n] = find(costunits_offshorewind1(:,1)==nnnnnn(i,1));
        lcoeunit=zeros(size(m,1),8);
        %     jopt=0;
        for j=m(1):m(end)
            j1=j-m(1)+1;
            if j-m(1)>0
                lcoeunit(j1,1)=lcoeunit(j1-1,1)+costunits_offshorewind1(j,5)/10^9; % electricity TWh / year
                lcoeunit(j1,2)=lcoeunit(j1-1,2)+costunits_offshorewind1(j,4)*initialcost_ratio_country(bb1(coun,1),1); % Cost,million dollar
                lcoeunit(j1,3)=lcoeunit(j1-1,3)+costunits_offshorewind1(j,3); % CP,kW
                lcoeunit(j1,4)=lcoeunit(j1-1,4)+costunits_offshorewind1(j,6); % S km2;
                lcoeunit(j1,5)=lcoeunit(j1-1,5)-costunits_offshorewind1(j,5)/10^9*fossilfuel_emissionfactor(bb1(coun,1),1)*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
                lcoeunit(j1,9)=lcoeunit(j1-1,9)+costunits_offshorewind1(j,3)*1000*module_price/10^6*initialcost_ratio_country(bb1(coun,1),1); % module cost, million USD
                lcoeunit(j1,21)=costunits_offshorewind1(j,2);% 步长
            else
                lcoeunit(j1,1)=costunits_offshorewind1(j,5)/10^9; % electricity TWh / year
                lcoeunit(j1,2)=costunits_offshorewind1(j,4)*initialcost_ratio_country(bb1(coun,1),1); % Cost,million dollar
                lcoeunit(j1,3)=costunits_offshorewind1(j,3); % CP,kW
                lcoeunit(j1,4)=costunits_offshorewind1(j,6); % S km2;
                lcoeunit(j1,5)=-costunits_offshorewind1(j,5)/10^9*fossilfuel_emissionfactor(bb1(coun,1),1)*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
                lcoeunit(j1,9)=costunits_offshorewind1(j,3)*1000*module_price/10^6*initialcost_ratio_country(bb1(coun,1),1); % module cost, million USD
                lcoeunit(j1,21)=costunits_offshorewind1(j,2);% 步长
            end
            lcoeunit(j1,6)=lcoeunit(j1,2)*(1+OMratio_substation*discount1yr(bb1(coun,1),1)); % cost
            lcoeunit(j1,10)=lcoeunit(j1,9)*(1+OMratio_substation*discount1yr(bb1(coun,1),1)); % cost
            
            % LCoE
            lcoeunit(j1,8)=(lcoeunit(j1,6)+lcoeunit(j1,7))/lcoeunit(j1,1)/degrat1yr(bb1(coun,1),1)/1000; % LCoE million USD2019/TWh->USD2019/kWh
            if lcoeunit(j1,8)<lcoemin
                jopt=j1;
                lcoemin=lcoeunit(j1,8);
            end
        end

            if lcoeunit(jopt,3)<=10^8 % power capacity kW
                jj=jopt;
                r=1;
            end
            if lcoeunit(jopt,3)>10^8
                [Index123,~] = find(lcoeunit(:,3)<=10^8);
                if ~isempty(Index123)
                    [mm,nn]=find(lcoeunit(Index123,8)==min(min(lcoeunit(Index123,8))));
                    jj=mm;
                    r=1;
                else
                    r= 10^8/lcoeunit(1,3);
                    lcoeunit(1,1:6)=lcoeunit(1,1:6).*r;
                    lcoeunit(1,9:10)=lcoeunit(1,9:10).*r;
                    jj=1;
                end
            end
        powers(nnnnnn(i,1),1:10)=lcoeunit(jj,1:10); % 20 for LCoE USD2019/kWh
        pro_ID_minlcoe(nnnnnn(i,1),1)=lcoeunit(jj,21); % step
        
        [m,n]=find(dist_id_pro120_2==mmmmmm(i,1) & ID_dist_offshore120_2<=pro_ID_minlcoe(nnnnnn(i,1),1));
        offshorewind_power_lcoemin(sub2ind(size(offshorewind_power_lcoemin), m+mm1m-1, n+nnn1m-1))= sum(sum(powers(nnnnnn(i,1),1)));%pro_offshore120(sub2ind(size(pro_offshore120), m, n));
        unitid_lcoe(sub2ind(size(unitid_lcoe), m+mm1m-1, n+nnn1m-1))= nnnnnn(i);%pro_offshore120(sub2ind(size(pro_offshore120), m, n));
        position_off(nnnnnn(i,1),1)=mean(m+mm1m-1);
        [~,Index]=min(abs(m-mean(m)));
        n1  = n(Index);
        [~,Index1]=min(abs(n1-mean(n1)));
        position_off(nnnnnn(i,1),2)=(n1(Index1)+nnn1m-1)/4;
        r_unitid_lcoe(nnnnnn(i,1),1)= r;%pro_offshore120(sub2ind(size(pro_offshore120), m, n));
        clear jopt
    end
    clear costunits_offshorewind1
    coun
end

save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_ID_county_5%_100GW_.mat','off_ID','-v7.3'); % power plant ID,country,pro,county
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\position_off_100GW_county_5%.mat','position_off','-v7.3'); % position of power plant
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\r_unitid_lcoe_100GW_county_5%.mat','r_unitid_lcoe','-v7.3'); % 可利用面积最终确定的比例
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powers_offshorewind1227_2_xz_100GW_county_5%.mat','powers','-v7.3'); %
lcoe_offshorewind=powers(:,8);
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\lcoe_offshorewind1227_2_xz_100GW_county_5%.mat','lcoe_offshorewind','-v7.3'); % USD2019/kWh
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\unitid_lcoe_offshorewind1227_2_xz_100GW_county_5%.mat','unitid_lcoe','-v7.3'); %
azx = sum(sum(powers(:,6:7)))*10^6/(sum(sum(powers(:,3)))*10^3); %USD/W
powers(:,10)*10^6./(powers(:,3)*10^3);
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\offshorewind_power_lcoemin1227_2_xz_100GW_county_5%.mat','offshorewind_power_lcoemin','-v7.3'); % power genration of each power plant

[m,n]=find(off_ID(:,2)==35);
off_ID(m,3)
sum(powers(off_ID(m,1),1))
sum(powers(off_ID(m,1),3)/10^9)



