tic
clear;
load('H:\global-PV-wind\Data\powerdemand_monhour2070_ele288_1112_2_max.mat')  % 288*34 % TWh/h 
powerdemand_pro2060_ele = sum(powerdemand_monhour2070_ele288)';
clear powerdemand_monhour2070_ele288
load('H:\global-PV-wind\Data\pro_CN_reg.mat') % 1.pro ID;2.region ID in China (1-7)
for i = 1:1:7
    [m,n] = find(pro_CN_reg(:,2)==i);
    powerdemand_CN2060(i,1) = sum(powerdemand_pro2060_ele(pro_CN_reg(m,1),1));
end

load('H:\global-PV-wind\Data\ID_pro3.mat') % 1.FID; 2.FIRST_ID_0; 3.ID_country120_0214; 4.FIRST_ID_1
load('H:\global-PV-wind\Data\GDP_PPP2015.mat'); % 2011 US dollar
load('H:\global-PV-wind\Data\GADM_pro120_xz.mat')  % 0-3638
for i = 1:192
    [m,n]=find(ID_pro(:,3)==i);
    ind = unique(ID_pro(m,1)+1);
    powerdemand_country2060_ele_noothers(i,1) = sum(powerdemand_pro2060_ele(ind,1));
    clear m
    clear n
    i
end
for i = 184
    [m,n]=find(ID_pro(:,3)==i & ID_pro(:,1)~=3401);
    ind = unique(ID_pro(m,1)+1);
    powerdemand_country2060_ele_noothers(i,1) = sum(powerdemand_pro2060_ele(ind,1));
    clear m
    clear n
    i
end
[m,n]=find(ID_pro(:,1)==3401);
ind = unique(ID_pro(m,1)+1);
powerdemand_country2060_ele_noothers_Alaska = sum(powerdemand_pro2060_ele(ind,1));

load('H:\global-PV-wind\Data\Area_country.mat') % km2
Q1 = prctile(Area_country,25);
[mmm,n]=find(Area_country>Q1);
UHV_station = zeros(192,1);
UHV_station(mmm,1) = 1; 

%%
load('H:\global-PV-wind\Data\GADM_country120_xz.mat')
[m,n]=find(GADM_pro120==3401);
GADM_country120(sub2ind(size(GADM_country120), m, n))=184.1; % Alaska
mmm = [mmm;184.1];

n11 = 1;
load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
for i = 1:size(mmm,1)
    if mmm(i)~=35
        [m,n] = find(GADM_country120==mmm(i));
        GDP_1 = sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n))));
        m_mean = round(mean(m,1));
        n_mean = round(mean(n,1));
        m_mean_all(i,1) = m_mean;
        n_mean_all(i,2) = n_mean;
        for dom=1:4
            if dom==1
                idx=find(m>m_mean & n<=n_mean); %  Upleft
            elseif dom==2
                idx=find(m>m_mean & n>n_mean); %  Upright
            elseif dom==3
                idx=find(m<=m_mean & n>n_mean); %  Downright
            elseif dom==4
                idx=find(m<=m_mean & n<=n_mean); %  Downleft
            end
            if mmm(i)~=184.1
                powerdemand = powerdemand_country2060_ele_noothers(mmm(i),1)*sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx)))))/GDP_1;
            end
            if mmm(i)==184.1
                powerdemand = powerdemand_country2060_ele_noothers_Alaska*sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx)))))/GDP_1;
            end
            if sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx)))))~=0
                position_c(n11,1) = round(sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx))).*m(idx)))./sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx))))));
                position_c(n11,2) = round(sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx))).*n(idx)))./sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx))))));
            end
            if sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m(idx), n(idx)))))==0
                position_c(n11,1) = round(mean(m(idx)));
                position_c(n11,2) = round(mean(n(idx)));
            end
            position_c(n11,3) = mmm(i); % country ID
            position_c(n11,4) = region_ID(floor(position_c(n11,3)),1);% region
            position_c(n11,5) = GADM_pro120(position_c(n11,1),position_c(n11,2)); % 0-3638,pro ID
            position_c(n11,6) = powerdemand;
            position_c(n11,7) = dom;
            n11 = n11+1;
        end
        i
    end
    
    if mmm(i)==35
            position_c(n11:n11+6,3) = mmm(i); % country ID
            position_c(n11:n11+6,4) = region_ID(floor(position_c(n11,3)),1);% region
            position_c(n11:n11+6,5) = -1; % 0-3638,pro ID
            position_c(n11:n11+6,6) = powerdemand_CN2060;
            position_c(n11:n11+6,7) = [1:1:7]';
            n11 = n11+7;
        i
    end    
end

[mmm,n]=find(UHV_station==0);
for i = 1:size(mmm,1)
    [m,n] = find(GADM_country120==mmm(i));
    GDP_1 = sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n))));
    m_mean = round(mean(m,1));
    n_mean = round(mean(n,1));
    powerdemand = powerdemand_country2060_ele_noothers(mmm(i),1);
    if sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n))))~=0
        position_c(n11,1) = round(sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n)).*m))./sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n)))));
        position_c(n11,2) = round(sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n)).*n))./sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n)))));
    end
    if sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n))))==0
        position_c(n11,1) = m_mean;
        position_c(n11,2) = n_mean;
    end
    position_c(n11,3) = mmm(i); % country ID
    position_c(n11,4) = region_ID(floor(position_c(n11,3)),1);% region
    position_c(n11,5) = GADM_pro120(position_c(n11,1),position_c(n11,2)); % 0-3638,pro ID
    position_c(n11,6) = powerdemand;
    position_c(n11,7) = 0;
    n11 = n11+1;
    i
end
UHV_Station_country = [[1:1:size(position_c,1)]' position_c];
save('H:\global-PV-wind\ANS\UHV_Station_country_all.mat', 'UHV_Station_country', '-v7.3')  % SubstatIon; row；col；country ID; region ID; pro ID(0-3638); power demand (TWh/year); REG(1-4)

%% 计算UHV station之间的距离：580*580
load('H:\global-PV-wind\ANS\UHV_Station_country_all.mat')  
UHV_Station_country(:,1)=[];
position_c = UHV_Station_country;
position_c(:,1) = 90-position_c(:,1)/120+1/240; % conversion of 行 to lat
position_c(:,2) = position_c(:,2)/30-180-1/60; % conversion of 列 to lon
for i = 1:size(position_c,1)
    for j = 1:size(position_c,1)
        if abs(position_c(i,2)-position_c(j,2))<=180
            distance_UHV_Station(i,j) = distance(position_c(i,1),position_c(i,2),position_c(j,1),position_c(j,2))/180*pi*6371; % km
        end
        if abs(position_c(i,2)-position_c(j,2))>180
            if position_c(i,2)<0
                a=180+position_c(i,2);
                b=180-position_c(j,2);
                distance_UHV_Station(i,j) = distance(position_c(i,1),a,position_c(j,1),b)/180*pi*6371; % km
            end
            if position_c(i,2)>0
                a=180-position_c(i,2);
                b=180+position_c(j,2);
                distance_UHV_Station(i,j) = distance(position_c(i,1),a,position_c(j,1),b)/180*pi*6371; % km
            end
        end
    end
    i
end
save('H:\global-PV-wind\ANS\distance_UHV_Station_all.mat', 'distance_UHV_Station', '-v7.3')  % km
