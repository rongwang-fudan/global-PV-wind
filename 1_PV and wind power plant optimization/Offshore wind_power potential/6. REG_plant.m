tic
clear;
load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_ID_county_5%_100GW_.mat'); % off_ID % 4列分别是 power ID; Country; pro; county
load('H:\Global PV and wind\Data\position_off_100GW_county_5%.mat'); % power ID
powerunit_offshorewind(:,1)=floor(position_off(:,1)); % conversion of lat to y
powerunit_offshorewind(:,2)=floor(position_off(:,2)/4); % conversion of lon to x
clear position_off

load('H:\Global PV and wind\Data\Area_country.mat')
Q1 = prctile(Area_country,25);
[mmm,n]=find(Area_country>Q1);
load('H:\Global PV and wind\Data\GADM_country120_xz2.mat')
load('H:\Global PV and wind\Data\GADM_pro120_xz2.mat')  % 0-3638
[m,n]=find(GADM_pro120==3401);
GADM_country120(sub2ind(size(GADM_country120), m, n))=184.1; % Alaska
mmm = [mmm;184.1];

load('H:\Global PV and wind\ANS\UHV_Station_country_all.mat')
REG_plant_off = ones(size(powerunit_offshorewind,1),3)*(-1);
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
load('H:\Global PV and wind\Data\pro_CN_reg.mat') % 
for i = 1:size(mmm,1)
        [m,n] = find(GADM_country120==mmm(i));
        m_mean = round(mean(m,1));
        n_mean = round(mean(n,1));
        m_mean_all(i,1) = m_mean;
        n_mean_all(i,2) = n_mean;
        [m1,n1] = find(off_ID(:,2)==mmm(i));
        if ~isempty(m1)
            for ii = 1:size(m1,1)
                if powerunit_offshorewind(m1(ii),1)>m_mean & powerunit_offshorewind(m1(ii),2)<=n_mean
                    REG_plant_off(m1(ii),2) = 1;
                    REG_plant_off(m1(ii),1) = mmm(i); % country ID
                    [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==1);
                    REG_plant_off(m1(ii),3) = UHV_Station_country(mma,1); % UHV Station的ID
                end
                if powerunit_offshorewind(m1(ii),1)>m_mean & powerunit_offshorewind(m1(ii),2)>n_mean
                    REG_plant_off(m1(ii),2) = 2;
                    REG_plant_off(m1(ii),1) = mmm(i); % country ID
                    [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==2);
                    REG_plant_off(m1(ii),3) = UHV_Station_country(mma,1); % UHV Station的ID
                end
                if powerunit_offshorewind(m1(ii),1)<=m_mean & powerunit_offshorewind(m1(ii),2)>n_mean
                    REG_plant_off(m1(ii),2) = 3;
                    REG_plant_off(m1(ii),1) = mmm(i); % country ID
                    [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==3);
                    REG_plant_off(m1(ii),3) = UHV_Station_country(mma,1); % UHV Station的ID
                end
                if powerunit_offshorewind(m1(ii),1)<=m_mean & powerunit_offshorewind(m1(ii),2)<=n_mean
                    REG_plant_off(m1(ii),2) = 4;
                    REG_plant_off(m1(ii),1) = mmm(i); % country ID
                    [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==4);
                    REG_plant_off(m1(ii),3) = UHV_Station_country(mma,1); % UHV Station的ID
                end
            end
        end
        i
end
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\REG_plant_off_all.mat', 'REG_plant_off', '-v7.3')  % 各电厂所在的国家，REG的ID，UHV Station的ID

