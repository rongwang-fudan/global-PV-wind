tic
clear;

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_w_pv_county.dat','-mat');
powerunit = powerunit_w;
clear powerunit_w

pos = zeros(21600,10800);
for i = 1:size(powerunit,1)
    pos(powerunit(i,1), powerunit(i,2))=i;
end

load('H:\Global PV and wind\Data\Area_country.mat')
Q1 = prctile(Area_country,25);
[mmm,n]=find(Area_country>Q1);
load('H:\Global PV and wind\Data\GADM_country120_xz2.mat')
load('H:\Global PV and wind\Data\GADM_pro120_xz2.mat')  % 0-3638
[m,n]=find(GADM_pro120==3401);
GADM_country120(sub2ind(size(GADM_country120), m, n))=184.1; % Alaska
mmm = [mmm;184.1];
Country_ID_withUHV = mmm;
save('H:\Global PV and wind\ANS\Country_ID_withUHV.mat', 'Country_ID_withUHV', '-v7.3')  % REG ID of power plant

load('H:\Global PV and wind\ANS\UHV_Station_country_all.mat')
% 1SubstatIon; 2 row；3 col；4 country ID; 5 region ID; 6 pro ID(0-3638); 7
% power demand
% demand (TWh/year); 8REG(1-4)
n11 = 1;
n12 = 1;
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
load('H:\Global PV and wind\Data\pro_CN_reg.mat') % 1. pro ID， 2. region ID in China (1-7)
for i = 1:size(mmm,1)
    if mmm(i)~=35
        [m,n] = find(GADM_country120==mmm(i));
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
            aa=unique(pos(sub2ind(size(pos), m(idx), n(idx))));
            aa(1)=[];
            if ~isempty(aa)
                Num_pant_PV(n12:n12+size(aa,1)-1,1) = dom;% UHV_Station_country id
                Num_pant_PV(n12:n12+size(aa,1)-1,2) = aa;% plant id
                Num_pant_PV(n12:n12+size(aa,1)-1,3) = size(aa,1); % power plant number
                Num_pant_PV(n12:n12+size(aa,1)-1,4) = mmm(i); % country ID
                [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==dom);
                Num_pant_PV(n12:n12+size(aa,1)-1,5) = UHV_Station_country(mma,1); % UHV Station的ID
                n11 = n11+1;
                n12 = n12+size(aa,1);
            end
            if isempty(aa)
                Num_pant_PV(n12,1) = dom;% UHV_Station_country id
                Num_pant_PV(n12,2) = nan;% plant id
                Num_pant_PV(n12,3) = 0; % power plant number
                Num_pant_PV(n12,4) = mmm(i); % country ID
                n11 = n11+1;
                n12 = n12+1;
            end
        end
        i
    end
    
    if mmm(i)==35
        for dom=1:1:7
            [mp,np] = find(pro_CN_reg(:,2)==dom);
            aa1 = [];
            for iii = 1:size(mp,1)
                [mmm22,nnn22] = find(GADM_pro120==pro_CN_reg(mp(iii),1));
                if ~isempty(mmm22)
                    aa=unique(pos(sub2ind(size(pos), mmm22, nnn22)));
                    aa(1)=[];
                    aa1 = [aa1;aa];
                end
            end
            aa = aa1;
            clear aa1
            if ~isempty(aa)
                Num_pant_PV(n12:n12+size(aa,1)-1,1) = dom;% UHV_Station_country id
                Num_pant_PV(n12:n12+size(aa,1)-1,2) = aa;% plant id
                Num_pant_PV(n12:n12+size(aa,1)-1,3) = size(aa,1); % power plant ID
                Num_pant_PV(n12:n12+size(aa,1)-1,4) = mmm(i); % country ID
                [mma,nna]=find(UHV_Station_country(:,4)==mmm(i) & UHV_Station_country(:,8)==dom);
                Num_pant_PV(n12:n12+size(aa,1)-1,5) = UHV_Station_country(mma,1); % UHV Station ID
                n11 = n11+1;
                n12 = n12+size(aa,1);
            end
            if isempty(aa)
                Num_pant_PV(n12,1) = dom;% UHV_Station_country id
                Num_pant_PV(n12,2) = nan;% plant id
                Num_pant_PV(n12,3) = 0; % plant number
                Num_pant_PV(n12,4) = mmm(i); % country ID
                n11 = n11+1;
                n12 = n12+1;
            end
        end
        i
    end
end

aa=unique(Num_pant_PV(:,2));
aa(find(isnan(aa)==1))=[];

powerdemand_plant_pv = ones(size(powerunit,1),1)*(-1);
REG_plant_pv = ones(size(powerunit,1),3)*(-1);
for i = 1:size(aa,1)
    [m,n]=find(Num_pant_PV(:,2)==aa(i));
    powerdemand_plant_pv(Num_pant_PV(m,2),1) = UHV_Station_country(Num_pant_PV(m,5),7)./Num_pant_PV(m,3);
    REG_plant_pv(Num_pant_PV(m,2),1) = Num_pant_PV(m,4); % country ID
    REG_plant_pv(Num_pant_PV(m,2),2) = Num_pant_PV(m,1); % ,1-4，China 1-7
    REG_plant_pv(Num_pant_PV(m,2),3) = Num_pant_PV(m,5); % UHV Station ID
end
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerdemand_plant_pv_all.mat', 'powerdemand_plant_pv', '-v7.3')  % power demand (TWh/year)
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\REG_plant_pv_all.mat', 'REG_plant_pv', '-v7.3')  % country，REG ID， UHV Station ID