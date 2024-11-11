% Author: Rong Wang
% Date: 2021.5.16

tic
clear;
trans = load('H:\global-PV-wind\Data\transmissionline0626.txt');
aaa = ones(15,1);
aaa(16,1) = 2;
aaa(17:39,1) = 3;
% 1-15行是1000kV交流电站，对应变电站； 1
% 16行是±1100kV直流电站，对应换流站； 2
% 17-39行±800kV直流电站，对应换流站； 3
% save transmissionline0626.txt -ascii trans
% image(unitid,'cdatamapping','scaled');
id=0; % id of transmission line
for i=1:size(trans,1)
    idx=find(trans(i,:)~=0);
    num=(size(idx,2)-1)/2;
    if num>2
        % first part
        id=id+1;
        lines(id,1:2)=trans(i,2:3);
        lines(id,3:4)=trans(i,6:7);
        lines(id,5)=i;
        lines(id,16)=aaa(i,1);
        % centeral part
        for j=1:(num-3)
            id=id+1;
            lines(id,1:4)=trans(i,(j*2+4):(j*2+7));
            lines(id,5)=i;
            lines(id,16)=aaa(i,1);
        end       
        % last part        
        id=id+1;
        lines(id,1:2)=trans(i,idx(end-1):idx(end));
        lines(id,3:4)=trans(i,4:5);
        lines(id,5)=i;
        lines(id,16)=aaa(i,1);
    else
        id=id+1;
        lines(id,1:4)=trans(i,2:5);
        lines(id,5)=i;
        lines(id,16)=aaa(i,1);
    end
end
load('H:\global-PV-wind\Data\region_ID_new0811.mat'); %
lines(:,11)=35;
lines(:,17)=region_ID(35);
lines(:,18)=0;
lines(1:id,1)=floor((90-lines(1:id,1))*120)+1; % conversion of lat to y
lines(1:id,2)=floor((lines(1:id,2)+180)*30)+1; % conversion of lon to x
lines(1:id,3)=floor((90-lines(1:id,3))*120)+1; % conversion of lat to y
lines(1:id,4)=floor((lines(1:id,4)+180)*30)+1; % conversion of lon to x

load('H:\global-PV-wind\ANS\position_substation_2_CN0811_all.mat'); % 2*2° position_xy  lat;lon
% 1 行； 2 列； 10 pro ID；11 国家ID；17 region ID；18 该序号所分配的power demand (TWh/year)
[m,n] = find(position_xy(:,11)==35);
position_xy(m,:)=[];
lines = [lines;position_xy];
id=size(lines,1); % id of transmission line

% N55-N15, E73-E138 lat1/120; lon1/30
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\position_off_100GW_county_5%.mat'); % power ID
powerunit_offshorewind(:,1)=floor(position_off(:,1)); % conversion of lat to y
powerunit_offshorewind(:,2)=floor(position_off(:,2)); % conversion of lon to x
clear position_off
load('H:\global-PV-wind\Data\GADM_pro120_xz2.mat')

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_ID_county_5%_100GW_.mat'); % off_ID % 4列分别是 power ID; Country; pro; county
num=size(powerunit_offshorewind,1);
for i=1:num
    mindistance=1000000;
    substation=0;
    %     coun = powerunit(i,5);
    [m,n]=find(lines(1:id,11)==off_ID(i,2));
    if ~isempty(m)
        for j=1:size(m,1)
            distance=abs(powerunit_offshorewind(i,1)-lines(m(j),1))+abs(powerunit_offshorewind(i,2)-lines(m(j),2))*4;
            if distance<mindistance % & lines(m(j),11)== powerunit(i,5)
                if lines(m(j),3)~=0
                    prov=GADM_pro120(lines(m(j),3),lines(m(j),4));
                    if prov~=-1 % >=1 && prov<=34
                        mindistance=distance;
                        substation=m(j);
                    end
                end
                if lines(m(j),3)==0
                    mindistance=distance;
                    substation=m(j);
                end
                
            end
        end
    else
        for j=1:id
            distance=abs(powerunit_offshorewind(i,1)-lines(j,1))+abs(powerunit_offshorewind(i,2)-lines(j,2))*4;
            if distance<mindistance % & lines(m(j),11)== powerunit(i,5)
                if lines(j,3)~=0
                    prov=GADM_pro120(lines(j,3),lines(j,4));
                    if prov~=-1 % >=1 && prov<=34
                        mindistance=distance;
                        substation=j;
                    end
                end
                if lines(j,3)==0
                    mindistance=distance;
                    substation=j;
                end
            end
        end
    end
    lines(id+i,1:2)=powerunit_offshorewind(i,1:2); % start point of the line
    lines(id+i,3:4)=lines(substation,1:2); % end point of the line
    lines(id+i,5)=lines(substation,5); % major line
    lines(id+i,6)=substation;
    lines(id+i,7)=i;
    lines(id+i,8)=mindistance;
    lines(id+i,11)=off_ID(i,2);% country
    %     i
end


t=lines;
lines=t(1:(id+i),:); % delete unused data
load('H:\global-PV-wind\Data\GADM_county120_xz2.mat')
% for j=1:(id+i)
%     lines(j,9)=floor(GADM_county120(lines(j,1),lines(j,2)));
%     lines(j,10)=GADM_pro120(lines(j,1),lines(j,2));
% end
lines(id+1:end,9) = off_ID(:,4); % county
lines(id+1:end,10) = off_ID(:,3); % pro
lines(id+1:end,11) = off_ID(:,2); % country
                                                                                                                                                                                                                  
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_100GW_county_5%.dat','lines');

% clear
% toc

%%
tic
clear
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powers_offshorewind1227_2_xz_100GW_county_5%.mat'); % USD2019/kWh
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\position_off_100GW_county_5%.mat'); % power ID
load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_ID_county_5%_100GW_.mat'); % power ID

% 运行5.2 transmission_offshorewind

load('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_100GW_county_5%.dat','-mat');
powerunit(:,1)=floor(position_off(:,1)); % conversion of lat to y
powerunit(:,2)=floor(position_off(:,2)); % conversion of lon to x
clear position_off
numpowerunit=size(powerunit,1);

[B,IX]=sort(powers(:,8),1); %IX时powers(:,20)最小值到最大值的序号
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
powerunit_IX_offshorewind = IX;
powerunit_num_IX_offshorewind = powerunit(IX,:);
lines_IX = lines(numlines+IX,:);
optpowerunit_offshorewind = powers(IX,:);
optpowerunit_offshorewind(:,31) = powerunit_IX_offshorewind;
off_pro_IX = off_ID(IX,:); % 4列分别是 power ID; Country; pro; county

save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\optpowerunit_offshorewind_100GW_county_5%.mat','optpowerunit_offshorewind'); % power ID
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_IX_offshorewind_100GW_county_5%.mat','powerunit_IX_offshorewind'); % power ID
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\tranmission_lines_IX_100GW_county_5%.mat','lines_IX'); %
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\powerunit_num_IX_offshorewind_100GW_county_5%.mat','powerunit_num_IX_offshorewind'); %
save('H:\global-PV-wind\Code\1_PV and wind power plant optimization\Offshore wind_power potential\ANS\off_pro_IX_100GW_county_5%.mat','off_pro_IX'); %


