tic
clear;
trans = load('H:\Global PV and wind\Data\transmissionline0626.txt');
aaa = ones(15,1);
aaa(16,1) = 2;
aaa(17:39,1) = 3;
% Column 1-15 are 1000kV alternating current power stations, Column 16 is ±1100kV direct  current power station,  Column 17-39 are ±800kV direct  current power stations 
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
load('H:\Global PV and wind\Data\region_ID_new0811.mat'); %
lines(:,11)=35;
lines(:,17)=region_ID(35);
lines(:,18)=0;
lines(1:id,1)=floor((55-lines(1:id,1))*120)+35*120+1; % conversion of lat to y
lines(1:id,3)=floor((55-lines(1:id,3))*120)+35*120+1; % conversion of lat to y
lines(1:id,2)=floor((lines(1:id,2)-73)*30)+(180+73)*30+1; % conversion of lon to x
lines(1:id,4)=floor((lines(1:id,4)-73)*30)+(180+73)*30+1; % conversion of lon to x

load('H:\Global PV and wind\ANS\position_substation_2_CN0811_all.mat'); % 2*2° position_xy  lat;lon
% 1 row； 2 col； 10 pro ID；11 country ID；17 region ID；18 power demand (TWh/year)
[m,n] = find(position_xy(:,11)==35);
position_xy(m,:)=[];

lines = [lines;position_xy];
id=size(lines,1); % id of transmission line
clear position_xy

load('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\powerunit_w_pv_county.dat','-mat');
powerunit = powerunit_w;
clear powerunit_w
load('H:\Global PV and wind\Data\GADM_pro120_xz2.mat')
load('H:\Global PV and wind\Data\GADM_county120_xz2.mat')
load('H:\Global PV and wind\Data\GADM_country120_xz2.mat')

num=size(powerunit,1);
% transmission line for substation
for i=1:num
    mindistance=1000000;
    substation=0;
    [m,n]=find(lines(1:id,11)==powerunit(i,5));
    if ~isempty(m)
        for j=1:size(m,1)
            distance=abs(powerunit(i,1)-lines(m(j),1))+abs(powerunit(i,2)-lines(m(j),2))*4;
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
            distance=abs(powerunit(i,1)-lines(j,1))+abs(powerunit(i,2)-lines(j,2))*4;
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
    lines(id+i,1:2)=powerunit(i,1:2); % start point of the line
    lines(id+i,3:4)=lines(substation,1:2); % end point of the line
    lines(id+i,5)=lines(substation,5); % major line
    lines(id+i,6)=substation;
    lines(id+i,7)=i;
    lines(id+i,8)=mindistance;
    lines(id+i,11)=powerunit(i,5);% country
    %     i
end

t=lines;
lines=t(1:(id+i),:); % delete unused data
for j=1:(id+i)
    lines(j,9)=floor(GADM_county120(lines(j,1),lines(j,2)));
    lines(j,10)=GADM_pro120(lines(j,1),lines(j,2));
end
save('H:\Global PV and wind\Code\1_PV and wind power plant optimization\PV_power potential\ANS\tranmission_lines_all.dat','lines');

