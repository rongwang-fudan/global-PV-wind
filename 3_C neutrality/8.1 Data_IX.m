tic
clear;
nn = 1;
for i = 1:192
    for j = 1:192
        if i ~=j
            index_ij_s(nn,1) = i;
            index_ij_s(nn,2) = j;
            nn = nn+1;
        end
    end
end

%% 2025
load('H:\Global PV and wind\ANS\r_power_coop_2025.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2025.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2025.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2025.mat')  % index_ij_IX
power_mineral_powertrans_tech2025xz = power_mineral_powertrans_tech2030*0;
Inv_pvwind_mineral_powertrans_tech2025xz = Inv_pvwind_mineral_powertrans_tech2030*0;
r_power_coop2025xz = r_power_coop2030*0;
power_mineral_powertrans_tech2025xz(1) = power_mineral_powertrans_tech2030(1);
Inv_pvwind_mineral_powertrans_tech2025xz(1) = Inv_pvwind_mineral_powertrans_tech2030(1);
r_power_coop2025xz(:,1) = r_power_coop2030(:,1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2025xz(m+1) = power_mineral_powertrans_tech2030(i+1);
    Inv_pvwind_mineral_powertrans_tech2025xz(m+1) = Inv_pvwind_mineral_powertrans_tech2030(i+1);
    r_power_coop2025xz(:,m+1) = r_power_coop2030(:,i+1);
    nn = nn+1;
    i
end
clear power_mineral_powertrans_tech2030
clear Inv_pvwind_mineral_powertrans_tech2030
clear r_power_coop2030

save('H:\Global PV and wind\ANS\r_power_coop_2025xz.mat','r_power_coop2025xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2025xz.mat','power_mineral_powertrans_tech2025xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2025xz.mat','Inv_pvwind_mineral_powertrans_tech2025xz','-v7.3')  %
clear r_power_coop2025xz

%% 2030
load('H:\Global PV and wind\ANS\r_power_coop_2030.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2030.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2030.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2030xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2030xz = power_mineral_powertrans_tech2030*0;
Inv_pvwind_mineral_powertrans_tech2030xz = Inv_pvwind_mineral_powertrans_tech2030*0;
r_power_coop2030xz = r_power_coop2030*0;
power_mineral_powertrans_tech2030xz(1) = power_mineral_powertrans_tech2030(1);
Inv_pvwind_mineral_powertrans_tech2030xz(1) = Inv_pvwind_mineral_powertrans_tech2030(1);
r_power_coop2030xz(:,1) = r_power_coop2030(:,1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2030xz(m+1) = power_mineral_powertrans_tech2030(i+1);
    Inv_pvwind_mineral_powertrans_tech2030xz(m+1) = Inv_pvwind_mineral_powertrans_tech2030(i+1);
    r_power_coop2030xz(:,m+1) = r_power_coop2030(:,i+1);
    nn = nn+1;
    i
end
clear power_mineral_powertrans_tech2030
clear Inv_pvwind_mineral_powertrans_tech2030
clear r_power_coop2030

save('H:\Global PV and wind\ANS\r_power_coop_2030xz.mat','r_power_coop2030xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2030xz.mat','power_mineral_powertrans_tech2030xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2030xz.mat','Inv_pvwind_mineral_powertrans_tech2030xz','-v7.3')  %
clear r_power_coop2030xz


%% 2035
load('H:\Global PV and wind\ANS\r_power_coop_2035.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2035.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2035.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2035xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2035xz = power_mineral_powertrans_tech2030*0;
Inv_pvwind_mineral_powertrans_tech2035xz = Inv_pvwind_mineral_powertrans_tech2030*0;
r_power_coop2035xz = r_power_coop2030*0;
power_mineral_powertrans_tech2035xz(1) = power_mineral_powertrans_tech2030(1);
Inv_pvwind_mineral_powertrans_tech2035xz(1) = Inv_pvwind_mineral_powertrans_tech2030(1);
r_power_coop2035xz(:,1) = r_power_coop2030(:,1);
nn = 1;
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2035xz(m+1) = power_mineral_powertrans_tech2030(i+1);
    Inv_pvwind_mineral_powertrans_tech2035xz(m+1) = Inv_pvwind_mineral_powertrans_tech2030(i+1);
    r_power_coop2035xz(:,m+1) = r_power_coop2030(:,i+1);
    nn = nn+1;
    i
end
clear power_mineral_powertrans_tech2030
clear Inv_pvwind_mineral_powertrans_tech2030
clear r_power_coop2030

save('H:\Global PV and wind\ANS\r_power_coop_2035xz.mat','r_power_coop2035xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2035xz.mat','power_mineral_powertrans_tech2035xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2035xz.mat','Inv_pvwind_mineral_powertrans_tech2035xz','-v7.3')  %
clear r_power_coop2035xz

%% 2040
load('H:\Global PV and wind\ANS\r_power_coop_2040.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2040.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2040.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2040xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2040xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2040xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2040xz = r_power_coop2040*0;
power_mineral_powertrans_tech2040xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2040xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2040xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2040xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2040xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2040xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2040xz.mat','r_power_coop2040xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2040xz.mat','power_mineral_powertrans_tech2040xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2040xz.mat','Inv_pvwind_mineral_powertrans_tech2040xz','-v7.3')  %
clear r_power_coop2040xz

%% 2045
load('H:\Global PV and wind\ANS\r_power_coop_2045.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2045.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2045.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2045xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2045xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2045xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2045xz = r_power_coop2040*0;
power_mineral_powertrans_tech2045xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2045xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2045xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2045xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2045xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2045xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2045xz.mat','r_power_coop2045xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2045xz.mat','power_mineral_powertrans_tech2045xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2045xz.mat','Inv_pvwind_mineral_powertrans_tech2045xz','-v7.3')  %
clear r_power_coop2045xz

%% 2050
load('H:\Global PV and wind\ANS\r_power_coop_2050.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2050.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2050.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2050xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2050xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2050xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2050xz = r_power_coop2040*0;
power_mineral_powertrans_tech2050xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2050xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2050xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2050xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2050xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2050xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2050xz.mat','r_power_coop2050xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2050xz.mat','power_mineral_powertrans_tech2050xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2050xz.mat','Inv_pvwind_mineral_powertrans_tech2050xz','-v7.3')  %
clear r_power_coop2050xz

%% 2055
load('H:\Global PV and wind\ANS\r_power_coop_2055.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2055.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2055.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2055xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2055xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2055xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2055xz = r_power_coop2040*0;
power_mineral_powertrans_tech2055xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2055xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2055xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2055xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2055xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2055xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2055xz.mat','r_power_coop2055xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2055xz.mat','power_mineral_powertrans_tech2055xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2055xz.mat','Inv_pvwind_mineral_powertrans_tech2055xz','-v7.3')  %
clear r_power_coop2055xz

%% 2060
load('H:\Global PV and wind\ANS\r_power_coop_2060.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2060.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2060.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2060xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2060xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2060xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2060xz = r_power_coop2040*0;
power_mineral_powertrans_tech2060xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2060xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2060xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2060xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2060xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2060xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2060xz.mat','r_power_coop2060xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2060xz.mat','power_mineral_powertrans_tech2060xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2060xz.mat','Inv_pvwind_mineral_powertrans_tech2060xz','-v7.3')  %
clear r_power_coop2060xz

%% 2065
load('H:\Global PV and wind\ANS\r_power_coop_2065.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2065.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2065.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2065xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2065xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2065xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2065xz = r_power_coop2040*0;
power_mineral_powertrans_tech2065xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2065xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2065xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2065xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2065xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2065xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2065xz.mat','r_power_coop2065xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2065xz.mat','power_mineral_powertrans_tech2065xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2065xz.mat','Inv_pvwind_mineral_powertrans_tech2065xz','-v7.3')  %
clear r_power_coop2065xz

%% 2070
load('H:\Global PV and wind\ANS\r_power_coop_2070.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2070xz.mat')  % index_ij_IX
power_mineral_powertrans_tech2070xz = power_mineral_powertrans_tech2040*0;
Inv_pvwind_mineral_powertrans_tech2070xz = Inv_pvwind_mineral_powertrans_tech2040*0;
r_power_coop2070xz = r_power_coop2040*0;
power_mineral_powertrans_tech2070xz(1) = power_mineral_powertrans_tech2040(1);
Inv_pvwind_mineral_powertrans_tech2070xz(1) = Inv_pvwind_mineral_powertrans_tech2040(1);
r_power_coop2070xz(:,1) = r_power_coop2040(:,1);
for i = 1:size(index_ij_IX,1)
    [m,n]=find(index_ij_s(:,1)==index_ij_IX(i,1) & index_ij_s(:,2)==index_ij_IX(i,2));
    power_mineral_powertrans_tech2070xz(m+1) = power_mineral_powertrans_tech2040(i+1);
    Inv_pvwind_mineral_powertrans_tech2070xz(m+1) = Inv_pvwind_mineral_powertrans_tech2040(i+1);
    r_power_coop2070xz(:,m+1) = r_power_coop2040(:,i+1);
    i
end
clear power_mineral_powertrans_tech2040
clear Inv_pvwind_mineral_powertrans_tech2040
clear r_power_coop2040

save('H:\Global PV and wind\ANS\r_power_coop_2070xz.mat','r_power_coop2070xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070xz.mat','power_mineral_powertrans_tech2070xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070xz.mat','Inv_pvwind_mineral_powertrans_tech2070xz','-v7.3')  %
clear r_power_coop2070xz
