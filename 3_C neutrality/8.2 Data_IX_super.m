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
load('H:\Global PV and wind\ANS\r_power_coop_2025super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2025super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2025super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2025super.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2025xzsuper.mat','r_power_coop2025xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2025xzsuper.mat','power_mineral_powertrans_tech2025xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2025xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2025xz','-v7.3')  %
clear r_power_coop2025xz

%% 2030
load('H:\Global PV and wind\ANS\r_power_coop_2030super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2030super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2030super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2030xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2030xzsuper.mat','r_power_coop2030xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2030xzsuper.mat','power_mineral_powertrans_tech2030xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2030xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2030xz','-v7.3')  %
clear r_power_coop2030xz


%% 2035
load('H:\Global PV and wind\ANS\r_power_coop_2035super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2035super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2035super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2035xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2035xzsuper.mat','r_power_coop2035xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2035xzsuper.mat','power_mineral_powertrans_tech2035xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2035xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2035xz','-v7.3')  %
clear r_power_coop2035xz

%% 2040
load('H:\Global PV and wind\ANS\r_power_coop_2040super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2040super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2040super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2040xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2040xzsuper.mat','r_power_coop2040xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2040xzsuper.mat','power_mineral_powertrans_tech2040xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2040xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2040xz','-v7.3')  %
clear r_power_coop2040xz

%% 2045
load('H:\Global PV and wind\ANS\r_power_coop_2045super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2045super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2045super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2045xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2045xzsuper.mat','r_power_coop2045xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2045xzsuper.mat','power_mineral_powertrans_tech2045xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2045xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2045xz','-v7.3')  %
clear r_power_coop2045xz

%% 2050
load('H:\Global PV and wind\ANS\r_power_coop_2050super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2050super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2050super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2050xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2050xzsuper.mat','r_power_coop2050xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2050xzsuper.mat','power_mineral_powertrans_tech2050xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2050xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2050xz','-v7.3')  %
clear r_power_coop2050xz

%% 2055
load('H:\Global PV and wind\ANS\r_power_coop_2055super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2055super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2055super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2055xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2055xzsuper.mat','r_power_coop2055xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2055xzsuper.mat','power_mineral_powertrans_tech2055xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2055xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2055xz','-v7.3')  %
clear r_power_coop2055xz

%% 2060
load('H:\Global PV and wind\ANS\r_power_coop_2060super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2060super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2060super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2060xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2060xzsuper.mat','r_power_coop2060xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2060xzsuper.mat','power_mineral_powertrans_tech2060xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2060xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2060xz','-v7.3')  %
clear r_power_coop2060xz

%% 2065
load('H:\Global PV and wind\ANS\r_power_coop_2065super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2065super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2065super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2065xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2065xzsuper.mat','r_power_coop2065xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2065xzsuper.mat','power_mineral_powertrans_tech2065xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2065xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2065xz','-v7.3')  %
clear r_power_coop2065xz

%% 2070
load('H:\Global PV and wind\ANS\r_power_coop_2070super.mat')  % r_power_coop2030
% plant*connection(0-max)
load('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070super.mat')  % power_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070super.mat')  % Inv_pvwind_mineral_powertrans_tech2030
load('H:\Global PV and wind\ANS\index_ij_IX_2070xzsuper.mat')  % index_ij_IX
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

save('H:\Global PV and wind\ANS\r_power_coop_2070xzsuper.mat','r_power_coop2070xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\power_mineral_powertrans_tech_2070xzsuper.mat','power_mineral_powertrans_tech2070xz','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_pvwind_mineral_powertrans_tech_2070xzsuper.mat','Inv_pvwind_mineral_powertrans_tech2070xz','-v7.3')  %
clear r_power_coop2070xz
