% 实际各国learning rate，当为nan和inf时和全球统一
tic
clear;
load('H:\global-PV-wind\ANS\3power_country_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat'); % TWh/year
load('H:\global-PV-wind\ANS\3unit_offshorewind_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat');
load('H:\global-PV-wind\ANS\learnprice_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat'); % 10 regions * 6



save('H:\global-PV-wind\ANS\3power_country_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','power_country'); % TWh/year
save('H:\global-PV-wind\ANS\3unit_offshorewind_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','unit_offshorewind_global');
save('H:\global-PV-wind\ANS\learnprice_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','learnprice_region'); % 10 regions * 6
