tic
clear;
load('H:\Global PV and wind\ANS\3power_country_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat'); % TWh/year
load('H:\Global PV and wind\ANS\3unit_offshorewind_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat');
load('H:\Global PV and wind\ANS\learnprice_region_county0811_all_country_5%_LRglobal_inilow_off_pro2070.mat'); 



save('H:\Global PV and wind\ANS\3power_country_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','power_country'); % TWh/year
save('H:\Global PV and wind\ANS\3unit_offshorewind_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','unit_offshorewind_global');
save('H:\Global PV and wind\ANS\learnprice_region_county0811_all_country_5%_LRglobal_inilow_off_pro_Cneu.mat','learnprice_region'); 
