tic
clear;
load('H:\Global PV and wind\ANS\B_UHV_STO_INT_county_all_withUHVcost_pro2_8_2040_2s_2070_test6xz.mat')  % B_utilize_trans_storage
B_utilize_trans_storage5(:,5) = B_utilize_trans_storage;
clear B_utilize_trans_storage
Plant_ID = [1:1:size(B_utilize_trans_storage5,1)]';
[B,IX] = sort(B_utilize_trans_storage5(:,5));
B_utilize_trans_storage_IX = B_utilize_trans_storage5(IX,:);
Plant_ID_IX = Plant_ID(IX,:,:);
clear B_utilize_trans_storage5
clear Plant_ID
[m,n]=find(B_utilize_trans_storage_IX(:,5)==Inf);
B_utilize_trans_storage_IX(m,:) = [];
Plant_ID_IX(m) = [];
clear B_utilize_trans_storage_IX
% 
load('H:\Global PV and wind\ANS\Ph_PVwind2040_8_2s_6_2a.mat')  % Ph_PVwind5
Ph_PVwind2040 = Ph_PVwind5(2:end,2);
clear Ph_PVwind5
load('H:\Global PV and wind\ANS\m_2040_8_2s_6_2a.mat')  % m
m_cho = m;
Plant_ID_IX(m_cho,:) = []; 
Ph_PVwind2040(m_cho,:) = [];
save('H:\Global PV and wind\ANS\Plant_ID_IX_2040_8.mat', 'Plant_ID_IX', '-v7.3')
[m,n]=find(Ph_PVwind2040==0);
Plant_ID_IX(m,:) = []; 
save('H:\Global PV and wind\ANS\Plant_ID_IX_2040_22_8.mat', 'Plant_ID_IX', '-v7.3')
