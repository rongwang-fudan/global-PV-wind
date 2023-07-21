tic
clear;
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2025_2040.mat')  % inv_cou,没有国家合作的投资
inv_cou2040 = inv_cou;
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2030_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2035_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2040_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2045_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2050_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2055_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2060_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2065_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];
load('H:\Global PV and wind\ANS\inv_cou_PW_noconn_2070_2040.mat')  % inv_cou
inv_cou2040 = [inv_cou2040 inv_cou];

discount = 0.03;
rr2(1,1) = 1/(1+discount)^(25-20);
rr2(2,1) = 1/(1+discount)^(30-20);
rr2(3,1) = 1/(1+discount)^(35-20);
rr2(4,1) = 1/(1+discount)^(40-20);
rr2(5,1) = 1/(1+discount)^(45-20);
rr2(6,1) = 1/(1+discount)^(50-20);
rr2(7,1) = 1/(1+discount)^(55-20);
rr2(8,1) = 1/(1+discount)^(60-20);
rr2(9,1) = 1/(1+discount)^(65-20);
rr2(10,1) = 1/(1+discount)^(70-20);
inv_cou2040_2 = inv_cou2040.*repmat(rr2',[192,1]);
inv_cou2040_22 = [zeros(size(inv_cou2040_2,1),1) inv_cou2040_2];
for i = 1:10
    inv_cou2040_dis2a(:,i) = sum(inv_cou2040_22(:,i:i+1),2)/2;
end
save('H:\Global PV and wind\ANS\inv_cou2040_dis2a_noconnect.mat','inv_cou2040_dis2a','-v7.3')  % TWh/year

%%
tic
clear;
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_8_2s_6_2a.mat')  % Inv_CCS2070_c_cou0
load('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_8_2s_6_2a.mat')  % Inv_others2070_c_cou0
discount = 0.03;
rr2(1,1) = 1/(1+discount)^(25-20);
rr2(2,1) = 1/(1+discount)^(30-20);
rr2(3,1) = 1/(1+discount)^(35-20);
rr2(4,1) = 1/(1+discount)^(40-20);
rr2(5,1) = 1/(1+discount)^(45-20);
rr2(6,1) = 1/(1+discount)^(50-20);
rr2(7,1) = 1/(1+discount)^(55-20);
rr2(8,1) = 1/(1+discount)^(60-20);
rr2(9,1) = 1/(1+discount)^(65-20);
rr2(10,1) = 1/(1+discount)^(70-20);
Inv_others2070_c_cou0_dis = Inv_others2070_c_cou0.*repmat(rr2',192,1);
Inv_others2070_c_cou0_dis2 = [zeros(size(Inv_others2070_c_cou0_dis,1),1) Inv_others2070_c_cou0_dis];
Inv_CCS2070_c_cou0_dis = Inv_CCS2070_c_cou0.*repmat(rr2',192,1);
Inv_CCS2070_c_cou0_dis2 = [zeros(size(Inv_CCS2070_c_cou0_dis,1),1) Inv_CCS2070_c_cou0_dis];
for i = 1:10
    Inv_CCS2070_c_cou0_dis2a(:,i) = sum(Inv_CCS2070_c_cou0_dis2(:,i:i+1),2)/2;
    Inv_others2070_c_cou0_dis2a(:,i) = sum(Inv_others2070_c_cou0_dis2(:,i:i+1),2)/2;
end
save('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_dis2a.mat','Inv_CCS2070_c_cou0_dis2a','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_dis2a.mat','Inv_others2070_c_cou0_dis2a','-v7.3')  % TWh/year

%%
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

rmb2us=6.8967; % RMB to USD2019
rmb2us2020=6.8996; % RMB to USD2019
rrr1 = rmb2us2020./rmb2us;
rrr = rrr1;

% load('H:\Global PV and wind\ANS\Inv_CCSall_dis2a_8_2040.mat')  %Inv_CCSall_dis2a
% load('H:\Global PV and wind\ANS\Inv_othersall_dis2a_8_2040.mat')  %Inv_othersall_dis2a
% load('H:\Global PV and wind\ANS\Inv_PWall_dis2a_8_2040.mat')  %Inv_PWall_dis2a
% load('H:\Global PV and wind\ANS\Inv_CCS0_dis2a_8_2040.mat')  %Inv_CCS0_dis2a
% load('H:\Global PV and wind\ANS\Inv_others_c0_dis2a_8_2040.mat')  %Inv_others_c0_dis2a
% aa = Inv_CCSall_dis2a+Inv_othersall_dis2a+Inv_PWall_dis2a;
% a1 = zeros(size(aa,1),1);
% for i = 1:5
%     a1(:,i+1) = aa(:,i)*2-a1(:,i);
% end
% discount = 0.03;
% rr2(1,1) = 0;
% rr2(2,1) = 1/(1+discount)^(30-20);
% rr2(3,1) = 1/(1+discount)^(40-20);
% rr2(4,1) = 1/(1+discount)^(50-20);
% rr2(5,1) = 1/(1+discount)^(60-20);
% rr2(6,1) = 1/(1+discount)^(70-20);
% a2 = a1./repmat(rr2',[size(aa,1),1]);
% a2(find(isnan(a2)==1))=0;
% % for i = 1:5
% %     Inv_dif2040(:,i) = sum(a2(:,i:i+1),2)/2;
% % end
% Inv_dif2040 = a2;

% 没有PV和wind时各国各年的CCS和other renewable的discounted to 2019之后的投资
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_dis2a.mat')  % Inv_CCS2070_c_cou0
load('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_dis2a.mat')  % Inv_others2070_c_cou0_dis2a
aa = Inv_CCS2070_c_cou0_dis2a+Inv_others2070_c_cou0_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
discount = 0.03;
rr2(1,1) = 0;
rr2(2,1) = 1/(1+discount)^(25-20);
rr2(3,1) = 1/(1+discount)^(30-20);
rr2(4,1) = 1/(1+discount)^(35-20);
rr2(5,1) = 1/(1+discount)^(40-20);
rr2(6,1) = 1/(1+discount)^(45-20);
rr2(7,1) = 1/(1+discount)^(50-20);
rr2(8,1) = 1/(1+discount)^(55-20);
rr2(9,1) = 1/(1+discount)^(60-20);
rr2(10,1) = 1/(1+discount)^(65-20);
rr2(11,1) = 1/(1+discount)^(70-20);
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv0_cou2040(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv0_cou2040 = a2;


% 没有国际合作时的CCS,other renewable和PV&wind的discounted to 2019之后的投资的变化
load('H:\Global PV and wind\ANS\Inv_CCS_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_CCS_dif_noconnect_cou_dis2a
load('H:\Global PV and wind\ANS\Inv_others_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_others_dif_noconnect_cou_dis2a
load('H:\Global PV and wind\ANS\inv_cou2040_dis2a_noconnect.mat')  % inv_cou2040_dis2a
aa = Inv_CCS_dif_noconnect_cou_dis2a+Inv_others_dif_noconnect_cou_dis2a+inv_cou2040_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv_dif2040_noconnect_cou(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv_dif2040_noconnect_cou = a2;

Inv0_coudif = Inv0_cou2040*rrr/10^6*5;
Inv_dif_noconnect_cou = Inv_dif2040_noconnect_cou*rrr/10^6*5;
NetInv = Inv0_coudif+Inv_dif_noconnect_cou;
sum(NetInv)
save('H:\Global PV and wind\ANS\NetInv_nodiscounted_5_2040peryear.mat','NetInv'); % trillion discounted net cost, trillion 2020$
