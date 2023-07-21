tic
clear;
rmb2us=6.8967; % RMB to USD2019
rmb2us2020=6.8996; % RMB to USD2019
rrr1 = rmb2us2020./rmb2us;
rrr = rrr1;

%% Baseline
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_2_2040baseline_cou.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_2_2040baseline_cou.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
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
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,1) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5; % trillion 2020$/10yr
end
clear Inv_PW
clear Inv_PWall_total_dis2a

%% Case A
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_2_2040CaseA_cou.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_2_2040CaseA_cou.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,2) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5;
end
clear Inv_PW
clear Inv_PWall_total_dis2a

%% Case B
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_2_2040CaseB_cou.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_2_2040CaseB_cou.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,3) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5;
end
clear Inv_PW
clear Inv_PWall_total_dis2a

%% Case C
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_5_2040.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_5_2040.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,4) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5;
end
clear Inv_PW
clear Inv_PWall_total_dis2a

%% Case D
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_6_2040xz.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_6_2040xz.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(size(aa,1),1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[size(aa,1),1]);
a2(find(isnan(a2)==1))=0;
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,5) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5;
end
clear Inv_PW
clear Inv_PWall_total_dis2a


%% Case E
load('H:\Global PV and wind\ANS\Inv_PWall_total_dis2a_7_2040.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_7_2040.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_total_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(size(aa,1),1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[size(aa,1),1]);
a2(find(isnan(a2)==1))=0;
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a
for i = 1:10
    Inv_PWaa_2040Cneutrality(i,6) = sum(sum(Inv_PW(:,i:i+1)))/2*rrr/10^6*5;
end
clear Inv_PW
clear Inv_PWall_total_dis2a

save('H:\Global PV and wind\ANS\Inv_PWaa_2040Cneutrality_2040.mat','Inv_PWaa_2040Cneutrality'); % trillion discounted net cost, trillion 2020$

Totalcost_21to70per = sum(Inv_PWaa_2040Cneutrality)/50;
Totalcost_21to30per = Inv_PWaa_2040Cneutrality(1,:)/10;


