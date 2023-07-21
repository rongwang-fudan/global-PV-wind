tic
clear;

load('H:\Global PV and wind\ANS\index_cou_2040xzsuper.mat');
reg=0;
a = learn_country(2:11,:);
UU = [0,0;2,7;3,8;4,9;5,10;1,6];
nnn = 1;
[m,n]=find(a(1,:)==0);
if ~isempty(m)
    for i2 = 1:size(m,2)
        [m2,n2]=find(a(:,n(i2))~=0);
        if ~isempty(m2)
            a(1,n(i2)) = a(m2(1),n(i2));
        end
    end
end
for i = 2:1:10
    [m,n]=find(a(i,:)==0);
    if ~isempty(m)
        for i2 = 1:size(n,2)
            a(i,n(i2)) = a(i-1,n(i2));
        end
    end
end
for i = 1:3
    if size(unique(a(1,i*2-1:i*2)),2)>1
        xxx = unique(a(1,i*2-1:i*2));
        xxx(xxx==unique(a(2,i*2)))=[];
        a(1,i*2-1:i*2) = xxx(1);
    end
    if size(unique(a(1:5,i*2)-a(6,i*2)),1)>1
        [mv,nv]=find(a(1:5,i*2)-a(6,i*2)~=0);
        a(5,i*2-1:i*2) = a(mv(end),i*2);
    end
    
    aw = a(5:10,i*2);
    if size(unique(aw),1)>1
        v1 = unique(aw);
        for i2 = 1:size(aw,1)-1
            if size(unique(aw(i2:i2+1)),1)>1
                Inn_2070A(nnn,1:2) = aw(i2:i2+1)';
                Inn_2070A(nnn,3) = i; % type
                Inn_2070A(nnn,4) = reg; % reg
                if UU(i2+1,2)==0
                    Inn_2070A(nnn,5) = UU(i2+1,1);
                    nnn = nnn+1;
                else
                    Inn_2070A(nnn,5) = UU(i2+1,1);
                    nnn = nnn+1;
                    Inn_2070A(nnn,1:2) = aw(i2:i2+1)';
                    Inn_2070A(nnn,3) = i; % type
                    Inn_2070A(nnn,4) = reg; % reg
                    Inn_2070A(nnn,5) = UU(i2+1,2);
                    nnn = nnn+1;
                end
            else if size(unique(aw(1:i2+1)),1)>1
                    aw2 = aw(1:i2+1)-aw(i2+1);
                    [m2,n2]=find(aw2~=0);
                    Inn_2070A(nnn,1) = aw(m2(end));
                    Inn_2070A(nnn,2) = aw(i2+1);
                    Inn_2070A(nnn,3) = i; % type
                    Inn_2070A(nnn,4) = reg; % reg
                    if UU(i2+1,2)==0
                        Inn_2070A(nnn,5) = UU(i2+1,1);
                        nnn = nnn+1;
                    else
                        Inn_2070A(nnn,5) = UU(i2+1,1);
                        nnn = nnn+1;
                        Inn_2070A(nnn,1) = aw(m2(end));
                        Inn_2070A(nnn,2) = aw(i2+1);
                        Inn_2070A(nnn,3) = i; % type
                        Inn_2070A(nnn,4) = reg; % reg
                        Inn_2070A(nnn,5) = UU(i2+1,2);
                        nnn = nnn+1;
                    end
                end
            end
        end
    end
end

Inn_2070A = Inn_2070A;
nnn = 1;
Inn_2070(1,:) = Inn_2070A(1,:);
nnum(1,1) = 1;
nnn = nnn+1;
for i = 2:1:size(Inn_2070A,1)
    [m,n]=find(Inn_2070(1:nnn-1,1)==Inn_2070A(i,1) & Inn_2070(1:nnn-1,2)==Inn_2070A(i,2));
    if ~isempty(m)
        Inn_2070(m,5+(nnum(m,1)-1)*3+1:5+(nnum(m,1)-1)*3+3) = Inn_2070A(i,3:5);
        nnum(m,1) = nnum(m,1)+1;
    else
        Inn_2070(nnn,1:5) = Inn_2070A(i,:);
        nnum(nnn,1) = 1;
        nnn = nnn+1;
    end
end
save('H:\Global PV and wind\ANS\Inn_2070_2040super.mat','Inn_2070')
% type, reg, period