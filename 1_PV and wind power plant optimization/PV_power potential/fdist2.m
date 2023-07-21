% Determine the number of grids as a function of the distance to the center

function [ x, y, num, xy, xy12 ]  =  fdist2( gridx , gridy , gridx2 , gridy2 )

global np0 npfraction

    x=round(mean(gridx,1));
    y=round(mean(gridy,1));
    z=floor(abs(gridx-x)+abs(gridy-y)*4)+1;
    zmin=min(z,[],1); idx1=find(z==zmin); % adjustment of the center
    x=gridx(idx1(1));
    y=gridy(idx1(1));
    z=floor(abs(gridx(:)-x)+abs(gridy(:)-y)*4)+1;
    z22=floor(abs(gridx2(:)-x)+abs(gridy2(:)-y)*4)+1;

    npoint=zeros(max(z,[],1),1);
    for i=1:size(z,1)
        npoint(z(i),1)=npoint(z(i),1)+1;
    end
    for i=1:(size(npoint,1)-1)
        npoint(i+1,1)=npoint(i,1)+npoint(i+1,1);
    end
    npoint(1,1)=1;
    
%     z2=npoint(:,1)./np0(1:size(npoint,1),1);
   if size(npoint,1)<=size(np0,1)
       z2=npoint(:,1)./np0(1:size(npoint,1),1);
   else
       z2=npoint(1:size(np0,1),1)./np0(:);
   end
    idx2=find(z2>npfraction);
    
    num=npoint(idx2(end),1);
    
    idx3=find(z<=idx2(end));
    xy=zeros(size(idx3,1),2);
    xy(:,1)=gridx(idx3);
    xy(:,2)=gridy(idx3);
    
    idx32=find(z22<=idx2(end));
    xy12=zeros(size(idx32,1),2);
    xy12(:,1)=gridx2(idx32);
    xy12(:,2)=gridy2(idx32);

