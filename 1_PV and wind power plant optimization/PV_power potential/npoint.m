
% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.5.21
% Determine the number of grids as a function of the distance to the center

function np0  =  npoint( gridx , gridy )

    x=round(mean(gridx,1));
    y=round(mean(gridy,1));
    z=floor(abs(gridx-x)+abs(gridy-y)*4)+1;
    zmin=min(z,[],1); idx1=find(z==zmin); % adjustment of the center
    x=gridx(idx1(1));
    y=gridy(idx1(1));
    z=floor(abs(gridx(:)-x)+abs(gridy(:)-y)*4)+1;
    np0=zeros(max(z,[],1),1);
    for i=1:size(z,1)
        np0(z(i),1)=np0(z(i),1)+1;
    end
    for i=1:(size(np0,1)-1)
        np0(i+1,1)=np0(i,1)+np0(i+1,1);
    end
    np0(1,1)=1;
    
    
