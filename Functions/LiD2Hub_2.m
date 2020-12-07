%transform the position to the hub referecne system(x=downstream,
%y=spanwise,z=vertical) from the LiDAR reference system (x=W-E,y=S-N)
%04/10/2018: created
 
function Pos_Hub=LiD2Hub_2(Pos_LiD,yaw,Turbine_loc,LiDAR_loc)
%Pos_hub: (3xN) x,y,z position of the points in hub referecnce sys
%yaw: yaw angle in geophysiscal reference sys
%Turbine_loc: turbine location
%LiDAR_loc=LiDAR location
    Turbine_loc=reshape(Turbine_loc,[],1);
    LiDAR_loc=reshape(LiDAR_loc,[],1);
    angle=-90-yaw;
    M_rot=[cosd(angle),sind(angle),0;-sind(angle),cosd(angle),0;0,0,1];
    Pos_Hub=M_rot*(Pos_LiD-repmat(Turbine_loc-LiDAR_loc,1,length(Pos_LiD(1,:)))); 
end