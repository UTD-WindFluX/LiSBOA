%calculate equivalent velocity
%09/18/2018: created
%06/12/2019: finalized
function U_x=Eq_velocity(Azi,Ele,RWS,WindDirection,max_angle)
    Azi=vert(Azi);
    Ele=vert(Ele);
    Dir=[cosd(90-Azi).*cosd(Ele),sind(90-Azi).*cosd(Ele),sind(Ele)]';
    angle=wrapTo360(270-WindDirection);
    w=[cosd(angle) sind(angle) 0];%wind direction 
    mat_cos=repmat(-w*Dir,length(RWS(:,1)),1); %local scalar product between beam direction and w
    U_x=RWS./mat_cos;%approximated 1D field
    U_x(abs(mat_cos)<cosd(max_angle))=nan;
end