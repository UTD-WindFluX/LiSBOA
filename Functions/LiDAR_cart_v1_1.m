%gives the location in a cartesian grid centered in the LiDAR of the exp
%points.
%11/21/2018:created
%08/12/2019: added LiDAR location

function [X,Y,Z]=LiDAR_cart_v1_1(LiDAR,Position)
    dR=nanmean(diff(round(LiDAR.Range)));
    dAZI=round(nanmean(diff(LiDAR.Azi)),2);
    dELE=round(nanmean(diff(LiDAR.Ele)),2);

    switch Position
        case 'o' 
            R=vert(LiDAR.Range);
            AZI=hor(LiDAR.Azi);
            ELE=hor(LiDAR.Ele);
        case 'SE'
            switch LiDAR.Scan_type
                case 'PPI'
                    R=vert(LiDAR.Range)-dR/2;
                    AZI=hor(LiDAR.Azi)+dAZI/2;
                    ELE=hor(LiDAR.Ele);
                case 'RHI'
                    R=vert(LiDAR.Range)+dR/2;
                    AZI=hor(LiDAR.Azi);
                    ELE=hor(LiDAR.Ele)-dELE/2;
            end
        case 'NE'
            switch LiDAR.Scan_type
                case 'PPI'
                    R=vert(LiDAR.Range)+dR/2;
                    AZI=hor(LiDAR.Azi)+dAZI/2;
                    ELE=hor(LiDAR.Ele);
                case 'RHI'
                    R=vert(LiDAR.Range)+dR/2;
                    AZI=hor(LiDAR.Azi);
                    ELE=hor(LiDAR.Ele)+dELE/2;
            end
        case 'NW'
            switch LiDAR.Scan_type
                case 'PPI'
                    R=vert(LiDAR.Range)+dR/2;
                    AZI=hor(LiDAR.Azi)-dAZI/2;
                    ELE=hor(LiDAR.Ele);
                case 'RHI'
                    R=vert(LiDAR.Range)-dR/2;
                    AZI=hor(LiDAR.Azi);
                    ELE=hor(LiDAR.Ele)+dELE/2;
            end
        case 'SW'
            switch LiDAR.Scan_type
                case 'PPI'
                    R=vert(LiDAR.Range)-dR/2;
                    AZI=hor(LiDAR.Azi)-dAZI/2;
                    ELE=hor(LiDAR.Ele);
                case 'RHI'
                    R=vert(LiDAR.Range)-dR/2;
                    AZI=hor(LiDAR.Azi);
                    ELE=hor(LiDAR.Ele)-dELE/2;
            end
    end
    X=R*(cosd(90-AZI).*cosd(ELE))+LiDAR.loc(1);
    Y=R*(sind(90-AZI).*cosd(ELE))+LiDAR.loc(2);
    Z=R*sind(ELE)+LiDAR.loc(3);
end