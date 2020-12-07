%Example of LiSBOA of LiDAR data
%12/04/2020: created, finalized
restoredefaultpath
clear all
close all
addpath('./Functions');

%% Inputs
source_met_data='./Data/0813_1543_MET_SCADA.mat';%source of met and SCADA data
source_LiDAR_data='./Data/Filt_737345.1042_737345.1459_154_154_154_600_0.001_-18_0.5.mat';%source of filtered LiDAR data

DWD=10;%[°] maximum range of wind direction
LiDAR_loc=[912.3    9708.4   -37.1];%[m] LiDAR location

D=77;%[m] turbine diamter GE 1.5sle

%flow characteristics
T=1200;%[s] total sampling time
tau=5;%[s] characteristic integral time scale
u_var=1.5;%[m/s] charactersitic velocity standard deviation

%origin of the Cartesian coordinates [m]
x0=LiDAR_loc(1);
y0=LiDAR_loc(2);
z0=LiDAR_loc(3);

%Cartesian domain size (in the domain aligned with wind direction and centered in [x0,y0,z0])
xmin=300;%[m]
xmax=1300;%[m];
ymin=-500;%[m]
ymax=450;%[m]
zmax=200;%[m]
zmin=-30;%[m]

%LiSBOA parameters
sigma=1/4;%smoothing parameter
m=5;%number of iterations
order=[2 3 4];%order of the central statistical moments
Dn_x=2.5*D;%[m] fundamental half-wavelength in x (streamwise)
Dn_y=0.5*D;%[m] fundamental half-wavelength in y (spanwise)
Dn_z=0.5*D;%[m] fundamental half-wavelength in z (vertical)
grid_factor=0.25;%ratio of the resolution of the Cartesian to the fundamental half-wavelength vector
R_max=3;%[sigmas] radius of influence
max_Dd=1;%maximum local data spacing (in the scaled frame of reference)
tol_dist=1/14;%tolerance of the spatial location of samples in the scaled frame of reference, points closer than this are assumed as unique
max_iter=5;%maximu number of iterations

%% Initilization
load(source_met_data);
load(source_LiDAR_data);
Orientation=wrapTo360(180+WD_avg);%[°] wind direction (0°=North,90°=East,180°=South,270°=West); the x axis will be aligned with the stremwise direction

%zeroing
X_all=[];
Y_all=[];
Z_all=[];
Ueq_all=[];
inc=0;

%% LiDAR data preconditioning
for ID_time=1:length(LiDAR_all)
    
    LiDAR=LiDAR_all{ID_time};
    
    if min(size(LiDAR.RWS))>1 && abs(angdiffd(interp1(t_sel,WD,nanmean(LiDAR.Time),'pchip'),WD_avg))<=DWD/2
        if strcmp(LiDAR.Scan_type,'PPI') || strcmp(LiDAR.Scan_type,'RHI')
            inc=inc+1;
            LiDAR.RWS=-LiDAR.RWS;%reverse radial velocity (for Halo LiDAR)
            LiDAR.RWS_filt=-LiDAR.RWS_filt;
            LiDAR.RWS_filt(1:6,:)=nan;%excluding first gates
            
            %SCADA/MET
            t_LiDAR(inc)=nanmean(LiDAR.Time);
            U_freestream(inc)=interp1(t_sel,WS,t_LiDAR(inc),'pchip');
            WindDirection(inc)=interp1(t_sel,WD,t_LiDAR(inc),'pchip');
            
            %position
            [LiDAR.X,LiDAR.Y,LiDAR.Z]=LiDAR_cart_v1_1(LiDAR,'o');
            
            %eq. X (for RHI only)
            LiDAR.X_eq=LiDAR.Range.*cosd(LiDAR.Ele)'*sign(cosd(90-nanmean(LiDAR.Azi))');
            
            %ref. system change
            Pos_LiD=[reshape(LiDAR.X,1,[]);reshape(LiDAR.Y,1,[]);reshape(LiDAR.Z,1,[])];
            Pos_Hub=LiD2Hub_2(Pos_LiD,Orientation,[x0;y0;z0],[0;0;0]);
            
            %interpolation of isolated nans
            L=~isnan(LiDAR.RWS_filt);
            out=Find_bound3D(L);
            LiDAR.RWS_inpainted=inpaint_nans(LiDAR.RWS_filt,0);
            LiDAR.RWS_inpainted(out)=nan;
            
            %equivalent velocity
            Ueq=Eq_velocity(LiDAR.Azi,LiDAR.Ele,LiDAR.RWS_inpainted,WindDirection(inc),60)/U_freestream(inc);
            
            %stacking
            X_all=[X_all;Pos_Hub(1,:)'];
            Y_all=[Y_all;Pos_Hub(2,:)'];
            Z_all=[Z_all;Pos_Hub(3,:)'];
            Ueq_all=[Ueq_all;reshape(Ueq,[],1)];
            
            disp(datestr(nanmean(LiDAR.Time)));
        end
    end
end

%3D reconstruction
[X2,Y2,Z2,excl,Dd,Avg_all,HOM_all]=LiSBOA(X_all,Y_all,Z_all,xmin,xmax,ymin,ymax,zmin,zmax,Dn_x,Dn_y,Dn_z,...
    sigma,R_max,grid_factor,tol_dist,max_Dd,'all',Ueq_all,order,max_iter);

%Petersen-Middleton enforced
U_mean=Avg_all{m+1};
U_mean(excl)=nan;
for ID_order=1:length(order)
    HOM{ID_order}=HOM_all{ID_order}{m+1};
    HOM{ID_order}(excl)=nan;
end

%% Plot
close all
mkfig('max');
scatter3(X_all,Y_all,Z_all,4,colorfcn_v2_1(Ueq_all,0.4,1.1,'coolwarm'),'filled');axis equal;
xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');
axis([xmin xmax ymin ymax zmin zmax]);
smart_colorbar_v2_1(0.4:0.1:1.1,array2cell(0.4:0.1:1.1),'coolwarm','$u_{eq}$');
title('Non-dimensional LiDAR data aigned with wind direction');TNR(13);
view([325 43]);

%permutation
U_plot=permute(U_mean,[3 1 2]);
X2_plot=permute(X2,[3 1 2]);
Y2_plot=permute(Y2,[3 1 2]);
Z2_plot=permute(Z2,[3 1 2]);
for ID_order=1:length(order)
    HOM_plot{ID_order}=permute(HOM{ID_order},[3 1 2]);
end

%plot mean flow
mkfig('max');
for ID_X=1:length(X2_plot(1,1,:))
    Pos_Hub=[hor(X2_plot(:,:,ID_X));hor(Y2_plot(:,:,ID_X));hor(Z2_plot(:,:,ID_X))];
    Pos_LiD=Hub2LiD_2(Pos_Hub,Orientation,[x0;y0;z0],[0;0;0]);
    X2_LiD=reshape(Pos_LiD(1,:),size(U_plot(:,:,1)));
    Y2_LiD=reshape(Pos_LiD(2,:),size(U_plot(:,:,1)));
    Z2_LiD=reshape(Pos_LiD(3,:),size(U_plot(:,:,1)));
    pcolor3(X2_LiD,Y2_LiD,Z2_LiD,U_plot(:,:,ID_X));
    hold on;shading flat;alpha 0.5;colormap coolwarm;
end
   view([244 47]);grid on;ax0.View=[30 30];
xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');c=colorbar;c.Label.String='$\overline{u}$';caxis([0.4 1.1]);
axis equal;plot3(x0,y0,z0,'.b','markersize',15);plot3(LiDAR_loc(1),LiDAR_loc(2),LiDAR_loc(3),'xk','markersize',15,'linewidth',2);
TNR(13);

%turbulence intensity
try
    mkfig('max');
    for ID_X=1:length(X2_plot(1,1,:))
        Pos_Hub=[hor(X2_plot(:,:,ID_X));hor(Y2_plot(:,:,ID_X));hor(Z2_plot(:,:,ID_X))];
        Pos_LiD=Hub2LiD_2(Pos_Hub,Orientation,[x0;y0;z0],[0;0;0]);
        X2_LiD=reshape(Pos_LiD(1,:),size(U_plot(:,:,1)));
        Y2_LiD=reshape(Pos_LiD(2,:),size(U_plot(:,:,1)));
        Z2_LiD=reshape(Pos_LiD(3,:),size(U_plot(:,:,1)));
        pcolor3(X2_LiD,Y2_LiD,Z2_LiD,HOM_plot{order==2}(:,:,ID_X).^0.5./U_plot(:,:,ID_X)*100);
        hold on;shading flat;alpha 0.5;colormap coolwarm;
    end
    view([244 47]);grid on;ax0.View=[30 30];
    xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');c=colorbar;
    c.Label.String='$\sqrt{\overline{u''^2}}/\overline{u}$ [$\%$]';
    axis equal;plot3(x0,y0,z0,'.b','markersize',15);plot3(LiDAR_loc(1),LiDAR_loc(2),LiDAR_loc(3),'xk','markersize',15,'linewidth',2);
    caxis([0 30]);
    TNR(13);
    
    %higher order moments
    for ID_order=1:length(order)
        if order(ID_order)>2
            mkfig('max');
            for ID_X=1:length(X2_plot(1,1,:))
                Pos_Hub=[hor(X2_plot(:,:,ID_X));hor(Y2_plot(:,:,ID_X));hor(Z2_plot(:,:,ID_X))];
                Pos_LiD=Hub2LiD_2(Pos_Hub,Orientation,[x0;y0;z0],[0;0;0]);
                X2_LiD=reshape(Pos_LiD(1,:),size(U_plot(:,:,1)));
                Y2_LiD=reshape(Pos_LiD(2,:),size(U_plot(:,:,1)));
                Z2_LiD=reshape(Pos_LiD(3,:),size(U_plot(:,:,1)));
                pcolor3(X2_LiD,Y2_LiD,Z2_LiD,HOM_plot{ID_order}(:,:,ID_X)./HOM_plot{order==2}(:,:,ID_X).^(0.5*order(ID_order)));
                hold on;shading flat;alpha 0.5;colormap coolwarm;
            end
            view([244 47]);grid on;ax0.View=[30 30];
            xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');c=colorbar;
            c.Label.String=['$\mu^{',num2str(order(ID_order)),'}/(\sqrt{\mu^2})^{',num2str(order(ID_order)),'}$'];
            axis equal;plot3(x0,y0,z0,'.b','markersize',15);plot3(LiDAR_loc(1),LiDAR_loc(2),LiDAR_loc(3),'xk','markersize',15,'linewidth',2);
            caxis(prctile(vert(HOM_plot{ID_order}(:,:,ID_X)./HOM_plot{order==2}(:,:,ID_X).^(0.5*order(ID_order))),[1 90]));
            TNR(13);
        end
    end  
catch
end
