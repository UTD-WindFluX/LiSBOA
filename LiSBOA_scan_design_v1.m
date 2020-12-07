%Scan design of volumetric scan through the LiSBOA
%23/11/2020: created
%12/02/2020: finalized

restoredefaultpath
clear all
close all
addpath('./Functions');

%% Inputs

%Pareto front parameters
vec_dtheta=[0.5 1 1.5 2 4];%[°] azimuth-angle resolutions tested
vec_dbeta=[0.5 1 2 3];%[°] elevation-angle resolutions tested
vec_sigma=[1/4 1/6 1/8 1/17];%smoothing paramters tested; use [1/4 1/6 1/8 1/17] for 3D calculations;

%known LiDAR parameters
Azi_lim=[30 60];%[°] azimuth-angle range
Ele_lim=[0 10];%[°] elevation-angle range
Range_lim=[100 1500];%[m] range limits
GL=25;%[m] gate length
acc_time=0.5;%[s] accumulation time
LiDAR_loc=[0 0 0];%[m] location of the LiDAR (x-y-z)

%flow characteristics
T=1200;%[s] total sampling time
tau=5;%[s] characteristic integral time scale
u_std=1.5;%[m/s] characteristic velocity standard deviation
Orientation=45;%[°] wind direction (0°=North,90°=East,180°=South,270°=West); the x axis will be aligned with the stremwise direction

%origin of the Cartesian coordinates
x0=800;%[m]
y0=800;%[m]
z0=80;%[m]

%Cartesian domain size (in the domain aligned with wind direction and centered in [x0,y0,z0])
xmin=0;%[m]
xmax=750;%[m];
ymin=-150;%[m]
ymax=150;%[m]
zmax=100;%[m]
zmin=-10;%[m]

%LiSBOA parameters
Dn_x=250;%[m] fundamental half-wavelength in x (streamwise)
Dn_y=50; %[m] fundamental half-wavelength in y (spanwise)
Dn_z=50; %[m] fundamental half-wavelength in z (vertical)
grid_factor=0.25;%ratio of the resolution of the Cartesian to the fundamental half-wavelength vector
R_max=3;%[sigmas] radius of influence
max_Dd=1;%maximum local data spacing (in the scaled frame of reference)
tol_dist=1/14;%tolerance of the spatial location of samples in the scaled frame of reference, points closer than this are assumed as unique

%% Initialization
Range=Range_lim(1):GL:Range_lim(2);

%% Pareto front
for ID_dbeta=1:length(vec_dbeta)
    dbeta=vec_dbeta(ID_dbeta);
    Ele=Ele_lim(1):dbeta:Ele_lim(2);
    
    for ID_dtheta=1:length(vec_dtheta)
        dtheta=vec_dtheta(ID_dtheta);
        
        %scan geometry
        Azi=Azi_lim(1):dtheta:Azi_lim(2);
        [A,E,R]=ndgrid(Azi,Ele,Range);
        X{ID_dbeta,ID_dtheta}=cosd(E).*cosd(90-A).*R+LiDAR_loc(1);
        Y{ID_dbeta,ID_dtheta}=cosd(E).*sind(90-A).*R+LiDAR_loc(2);
        Z{ID_dbeta,ID_dtheta}=sind(E).*R+LiDAR_loc(3);
        
        %ref. system change (aligned with wind direction, centered in [x0, y0, z0])
        Pos_LiD=[hor(X{ID_dbeta,ID_dtheta});hor(Y{ID_dbeta,ID_dtheta});hor(Z{ID_dbeta,ID_dtheta})];
        Pos_Hub=LiD2Hub_2(Pos_LiD,Orientation,[x0;y0;z0],LiDAR_loc(:));
        X_Hub=reshape(Pos_Hub(1,:),size(X{ID_dbeta,ID_dtheta}));
        Y_Hub=reshape(Pos_Hub(2,:),size(Y{ID_dbeta,ID_dtheta}));
        Z_Hub=reshape(Pos_Hub(3,:),size(Z{ID_dbeta,ID_dtheta}));
        
        sampling_time=(acc_time*length(Azi))*length(Ele);%[s] sampling time of each scan
        
        L(ID_dbeta,ID_dtheta)=floor(T/sampling_time);%number of scan repetitions
        
        %Dd
        for ID_sigma=1:length(vec_sigma)
            sigma=vec_sigma(ID_sigma);
            [X2{ID_dbeta,ID_dtheta},Y2{ID_dbeta,ID_dtheta},Z2{ID_dbeta,ID_dtheta},excl{ID_dbeta,ID_dtheta,ID_sigma},Dd{ID_dbeta,ID_dtheta,ID_sigma}]=...
                LiSBOA(vert(X_Hub),vert(Y_Hub),vert(Z_Hub),xmin,xmax,ymin,ymax,zmin,zmax,Dn_x,Dn_y,Dn_z,...
                sigma,R_max,grid_factor,tol_dist,max_Dd,'distance',[],[],[]);
            
            epsilon1{ID_dbeta}(ID_dtheta,ID_sigma)=sum(vert(excl{ID_dbeta,ID_dtheta,ID_sigma}))./numel(X2{ID_dbeta,ID_dtheta})*100;%data loss percentage
        end
        
        p=1:(L(ID_dbeta,ID_dtheta)-1);
        epsilon2{ID_dbeta}(ID_dtheta)=(1/L(ID_dbeta,ID_dtheta)+2/L(ID_dbeta,ID_dtheta)^2*(sum(exp(-p*sampling_time/tau).*(L(ID_dbeta,ID_dtheta)-p)))).^0.5*u_std;%STD of the mean
        
        disp(['dbeta = ',num2str(dbeta),'°, dtheta = ',num2str(dtheta),'° done']);
    end
    
end

%% Plot
mkfig('max');
col=coolwarm(length(vec_dtheta));
markers={'o','^','s','p'};

N_plot1=floor(length(vec_dbeta)^0.5);
N_plot2=ceil(length(vec_dbeta)/N_plot1);
for ID_dbeta=1:length(vec_dbeta)
    ax_Pareto(ID_dbeta)=subplot(N_plot1,N_plot2,ID_dbeta);
    for ID_sigma=1:length(vec_sigma)
        scatter(epsilon1{ID_dbeta}(:,ID_sigma),epsilon2{ID_dbeta},50,col,'filled','marker',markers{ID_sigma},...
            'markeredgecolor','k','markerfacealpha',0.5,'linewidth',0.75);hold on;
    end
    
    if ID_dbeta==1
        clear pl vec_legend
        for ID_sigma=1:length(vec_sigma)
            pl(ID_sigma)=plot(1000,1,'ok','markersize',7,'marker',markers{ID_sigma},'linewidth',1);hold on;
            vec_legend{ID_sigma}=['$\sigma$ = 1/',num2str(1/vec_sigma(ID_sigma))];
        end
        legend(pl,vec_legend);
    end
    grid on;xlabel('$\epsilon^I$ [\%]');ylabel('$\epsilon^{II}$ [m s$^{-1}$]');
    
    smart_colorbar_disc(array2cell(vec_dtheta),'coolwarm','$\Delta \theta$ [$^\circ$]');box on
    xlim([0 100]);
    ylim([0 max([epsilon2{:}])]);
    
    title(['$\Delta \beta = ',num2str(vec_dbeta(ID_dbeta)),'^\circ$']);
    xticks([0 25 50 75 100]);
    TNR(13)
end
i=inputdlg('Optimal elevation-angle resolution');

try
    ID_dbeta_sel=find(str2double(i{1})==vec_dbeta);
    axes(ax_Pareto(ID_dbeta_sel));
    ax_Pareto(ID_dbeta_sel).Color=[0.9 1 0.9];
    msgbox('Select optimal case in the highlighted plot and press Enter');
    
    [x,y]=ginput;
    dist=((epsilon1{ID_dbeta_sel}-x).^2/max(vert(epsilon1{ID_dbeta_sel}))^2+(repmat(vert(epsilon2{ID_dbeta_sel}),1,length(vec_sigma))-y).^2/max(epsilon2{ID_dbeta_sel})^2).^0.5;
    [ID_dtheta_sel,ID_sigma_sel]=find(dist==min(vert(dist)));
    plot(epsilon1{ID_dbeta_sel}(ID_dtheta_sel,ID_sigma_sel),epsilon2{ID_dbeta_sel}(ID_dtheta_sel),'or','markersize',15);
    
    %plot scan geometry
    mkfig('hor');
    subplot(1,2,1);
    
    X_plot=X{ID_dbeta_sel,ID_dtheta_sel};
    Y_plot=Y{ID_dbeta_sel,ID_dtheta_sel};
    Z_plot=Z{ID_dbeta_sel,ID_dtheta_sel};
    
    scatter3(vert(X_plot),vert(Y_plot),vert(Z_plot),2,'k','filled','markerfacealpha',0.2);axis equal;grid on;box on;
    xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');hold on;
    
    Pos_LiD=[reshape(X_plot,1,[]);reshape(Y_plot,1,[]);reshape(Z_plot,1,[])];
    Pos_Hub=LiD2Hub_2(Pos_LiD,Orientation,[x0;y0;z0],LiDAR_loc(:));
    X_Hub=reshape(Pos_Hub(1,:),size(X_plot));
    Y_Hub=reshape(Pos_Hub(2,:),size(Y_plot));
    Z_Hub=reshape(Pos_Hub(3,:),size(Z_plot));
    X_plot(X_Hub<xmin)=nan;X_plot(X_Hub>xmax)=nan;
    Y_plot(Y_Hub<ymin)=nan;Y_plot(Y_Hub>ymax)=nan;
    Z_plot(Z_Hub<zmin)=nan;Z_plot(Z_Hub>zmax)=nan;
    
    scatter3(vert(X_plot),vert(Y_plot),vert(Z_plot),2,'r','filled');axis equal;grid on;box on;
    xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');
    plot3(x0,y0,z0,'.b','markersize',15);plot3(LiDAR_loc(1),LiDAR_loc(2),LiDAR_loc(3),'xk','markersize',15,'linewidth',2);
    
    drawbox(xmin,xmax,ymin,ymax,zmin,zmax,Orientation,[x0;y0;z0],LiDAR_loc(:),'b',0.1,'k',0.75);
    
    l=legend('All points','Points in the domain','location','northeast');
    ax0=gca;
    
    %plot local spacing
    subplot(1,2,2);
    Dd_sel=Dd{ID_dbeta_sel,ID_dtheta_sel,ID_sigma_sel};
    Dd_sel(excl{ID_dbeta_sel,ID_dtheta_sel,ID_sigma_sel})=nan;
    
    Dd_plot=permute(Dd_sel,[3 1 2]);
    X2_plot=permute(X2{ID_dbeta_sel,ID_dtheta_sel},[3 1 2]);
    Y2_plot=permute(Y2{ID_dbeta_sel,ID_dtheta_sel},[3 1 2]);
    Z2_plot=permute(Z2{ID_dbeta_sel,ID_dtheta_sel},[3 1 2]);
    
    for ID_X=1:length(X2{ID_dbeta_sel,ID_dtheta_sel}(1,1,:))
        Pos_Hub=[hor(X2_plot(:,:,ID_X));hor(Y2_plot(:,:,ID_X));hor(Z2_plot(:,:,ID_X))];
        Pos_LiD=Hub2LiD_2(Pos_Hub,Orientation,[x0;y0;z0],[0;0;0]);
        X2_LiD=reshape(Pos_LiD(1,:),size(Dd_plot(:,:,1)));
        Y2_LiD=reshape(Pos_LiD(2,:),size(Dd_plot(:,:,1)));
        Z2_LiD=reshape(Pos_LiD(3,:),size(Dd_plot(:,:,1)));
        pcolor3(X2_LiD,Y2_LiD,Z2_LiD,Dd_plot(:,:,ID_X));
        hold on;shading flat;alpha 0.5;colormap coolwarm;
    end
    view([30 30]);zlim([zmin zmax]);grid on;ax0.View=[30 30];
    xlabel('$x$ [m]');ylabel('$y$ [m]');zlabel('$z$ [m]');c=colorbar;c.Label.String='$\Delta \tilde{d}$';caxis([0 1]);c.Location='east';c.Position(1)=c.Position(1)+0.07;
    axis equal;axis([ax0.XLim ax0.YLim ax0.ZLim]);plot3(x0,y0,z0,'.b','markersize',15);plot3(LiDAR_loc(1),LiDAR_loc(2),LiDAR_loc(3),'xk','markersize',15,'linewidth',2);
    ax=axes('Position',[0 0.6 0.1 0.1],'XColor','none','YColor','none','Color','none');
    t=text(0.01,1,{['Azimuth resolution = ',num2str(vec_dtheta(ID_dtheta_sel)),'$^\circ$'];...
        ['Azimuth range = [',vec2str(Azi_lim,', '),']$^\circ$'];...
        ['Elevation resolution = ',num2str(vec_dbeta(ID_dbeta_sel)),'$^\circ$'];...
        ['Elevation range = [',vec2str(Ele_lim,', '),']$^\circ$'];...
        ['Gate length = ',num2str(GL),' m'];['Detection range = [',vec2str(Range_lim,','),'] m'];...
        ['Accumulation time = ',num2str(acc_time),' s'];...
        ['$\sigma$ = ',num2str(vec_sigma(ID_sigma_sel))]});
    t.FontSize=10;
    TNR(13);
    
catch
    error('Invalid value.');
end
