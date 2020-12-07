%LiSBOA for scan design and statistics reconstruction, 2D and 3D

%_________________________________________________________________________________
%Inputs
%x_exp, y_exp, z_exp: location of the samples in x, y, z
%xmin, xmax, ymin, ymax, zmin, zmax: domain boundaries
%Dn_x, Dn_y, Dn_y: x,y,z components of the fundemental half-wavelenght vector
%sigma: smoothing parameter
%R_max: influence radius (in sigmas, 3 recommended)
%grid factor: %ratio of the resolution of the Cartesian to the fundamental half-wavelength vector (0.25 recommended)
%tol_dist: %tolerance of the spatial location of samples in the scaled frame of reference, points closer than this are assumed as unique (<0.1 recommended)
%maximum local data spacing (in the scaled frame of reference)
%option: 
%%%%%%%%'distance'-> only calculates local data spacing, for Pareto fronts;
%%%%%%%%'all'->calculates also statistics (mean and HOM)
%F: input scalar field 
%order: order of the central statistical moment (2 for variance, 3 for skewness, etc..)
%max_iter: maximum number of iterations

%_________________________________________________________________________________
%Outputs
%X2, Y2, Z2: Cartesian grid
%excl: excluded points due to Peterson-Middleton violation
%Dd: local data spacing
%Avg_all: mean field for all iterations
%HOM_all: central statistical moment for all iterations
%MAD: mean absolute difference of the mean field at each iteration (as a convergence monitor)
%select: structure containing the index of the selected point for each cell (for more in-depth postprocessing)

%_________________________________________________________________________________
%Versioning
%11/23/2020: created, finalized
%12/04/2020: 2D option added, finalzied

function [X2,Y2,Z2,excl,Dd,Avg_all,HOM_all,MAD,select,w]=LiSBOA(x_exp,y_exp,z_exp,xmin,xmax,ymin,ymax,zmin,zmax,Dn_x,Dn_y,Dn_z,sigma,R_max,grid_factor,tol_dist,max_Dd,option,F,order,max_iter)

    %% Initialization
    
    %nans removal
    if strcmp(option,'all')
        reals=find(~isnan(F).*~isnan(x_exp).*~isnan(y_exp).*~isnan(z_exp));
        F=F(reals);
    else
        reals=find(~isnan(x_exp).*~isnan(y_exp).*~isnan(z_exp));
    end
    x_exp=x_exp(reals);
    y_exp=y_exp(reals);
    z_exp=z_exp(reals);

    %rescaling of maximum radii
    R_max=R_max*sigma;

    A=pi*R_max^2;%area of integration (for 2D only)
    V=4/3*pi*R_max^3;%volume of integration (for 3D only)

    %data centroid
    xc=nanmean(x_exp);
    yc=nanmean(y_exp);
    zc=nanmean(z_exp);

    %sampling points location normalization
    x=(x_exp-xc)/Dn_x;
    y=(y_exp-yc)/Dn_y;
    z=(z_exp-zc)/Dn_z;

    %numerical grid
    dx=grid_factor*Dn_x; 
    dy=grid_factor*Dn_y;
    dz=grid_factor*Dn_z;
    if isinf(dx)
        dx=10^99;
    end
    if isinf(dy)
        dy=10^99;
    end
    if isinf(dz)
        dz=10^99;
    end

    %grid
    X_bin=(((xmin-dx/2):dx:(xmax+dx))-xc)/Dn_x;
    Y_bin=(((ymin-dy/2):dy:(ymax+dy))-yc)/Dn_y;
    Z_bin=(((zmin-dz/2):dz:(zmax+dz))-zc)/Dn_z;

    X_vec=mid(X_bin);
    Y_vec=mid(Y_bin);
    Z_vec=mid(Z_bin);

    if isempty(X_vec)
        X_vec=0;
    end
    if isempty(Y_vec)
        Y_vec=0;
    end
    if isempty(Z_vec)
        Z_vec=0;
    end

    [X,Y,Z]=meshgrid(X_vec,Y_vec,Z_vec);
    [Ny,Nx,Nz]=size(X);

    %display grid
    X2=X*Dn_x+xc;
    Y2=Y*Dn_y+yc;
    Z2=Z*Dn_z+zc;

    %variable initialization
    Ddp=inf(Ny,Nx,Nz);    
    selectp=cell(Ny,Nx,Nz);
    wp=cell(Ny,Nx,Nz);

%     %parallel computation
%     p=gcp('nocreate');
%     if isempty(p)
%         parpool
%     end
%     ppm = ParforProgMon('Calculation of LiSBOA weights: ', Nx, 1, 1000, 200);

    %% Main
    %weights and spacing
    progbar = waitbar(0,'Calculating LiSBOA weights...');
    for ID_X=1:Nx
        x0=(xmin+dx*(ID_X-1)-xc)/Dn_x;
        waitbar(ID_X/Nx,progbar,'Calculating LiSBOA weights...');
        for ID_Y=1:Ny
            y0=(ymin+dy*(ID_Y-1)-yc)/Dn_y;
            for ID_Z=1:Nz
                z0=(zmin+dz*(ID_Z-1)-zc)/Dn_z;

                rSq=((x0-x).^2+(y0-y).^2+(z0-z).^2);
                selectp{ID_Y,ID_X,ID_Z}=rSq<R_max^2; 
                xs=x(selectp{ID_Y,ID_X,ID_Z});
                ys=y(selectp{ID_Y,ID_X,ID_Z});
                zs=z(selectp{ID_Y,ID_X,ID_Z});
                
                if strcmp(option,'all')
                    %weights of first iteration
                    f_x=exp(-(x0-xs).^2/(2*sigma^2));
                    f_y=exp(-(y0-ys).^2/(2*sigma^2));
                    f_z=exp(-(z0-zs).^2/(2*sigma^2));
                    wp{ID_Y,ID_X,ID_Z} =  f_x.*f_y.*f_z;
                else 
                    wp{ID_Y,ID_X,ID_Z}=nan;
                end
 
                %exclude overlapping points
                Pos_all=round2([vert(xs-x0),vert(ys-y0),vert(zs-z0)],tol_dist);
                Pos_uni=unique(Pos_all,'rows');
                M_uni=length(Pos_uni(:,1));

                if M_uni>1
                    %Koch with uni filter
                    if dx==10^99||dy==10^99||dz==10^99
                        Ddp(ID_Y,ID_X,ID_Z)=A^(1/2)/(M_uni^(1/2)-1);
                    else
                        Ddp(ID_Y,ID_X,ID_Z)=V^(1/3)/(M_uni^(1/3)-1);
                    end
                end
            end
        end
%         ppm.increment();  %#ok<PFBNS>
    end
%     ppm.delete()
    close(progbar);
    Dd=Ddp;
    
    if max_Dd~=1;warning('Maximum Dd is not equal to the recommended 1');end

    select=selectp;
    w=wp;

    %filter
    undersampled=Dd>max_Dd;
    x_un=hor(X(undersampled));
    y_un=hor(Y(undersampled));
    z_un=hor(Z(undersampled));
    
    tic
    if sum(vert(undersampled))>0
        for ID_X=1:Nx
            x0=(xmin+dx*(ID_X-1)-xc)/Dn_x;
            for ID_Y=1:Ny
                y0=(ymin+dy*(ID_Y-1)-yc)/Dn_y;
                for ID_Z=1:Nz
                    z0=(zmin+dz*(ID_Z-1)-zc)/Dn_z;
                    vec_distSq=(x0-x_un).^2+(y0-y_un).^2+(z0-z_un).^2;
                    DistSq2(ID_Y,ID_X,ID_Z)=min(vec_distSq);
                end
            end
        end
        edges=reshape(DistSq2<R_max^2,size(Dd));
        excl=edges>0;
    else
        excl=zeros(size(Dd))>inf;
    end
        
    if strcmp(option,'all')

        %iterative statistics
        Avgp_old=zeros(Ny,Nx,Nz);
        MAD(1)=Inf;
        Err_exp=F;

        for Iter=1:max_iter+2
            for ID_order=1:length(order)
                HOMp{ID_order}=nan(Ny,Nx,Nz);
            end
            for ID_X=1:Nx
                for ID_Y=1:Ny
                    for ID_Z=1:Nz
                        if excl(ID_Y,ID_X,ID_Z)==0
                                                       
                            %extract data
                            ff=Err_exp(select{ID_Y,ID_X,ID_Z});
                            ww=w{ID_Y,ID_X,ID_Z};

                            %nan removal
                            ww(isnan(ff))=0;
                            ff(isnan(ff))=0;
                            ww=ww/sum(ww+eps);
                            WM=sum(ff.*ww);%weighted mean

                            %weighted HOM
                            if Iter>1
                                for ID_order=1:length(order)
                                    HOMp{ID_order}(ID_Y,ID_X,ID_Z)=sum(ff.^order(ID_order).*ww);%central statistical moment
                                end
                                if isnan(WM);WM=0;end
                            end                            
                        else
                            WM=0;
                        end
                        Avgp(ID_Y,ID_X,ID_Z)=Avgp_old(ID_Y,ID_X,ID_Z)+WM;%average update                        
                    end
                end
            end
            if Iter<=max_iter+1
                Avg_all{Iter}=Avgp;%average save
                Avgp_old=Avgp;

                %correction
                Avg_filt=Avg_all{Iter};
                Avg_filt(excl)=nan;           
                if Nx==1
                    F_interp=interp2(permute(Y,[3 2 1]),permute(Z,[3 2 1]),permute(Avg_filt,[3 2 1]),y,z,'linear');
                elseif Ny==1
                    F_interp=interp2(permute(X,[3 1 2]),permute(Z,[3 1 2]),permute(Avg_filt,[3 1 2]),x,z,'linear'); 
                elseif Nz==1
                    F_interp=interp2(X,Y,Avg_filt,x,y,'linear');
                else
                    F_interp=interp3(X,Y,Z,Avg_filt,x,y,z,'linear');
                end
                Err_exp=F-F_interp;

                %residual
                if Iter>1
                    MAD(Iter)=mean(abs(reshape(Avg_all{Iter}-Avg_all{Iter-1},[],1)),'omitnan');
                end
                disp(['Iteration ',num2str(Iter-1),' - MAD = ',num2str(MAD(Iter))]);
            end
            if Iter>1
                for ID_order=1:length(order)
                    HOM_all{ID_order}{Iter-1}=HOMp{ID_order};%HOM save
                end
            end
        end
    else
        Avg_all=nan;
        HOM_all=nan;
        MAD=nan;
        select=nan;
    end
end

