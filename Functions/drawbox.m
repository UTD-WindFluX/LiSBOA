function drawbox(xmin,xmax,ymin,ymax,zmin,zmax,Orientation,pos_Hub,pos_LiD,facecolor,facealpha,edgecolor,linewidth)

    vertex_Hub=[xmin xmin xmax xmax xmin xmin xmax xmax;ymin ymax ymax ymin ymin ymax ymax ymin;zmin zmin zmin zmin zmax zmax zmax zmax];

    vertex_LiD=Hub2LiD_2(vertex_Hub,Orientation,pos_Hub,pos_LiD);
    
    patch(vertex_LiD(1,1:4),vertex_LiD(2,1:4),vertex_LiD(3,1:4),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);
    hold on;
    patch(vertex_LiD(1,5:8),vertex_LiD(2,5:8),vertex_LiD(3,5:8),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);

    patch(vertex_LiD(1,[1 4 8 5]),vertex_LiD(2,[1 4 8 5]),vertex_LiD(3,[1 4 8 5]),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);
    patch(vertex_LiD(1,[2 3 7 6]),vertex_LiD(2,[2 3 7 6]),vertex_LiD(3,[2 3 7 6]),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);
    
    patch(vertex_LiD(1,[4 3 7 8]),vertex_LiD(2,[4 3 7 8]),vertex_LiD(3,[4 3 7 8]),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);
    patch(vertex_LiD(1,[1 2 6 5]),vertex_LiD(2,[1 2 6 5]),vertex_LiD(3,[1 2 6 5]),edgecolor,'linewidth',linewidth,'facecolor',facecolor,'facealpha',facealpha);
end