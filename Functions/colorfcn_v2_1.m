function col=colorfcn_v2_1(x,xmin,xmax,color_map)
    cols=feval(color_map,100);
    x(x>xmax)=xmax;x(x<xmin)=xmin;
    for i=1:3 
        col(:,i)=interp1(0:1/(length(cols(:,1))-1):1,cols(:,i),(x(:)-xmin)/(xmax-xmin));
    end
end