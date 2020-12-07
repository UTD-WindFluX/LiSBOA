function s=vec2str(v,sep,format)

    if nargin==2
        format='%g';
    end
    s='';
    for i=1:length(v)
        
       s=[s,sep,num2str(v(i),format)]; 
    end
    s=s(length(sep)+1:end);
end