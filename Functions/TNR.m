function TNR(fontsize_fig)
    objs=get(gcf,'Children');
    for ID_obj=1:length(objs)
        if strcmp(objs(ID_obj).Type,'uimenu')==0
            
            c=objs(ID_obj).Children;
            for ID_c=1:length(c)
                try
                    c(ID_c).FontName='TimesNewRoman';
                    c(ID_c).Interpreter='latex';
                    c(ID_c).FontSize=fontsize_fig;
                catch
                end

            end
        
            l=[]; %#ok<*NASGU>
            try l=get(objs(ID_obj),'Layer');catch err;end
            if ~isempty(l)
                set(objs(ID_obj),'FontName','Times New Roman','FontSize',fontsize_fig);
            end
            l=[]; %#ok<*NASGU>
            try l=get(objs(ID_obj),'Layer');catch err;end
            if ~isempty(l);objs(ID_obj).Layer='top';end

            l=[];
            try l=get(objs(ID_obj),'XLabel');catch err;end
            if ~isempty(l);l.Interpreter='latex';end

            l=[];
            try l=get(objs(ID_obj),'YLabel');catch err;end
            if ~isempty(l);l.Interpreter='latex';end

            l=[];
            try l=get(objs(ID_obj),'ZLabel');catch err;end
            if ~isempty(l);l.Interpreter='latex';end

            l=[];
            try l=get(objs(ID_obj),'Text');catch err;end
            if ~isempty(l) && sum(l)~=0;l.Interpreter='latex';end

            l=[];
            try l=get(objs(ID_obj),'Title');catch err;end
            if ~isempty(l)
                if ~isempty(objs(ID_obj).Title.String)
                    if iscell(objs(ID_obj).Title.String)
                        for ID_line=1:length(objs(ID_obj).Title.String)
                            objs(ID_obj).Title.String{ID_line}=[objs(ID_obj).Title.String{ID_line}];
                        end
                    else
                        objs(ID_obj).Title.String=[objs(ID_obj).Title.String];
                    end
                end
                l.Interpreter='latex';
            end

            l=[];
            try l=get(objs(ID_obj),'Interpreter');catch err;end
            if ~isempty(l);objs(ID_obj).Interpreter='latex';end
        end
    end
    
    objs=findall(gcf,'type','annotation');
    
    for ID_obj=1:length(objs)
        c=objs(ID_obj).Children;
        for ID_c=1:length(c)
            try
                c(ID_c).FontName='TimesNewRoman';
                c(ID_c).Interpreter='latex';
                c(ID_c).FontSize=fontsize_fig;
            catch
            end
               
        end
    end
    
    warn = warning('query','last');
    if ~isempty(warn)
        id=warn.identifier;
        if ~isempty(id)
            warning('off',id);
        end
    end
end