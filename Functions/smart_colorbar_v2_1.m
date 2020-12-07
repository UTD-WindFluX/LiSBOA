function c=smart_colorbar_v2_1(values,labels,colormap_name,label)
    ax_ref=gca;
    ax=axes('Visible','off');
   
    c=colorbar; 
    c.Ticks=(values-min(values))/(max(values)-min(values));
    c.TickLabels=labels;
    c.Label.String=label;
    colormap(c,colormap_name)
    ax.Position=ax_ref.Position;
    ax.XLim=ax_ref.XLim;
    ax.YLim=ax_ref.YLim;
    ax.XScale=ax_ref.XScale;
    ax.YScale=ax_ref.YScale;
    axes(ax_ref);
end