%discrete colorbar

function c=smart_colorbar_disc(ticks,colormap_name,label)
    ax_ref=gca; 
    ax=axes('Visible','off');
    N=length(ticks);
    c=colorbar; 
    c.Ticks=1/N/2:1/N:1;
    c.TickLabels=ticks;
    c.Label.String=label;
    colormap(c,[colormap_name,'(',num2str(N),')'])
    ax.Position=ax_ref.Position;
    ax.XLim=ax_ref.XLim;
    ax.YLim=ax_ref.YLim;
    ax.XScale=ax_ref.XScale;
    ax.YScale=ax_ref.YScale;
    axes(ax_ref);
end