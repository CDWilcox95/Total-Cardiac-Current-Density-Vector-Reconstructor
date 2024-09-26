function ax2=PlotBSM(ax, Vq,R, s0, event_frames)
    title(ax, "Body Surface Map");
    xlabel(ax,"X (m)");   ylabel(ax,"Y (m)");   zlabel(ax,"Z (m)");
    clbr=colorbar;
    clbr.Label.String = 'Voltage (mV)';

    for s=event_frames
        Vqs=Vq(:,:,s);
        colormap('jet');
        bg=mean(Vqs(:));
        MapToCylinder(R, Vqs, bg);
        clim([0.75*min(Vq(:)), 0.75*max(Vq(:))]);
%         delete(ax);
    end

    ax2=ax;
end