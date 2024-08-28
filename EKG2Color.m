function [wave_color, ekg_color, color]=EKG2Color(EKG_Wave)

    
    if EKG_Wave=="p"
        wave_color="ro";
        ekg_color="r-";
        color="r";
    elseif EKG_Wave=="qrs"
        wave_color="bo";
        ekg_color="b-";
        color="b";
    elseif EKG_Wave=="t"
        wave_color="ko";
        ekg_color="k-";
        color="k";
    else
        fprintf("Region is not Labeled.  Plotting in default (blue)\n");
        wave_color="bo";
        ekg_color="b-";
        color="y";
    end


end