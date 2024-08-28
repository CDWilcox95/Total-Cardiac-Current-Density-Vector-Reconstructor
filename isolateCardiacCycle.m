function [first_idx, end_idx, int_cycles]=isolateCardiacCycle(mappedEKGWave, num_cycles)
    temp_wave=mappedEKGWave;    
    int_cycles=zeros(num_cycles-1,1);
    
    first_idx=-1; end_idx=-1;  cycle_complete=false; lbl0="";  cycle_cnt=0; 
    for i=1:length(mappedEKGWave)

        lbl=mappedEKGWave(i);
        if lbl=="p" && first_idx==-1
            first_idx=i;lbl0="p";
        end

        if lbl0=="t" && lbl=="p"
            cycle_complete=true;
            cycle_cnt=cycle_cnt+1;
            if cycle_cnt<num_cycles
                int_cycles(cycle_cnt)=i;
            end
        end

        if cycle_complete && cycle_cnt==num_cycles
            end_idx=i;
            break;
        end
        lbl0=lbl;
    end

    if end_idx==-1
        end_idx=length(mappedEKGWave);
    end

end