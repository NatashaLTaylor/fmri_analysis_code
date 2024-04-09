function dsmtx = dsmtx_hrf(onsets,tr)

    hrf = spm_hrf(tr);
    
    nOnsets = size(onsets,1);
    
    onsets_tr = onsets./(tr);
    
    time_temp = zeros(round(max(onsets_tr)),1);
    
    for n = 1:nOnsets
        time_temp(round(onsets_tr(n)),1) = 1;
    end

    dsmtx = conv(hrf,time_temp);

end