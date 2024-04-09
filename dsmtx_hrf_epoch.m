function dsmtx = dsmtx_hrf_epoch(onsets,durations,tr)
%onset defines the start time of the event in secs
%duration defines the duration for the event in secs
%tr is the number of repeats
 

    hrf = spm_hrf(tr);

    nOnsets = size(onsets,1);
    onsets_tr = onsets./(tr);

    durations_tr = durations./(tr);

    time_temp = zeros(round(max(onsets_tr)),1);

    

    %for n = 1:nOnsets
      for n=1:nOnsets
        if onsets_tr(n) < 0.5
          time_temp(1:round(1+(durations_tr(n)))) =1;
       else 
         time_temp(round(onsets_tr(n)):round(onsets_tr(n)+durations_tr(n)),1) = 1;
       end
      end    
%       for n=2:nOnsets
%           time_temp(round(onsets_tr(n)):round(onsets_tr(n)+durations_tr(n)),1) = 1;
%       end 
      
%    for n = 1:nOnsets
%         time_temp(round(onsets_tr(n)):round(onsets_tr(n)+durations_tr(n)),1) = 1;
%     end
    
    dsmtx = conv(hrf,time_temp);

 end