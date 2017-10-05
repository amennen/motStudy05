v=1;

for iTrial = 1:216
    
    zscoreLen = double(iTrial);
    zscoreLen1 = double(iTrial - 1);
    zscoreConst = 1.0/zscoreLen;
    zscoreConst1 = 1.0/zscoreLen1;
    
    patterns.realtimeMean(1,v) = mean(z0.patterns.raw_sm_filt(1:iTrial,v),1);
    
    patterns.realtimeStd(1,v) = std(z0.patterns.raw_sm_filt(1:iTrial,v),1,1); %flad to use N instead of N-1
    patterns.realtimeMean(1,v) = (z0.patterns.realtimeMean(1,v).*zscoreLen1 + z0.patterns.raw_sm_filt(iTrial,v)).*zscoreConst;
    if iTrial > 64
        if iTrial == 70
            dasd
        end
        patterns.raw_sm_filt_z(iTrial,v) = (z0.patterns.raw_sm_filt(iTrial,v) - patterns.realtimeMean(1,v))./patterns.realtimeStd(1,v);
    else
        patterns.raw_sm_filt_z(iTrial,v) = (z0.patterns.raw_sm_filt(iTrial,v) - patterns.realtimeMean(1,v))./z0.patterns.lastStd(1,v);
        
    end
    
end