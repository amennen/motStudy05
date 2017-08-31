subjectVec = 1:23;
bad = [2 7 15];
subjectVec(bad) = [];
nResp = zeros(1,length(subjectVec));
for s = 1:length(subjectVec)
    SESSION = 19;
    base_path = pwd;
    subjectNum = subjectVec(s);
    behavioral_dir = [base_path '/' 'BehavioralData/' num2str(subjectNum) '/'];
    r = dir(fullfile(behavioral_dir, ['EK' num2str(SESSION) '_' 'SUB'  '*.mat']));
    r = load(fullfile(behavioral_dir,r(end).name));
    trials = table2cell(r.datastruct.trials);
    resp = cell2mat(trials(:,12));
    NR = find(isnan(resp));
    resp(NR) = -1;
    
    cond = cell2mat(trials(:,9));
    TL = zeros(1,length(resp));
    Targ = find(cond<3);
    Lure = find(cond>2);
    resp_targ = resp(Targ);
    resp_lure = resp(Lure);
    nResp(s) = length(find(resp_lure==-1));
    bins = -1:5;
%         figure;
%        
%         subplot(1,2,1)
%         nT = hist(resp_targ,bins);
%         distT = nT/(16*4);
%         title('Target Resp')
%         bar(bins,distT);
%         ylim([0 1])
%         subplot(1,2,2)
%         nL = hist(resp_lure,bins);
%         distL = nL/(16*4);
%         title('Lure Resp')
%         bar(bins,distL);
%         ylim([0 1])
%         
%         title(['Subject ' num2str(subjectVec(s))]);
end

m_r = mean(nResp);
std_r = std(nResp);

subjectVec
nResp
m_r + 2*std_r
