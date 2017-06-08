% plot example classifier evidence and dot speed for RT and YC group 21, 29
projectName = 'motStudy02';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 15;
nblock = 3;
svec = [8 12 14 15 16 18 20:22 26 27 28 29 30];
RT = [8 12 14 15 18 21 22];
YC = [16 20 26 27 28 29 30];
RT_m = [8 12 14 15 18 21 22];
YC_m = [16 28 20 26 27 29 30];
iRT = find(ismember(svec,RT));
iYC = find(ismember(svec,YC));
%svec = 8:13;
nsub = length(svec);
sepbystim = zeros(nstim,nTRs*3);
speedbystim = zeros(nstim,nTRs*3);
MOT_PREP = 5;
colors = [207 127 102;130 161 171; 207 64 19]/255;

%colors = [110 62 106;83 200 212; 187 124 181]/255;
plotstim = 1; %if you want trial by trial plots
plotmixedstim = 0; %if you want trial by trial plots with mixed stimuli
allplotDir = ['/Data1/code/' projectName '/' 'Plots' '/' ];


subjectNum = 21;
stim = 8;

for iblock = 1:nblock
    
    blockNum = iblock;
    SESSION = 19 + blockNum;
    %blockNum = SESSION - 20 + 1;
    behavioral_dir = ['/Data1/code/' projectName '/' 'code' '/BehavioralData/' num2str(subjectNum) '/'];
    behavioral_dir = [fileparts(which('mot_realtime01.m')) '/BehavioralData/' num2str(subjectNum) '/'];
    save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/']; %this is where she sets the save directory!
    classOutputDir = fullfile(save_dir,['motRun' num2str(blockNum)], 'classOutput/');
    runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
    fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime01_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
    %get hard speed
    prep = dir([behavioral_dir 'mot_realtime01_' num2str(subjectNum) '_' num2str(MOT_PREP)  '*.mat']);
    prepfile = [behavioral_dir prep(end).name];
    lastRun = load(prepfile);
    hardSpeed(s) = 30 - lastRun.stim.tGuess(end);
    plotDir = ['/Data1/code/' projectName '/' 'Plots' '/' num2str(subjectNum) '/'];
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end
    matlabOpenFile = [behavioral_dir '/' fileSpeed(end).name];
    d = load(matlabOpenFile);
    allSpeed = d.stim.motionSpeed; %matrix of TR's
    speedVector = reshape(allSpeed,1,numel(allSpeed));
    allMotionTRs = convertTR(d.timing.trig.wait,d.timing.plannedOnsets.motion,d.config.TR); %row,col = mTR,trialnumber
    TRvector = reshape(allMotionTRs,1,numel(allMotionTRs));
    for i=1:length(TRvector)
        fileTR = TRvector(i); %+ 2;%take out so not shifted
        [~, tempfn{fileTR}] = GetSpecificClassOutputFile(classOutputDir,fileTR);
        tempStruct = load(fullfile(classOutputDir, tempfn{fileTR}));
        categsep(i) = tempStruct.classOutput;
    end
    run = dir([runHeader 'motpatternsdata_' num2str(SESSION) '*']);
    run = load(fullfile(runHeader,run(end).name));
    zcategsep = run.patterns.categsep(TRvector - 10 + 2); %minus 10 because we take out those 10
    %categsep =
    sepbytrial = reshape(categsep,15,10);
    sepbytrial = sepbytrial'; %results by trial number, TR number
    fbsepbytrial = sepbytrial(:,5:end);
    
    % sepbytrial = sepbytrial(:,5:end);%take only the ones once fb starts
    sepvec = reshape(sepbytrial,1,numel(sepbytrial));
    fbsepvec = reshape(fbsepbytrial, 1, numel(fbsepbytrial));
    
    speedbytrial = reshape(speedVector,nTRs,nstim);
    speedbytrial = speedbytrial';
    [~,indSort] = sort(d.stim.id);
    sepinorder = sepbytrial(indSort,:);
    speedinorder = speedbytrial(indSort,:);
    %test if fb only
    fbsepinorder = sepinorder(:,5:end);
    fbspeedinorder = speedinorder(:,5:end);
    nTRs2 = 11; %change back to 11 and sep... afterwards
    sepbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = sepinorder;
    speedbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = speedinorder;
    fbsepbystim(:,(iblock-1)*nTRs2 + 1: iblock*nTRs2 ) = fbsepinorder;
    fbspeedbystim(:,(iblock-1)*nTRs2 + 1: iblock*nTRs2 ) = fbspeedinorder;
    sepmixed(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = sepbytrial;
    speedmixed(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = speedbytrial;
    
    rep = 1:10;
    
    allspeeds = [allspeeds speedVector];
    allsep = [allsep sepvec];
    fbsep = [fbsep fbsepvec];
    
end

newspeedbystim = reshape(speedbystim,1,numel(speedbystim));
newsepbystim = reshape(sepbystim,1,numel(sepbystim));
[good] = find(newsepbystim > 0.05 & newsepbystim < 0.15);
goodSpeeds = newspeedbystim(good);

fbnewspeedbystim = reshape(fbspeedbystim,1,numel(fbspeedbystim));
fbnewsepbystim = reshape(fbsepbystim,1,numel(fbsepbystim));
[fbgood] = find(fbnewsepbystim > 0.05 & fbnewsepbystim < 0.15);
fbgoodSpeeds = fbnewspeedbystim(fbgood);

fbnewspeedbystim = reshape(fbspeedbystim,1,numel(fbspeedbystim));
fbnewsepbystim = reshape(fbsepbystim,1,numel(fbsepbystim));

thisfig = figure;
clf;
x = 1:nTRs*nblock;
subplot(2,1,1)
hAx(1) = plot(x,sepbystim(stim,:),'LineStyle', '-', 'Color', colors(2,:), 'LineWidth', 5)
hold on;
for rep = 1:2
    line([rep*nTRs+.5 rep*nTRs + .5], [-10 15], 'color', 'k', 'LineWidth', 2);
end
line([0 46], [0.1 0.1], 'color', [140 136 141]/255, 'LineWidth', 2.5,'LineStyle', '--');
ylim([-1 1])
xlim([1 45])

subplot(2,1,2)
xlim([1 45])
ylabel('Retrieval Evidence')
hAx(2) = plot(x,speedbystim(stim,:), 'LineStyle', '-', 'Color', colors(1,:), 'LineWidth', 5)
%[hAx,hLine1, hLine2] = plotyy(x,sepbystim(stim,:),x,speedbystim(stim,:));
xlabel('TR Number (2s)')
hold on;
for rep = 1:2
    line([rep*nTRs+.5 rep*nTRs + .5], [-10 15], 'color', 'k', 'LineWidth', 2);
end
ylim([0 3])
ylabel('Dot Speed')
%ylabel(hAx(2), 'Dot Speed', 'Color', 'k')
%ylabel(hAx(1), 'Category Evidence', 'Color', 'k')
%ylim(hAx(2),[-0.5 10])
%ylim(hAx(1), [-1 1])
%linkaxes([hAx(1) hAx(2)], 'x');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
xlim([1 45])

print(thisfig, sprintf('%sEXsubj%istim%i.pdf', allplotDir,subjectNum,stim), '-dpdf')


