% want to plot the dot speed and category separation timecourse
% need: speed (from behavioral file)
% category separation (can also pull from behavioral file)
close all;
clear all;
projectName = 'motStudy04';
allspeeds = [];
allsep = [];
nstim = 10;
nTRs = 25;
nblock = 3;
svec = [3];

nsub = length(svec);
sepbystim = zeros(nstim,nTRs*nblock);
speedbystim = zeros(nstim,nTRs*nblock);
MOT_PREP = 5;
colors = [207 127 102;130 161 171; 207 64 19]/255;

%colors = [110 62 106;83 200 212; 187 124 181]/255;
plotstim = 1; %if you want trial by trial plots
allplotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' ];

for s = 1:nsub
    subjectNum = svec(s);
    allsep = [];
    fbsep = [];
    allspeeds = [];
    for iblock = 1:nblock
        blockNum = iblock;
        SESSION = 20 + blockNum;
        behavioral_dir = [fileparts(which('mot_realtime04MB.m')) '/BehavioralData/' num2str(subjectNum) '/'];
        save_dir = ['/Data1/code/' projectName '/data/' num2str(subjectNum) '/'];
        classOutputDir = fullfile(save_dir,['motRun' num2str(blockNum)], 'classOutput/');
        runHeader = fullfile(save_dir,[ 'motRun' num2str(blockNum) '/']);
        fileSpeed = dir(fullfile(behavioral_dir, ['mot_realtime04_' num2str(subjectNum) '_' num2str(SESSION)  '*.mat']));
        d = load(fullfile(behavioral_dir, fileSpeed(end).name));
        %get hard speed
        prep = dir([behavioral_dir 'mot_realtime04_' num2str(subjectNum) '_' num2str(MOT_PREP)  '*.mat']);
        prepfile = [behavioral_dir prep(end).name];
        lastRun = load(prepfile);
        hardSpeed(s) = 30 - lastRun.stim.tGuess(end);
        
        plotDir = ['/Data1/code/' projectName '/' 'Plots2' '/' num2str(subjectNum) '/'];
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        
        % get the speed for every trial in that block
        allEv = zeros(nstim,25);
        allSpeed = zeros(nstim,25);
        for i = 1:nstim
            allEv(i,2:end) = d.rtData.RTVEC{i}(1:nTRs-1);
            allSpeed(i,:) = d.rtData.allSpeeds{i}(1:nTRs);
        end

        [~,indSort] = sort(d.stim.id);
        sepinorder = allEv(indSort,:);
        speedinorder = allSpeed(indSort,:);
        
       
        sepbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = sepinorder;
        speedbystim(:,(iblock-1)*nTRs + 1: iblock*nTRs ) = speedinorder;
        
    end
    
    if plotstim
        for stim = 1:nstim
            thisfig = figure(stim*50);
            clf;
            x = 1:nTRs*nblock;
            [hAx,hLine1, hLine2] = plotyy(x,sepbystim(stim,:),x,speedbystim(stim,:));
            xlabel('TR Number (2s)')
            ylabel(hAx(2), 'Dot Speed', 'Color', 'k')
            ylabel(hAx(1), 'Category Evidence', 'Color', 'k')
            ylim(hAx(2),[-20 20])
            ylim(hAx(1), [-1 1])
            xlim([0.5 80])
            set(hLine2, 'LineStyle', '-', 'Color', colors(2,:), 'LineWidth', 5)
            set(hLine1, 'LineStyle', '-', 'Color', colors(1,:), 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 25)
            linkaxes([hAx(1) hAx(2)], 'x');
            title(sprintf('Subject: %i Stimulus ID: %i',subjectNum,stim));
            set(findall(gcf,'-property','FontSize'),'FontSize',20)
            set(findall(gcf,'-property','FontColor'),'FontColor','k')
            set(hAx(1), 'FontSize', 12)
            set(hAx(2), 'YColor', colors(2,:), 'FontSize', 16, 'YTick', [-20:5:5]); %'YTickLabel', {'0', '1', '2', '3', '4', '5})
            set(hAx(1), 'YColor', colors(1,:), 'FontSize', 16, 'YTick', [-1:.5:1], 'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'});
            hold on;
            legend('Ev', 'Dot Speed')
            for rep = 1:2
                line([rep*nTRs+.5 rep*nTRs + .5], [-10 15], 'color', 'k', 'LineWidth', 2);
            end
            line([0 80], [0.15 0.15], 'color', [140 136 141]/255, 'LineWidth', 2.5,'LineStyle', '--');
            %         savefig(sprintf('%sstim%i.fig', plotDir,stim));
            %print(thisfig, sprintf('%sstim%i.pdf', plotDir,stim), '-dpdf')
        end
    end
end

