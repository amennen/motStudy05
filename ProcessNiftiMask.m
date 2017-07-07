% ProcessNiftiMask
% Written by ACM April 2017

% Collect:
% - t1-weighted mp-rage
% - spin echo series
% - example functional scan
% Complete:
% 1. Register t1 to standard space
% 2. Combine spin echos with topup (could take a longer amount of time)
% 3. Register epi example to t1, including field map corrections
% 4. Calculate inverse transformation matrices
% 5. Apply to anatomical mask

% setup paths
biac_dir = '/Data1/packages/BIAC_Matlab_R2014a/';
bxhpath='/opt/BXH/1.11.1/bin/';
fslpath='/opt/fsl/5.0.9/bin/';
%add necessary package
if ~exist('readmr','file')
    addpath(genpath(biac_dir));
    addpath([biac_dir '/mr/']);
    addpath([biac_dir '/general/'])
end
multipath = '/Data1/code/multibandutils/';1
addpath(genpath(multipath));
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

% inputs (eventually function)


subjNum = 10;
%subjDate = '4-5-17';
subjDate = NaN;
runNum = 1;
highresScan = 5;
functionalScan = 6;

% set up this subject's paths
projectName = 'motStudy05';
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
save_dir = ['/Data1/code/' projectName '/data/' num2str(subjNum) '/']; %this is where she sets the save directory!
process_dir = [save_dir 'reg' '/'];
roi_dir = ['/Data1/code/' projectName '/data/'];
code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
if ~isnan(subjDate)
    subjDate = '5-10-17';
    subjectName = [datestr(subjDate,5) datestr(subjDate,7) datestr(subjDate,11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(subjDate,10) datestr(subjDate,5) datestr(subjDate,7) '.' subjectName '.' subjectName '/'];
else
    subjectName = [datestr(now,5) datestr(now,7) datestr(now,11) num2str(runNum) '_' projectName];
    dicom_dir = ['/Data1/subjects/' datestr(now,10) datestr(now,5) datestr(now,7) '.' subjectName '.' subjectName '/'];
end
addpath(genpath(code_dir));
if ~exist(process_dir)
    mkdir(process_dir)
end
cd(process_dir)

%% Process t1-weighted MPRAGE and check brain extraction!
highresFN = 'highres';
highresFN_RE = 'highres_re';
highresfiles_genstr = sprintf('%s001_0000%s_0*',dicom_dir,num2str(highresScan,'%2.2i')); %general string for ALL mprage files**
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,highresfiles_genstr,highresFN));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,highresFN,highresFN_RE));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --2niftihdr -s %s.bxh %s',bxhpath,highresFN_RE,highresFN_RE))
unix(sprintf('%sbet %s.nii.gz %s_brain.nii.gz -R',fslpath,highresFN_RE,highresFN_RE)) 
%unix(sprintf('%sbet %s.nii.gz %s_brain.nii.gz -r 90 -R',fslpath,highresFN_RE,highresFN_RE)) 

% for dcm2niix the command would be 'dcm2niix dicomdir -f test -o dicomdir -s y dicomdir/001_000007_000008.dcm'
fprintf('%sfslview %s.nii.gz\n',fslpath,highresFN_RE)
fprintf('%sfslview %s_brain.nii.gz', fslpath,highresFN_RE)

%% Register standard to nifti
% Register to standard=2
unix(sprintf('%sflirt -in %s_brain.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -out highres2standard -omat highres2standard.mat -cost corratio -dof 12 -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -interp trilinear',fslpath,highresFN_RE));
unix(sprintf('%sfnirt --iout=highres2standard_head --in=%s.nii.gz --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2highres_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10', fslpath,highresFN_RE));
unix(sprintf('%sapplywarp -i %s_brain.nii.gz -r $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -o highres2standard -w highres2standard_warp',fslpath,highresFN_RE));
%compute inverse transform (standard to highres)
unix(sprintf('%sconvert_xfm -inverse -omat standard2highres.mat highres2standard.mat', fslpath));
unix(sprintf('%sinvwarp -w highres2standard_warp -o standard2highres_warp -r %s_brain.nii.gz',fslpath,highresFN_RE));
%% Process example epi file
fileN = 6; % we can choose 10 later
functionalFN = 'exfunc';
functionalFN_RE = 'exfunc_re';
exfunc_str = sprintf('%s001_0000%s_0000%s.dcm',dicom_dir,num2str(functionalScan,'%2.2i'),num2str(fileN,'%2.2i')); %general string for ALL mprage files**
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,exfunc_str,functionalFN));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,functionalFN,functionalFN_RE));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,functionalFN_RE,functionalFN_RE))

% now register to highres!
t1 = GetSecs;
exfunc2highres_mat='example_func2highres';
highres2exfunc_mat='highres2example_func';
unix(sprintf('%sepi_reg --epi=%s --t1=%s --t1brain=%s_brain --out=%s',fslpath,functionalFN_RE,highresFN_RE,highresFN_RE,exfunc2highres_mat))
timefunc2highres = GetSecs-t1;
unix(sprintf('%sconvert_xfm -inverse -omat %s.mat %s.mat',fslpath,highres2exfunc_mat,exfunc2highres_mat));

% now register mask to all data
roi_name = 'retrieval';
unix(sprintf('%sapplywarp -i %s%s.nii.gz -r %s.nii.gz -o %s_exfunc.nii.gz -w standard2highres_warp.nii.gz --postmat=%s.mat',fslpath,roi_dir,roi_name,functionalFN_RE,roi_name,highres2exfunc_mat));
% check after here that the applied warp is binary and in the right
% orientation so we could just apply to nifti files afterwards
if exist(sprintf('%s_exfunc.nii.gz',roi_name),'file')
    unix(sprintf('gunzip %s_exfunc.nii.gz',roi_name));
end

% brain extract functional scan to make sure we stay inside the brain of
% the subject!
unix(sprintf('%sbet %s.nii.gz %s_brain -R',fslpath,functionalFN_RE,functionalFN_RE)); % check that this is okay!
%CHECK OKAY
fprintf('%sfslview %s_brain.nii.gz', fslpath,functionalFN_RE)

%%
% now unzip and convert to load into matlab
%unzip, if necessary
if exist(sprintf('%s.nii.gz',functionalFN_RE),'file')
    unix(sprintf('gunzip %s.nii.gz',functionalFN_RE));
end
if exist(sprintf('%s_brain.nii.gz',functionalFN_RE),'file')
    unix(sprintf('gunzip %s_brain.nii.gz',functionalFN_RE));
end

%% Now create mask
niftiFN = sprintf('%s_exfunc.nii',roi_name);
maskData = readnifti(niftiFN);
funcData = readnifti(sprintf('%s_brain.nii',functionalFN_RE));

mask = zeros(size(maskData));
maskLogical = logical(maskData);
brainLogical = logical(funcData);
allinMask = find(maskLogical);
allinBrain = find(brainLogical);
mask_indices = allinMask(find(ismember(allinMask,allinBrain))); %these are the good mask indices that are only brain
[gX gY gZ] = ind2sub(size(mask),mask_indices);
mask_brain = zeros(size(mask,1),size(mask,2),size(mask,3));
for j=1:length(mask_indices)
    mask_brain(gX(j),gY(j),gZ(j)) = 1;
end

% now save mask
save(fullfile(process_dir, [roi_name '_mask']),'mask_brain')
sprintf('done')