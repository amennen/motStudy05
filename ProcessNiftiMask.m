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
multipath = '/Data1/code/multibandutils/';
addpath(genpath(multipath));
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

% inputs (eventually function)
subjNum = 8;
%subjDate = '4-5-17';
subjDate = NaN;
runNum = 1;
highresScan = 2;
APScan = 3;
PAScan = 4;
functionalScan = 5;

% set up this subject's paths
projectName = 'motStudy04';
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

%% Process t1-weighted MPRAGE
highresFN = 'highres';
highresFN_RE = 'highres_re';
highresfiles_genstr = sprintf('%s001_0000%s_0*',dicom_dir,num2str(highresScan,'%2.2i')); %general string for ALL mprage files**
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,highresfiles_genstr,highresFN));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,highresFN,highresFN_RE));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --2niftihdr -s %s.bxh %s',bxhpath,highresFN_RE,highresFN_RE))
unix(sprintf('%sbet %s.nii.gz %s_brain.nii.gz -R',fslpath,highresFN_RE,highresFN_RE)) 
% for dcm2niix the command would be 'dcm2niix dicomdir -f test -o dicomdir -s y dicomdir/001_000007_000008.dcm'

% Register to standard=
unix(sprintf('%sflirt -in %s_brain.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -out highres2standard -omat highres2standard.mat -cost corratio -dof 12 -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -interp trilinear',fslpath,highresFN_RE));
unix(sprintf('%sfnirt --iout=highres2standard_head --in=%s.nii.gz --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2highres_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10', fslpath,highresFN_RE));
unix(sprintf('%sapplywarp -i %s_brain.nii.gz -r $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -o highres2standard -w highres2standard_warp',fslpath,highresFN_RE));
%compute inverse transform (standard to highres)
unix(sprintf('%sconvert_xfm -inverse -omat standard2highres.mat highres2standard.mat', fslpath));
unix(sprintf('%sinvwarp -w highres2standard_warp -o standard2highres_warp -r %s_brain.nii.gz',fslpath,highresFN_RE));

%% Process each spin echo file
onlyFirstTR = 1;
if onlyFirstTR % only use first TR for each spin echo to save time
    APname = 'SE_AP1';
    AP_re = [APname '_re'];
    AP_genstr = sprintf('%s001_0000%s_0000{01..48}.dcm',dicom_dir,num2str(APScan,'%2.2i'));
    unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,AP_genstr,APname));
    unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,APname,AP_re));
    unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,AP_re,AP_re))
    
    % now do the same thing with PA
    PAname = 'SE_PAB';
    PA_re = [PAname '_re'];
    PA_genstr = sprintf('%s001_0000%s_0000{01..48}.dcm',dicom_dir,num2str(PAScan,'%2.2i'));
    unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,PA_genstr,PAname));
    unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,PAname,PA_re));
    unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,PA_re,PA_re))
    
    % use topup to calculate differences
    % now combine them into a single image
    time1 = GetSecs;
    fieldmapfn = 'all_SE';
    unix(sprintf('%sfslmerge -t %s.nii.gz %s.nii.gz %s.nii.gz', fslpath,fieldmapfn,AP_re,PA_re))
    
    % now run topup!
    textfile = 'acqparamsONE.txt';
    cnffile = 'b02b0.cnf';
    unix(sprintf('%stopup --imain=%s.nii.gz --datain=%s%s --config=%s%s --out=topup_output --iout=topup_iout --fout=topup_fout --logout=topup_logout',fslpath,fieldmapfn,multipath,textfile,multipath,cnffile))
    
    % create magnitude image from topup
    unix(sprintf('%sfslmaths topup_iout -Tmean magnitude',fslpath))
    % create brain-extracted magnitude image
    unix(sprintf('%sbet magnitude magnitude_brain -R',fslpath)) % check that this is okay afterwards!!!!
    unix(sprintf('%sfslmaths topup_fout.nii.gz -mul 6.28 fieldmap_rads',fslpath))
    time2 = GetSecs;
    topuptime = time2-time1;
else
    APname = 'SE_AP3';
    AP_re = [APname '_re'];
    AP_genstr = sprintf('%s001_00000%s_0*',dicom_dir,num2str(APScan));
    unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,AP_genstr,APname));
    unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,APname,AP_re));
    unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,AP_re,AP_re))
    
    % now do the same thing with PA
    PAname = 'SE_PA3';
    PA_re = [PAname '_re'];
    PA_genstr = sprintf('%s001_00000%s_0*',dicom_dir,num2str(PAScan));
    unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,PA_genstr,PAname));
    unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,PAname,PA_re));
    unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,PA_re,PA_re))
    
    % use topup to calculate differences
    % now combine them into a single image
    time1 = GetSecs;
    fieldmapfn = 'all_SE';
    unix(sprintf('%sfslmerge -t %s.nii.gz %s.nii.gz %s.nii.gz', fslpath,fieldmapfn,AP_re,PA_re))
    
    % now run topup!
    textfile = 'acqparams.txt';
    cnffile = 'b02b0.cnf';
    unix(sprintf('%stopup --imain=%s.nii.gz --datain=%s%s --config=%s%s --out=topup_output --iout=topup_iout --fout=topup_fout --logout=topup_logout',fslpath,fieldmapfn,multipath,textfile,multipath,cnffile))
    
    % create magnitude image from topup
    unix(sprintf('%sfslmaths topup_iout -Tmean magnitude3',fslpath))
    % create brain-extracted magnitude image
    unix(sprintf('%sbet magnitude3 magnitude3_brain -r 100',fslpath)) % check that this is okay afterwards!!!!
    unix(sprintf('%sfslmaths topup_fout.nii.gz -mul 6.28 fieldmap_rads3',fslpath))
    time2 = GetSecs;
    topuptime = time2-time1;
end
%% Process example epi file
fileN = 8; % we can choose 10 later
functionalFN = 'exfunc';
functionalFN_RE = 'exfunc_re';
exfunc_str = sprintf('%s001_0000%s_0000%s.dcm',dicom_dir,num2str(functionalScan,'%2.2i'),num2str(fileN,'%2.2i')); %general string for ALL mprage files**
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,exfunc_str,functionalFN));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,functionalFN,functionalFN_RE));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,functionalFN_RE,functionalFN_RE))

% now register to highres!
t1 = GetSecs;
exfunc2highres_mat='example_func2highresTESTSE1';
highres2exfunc_mat='highres2example_funcTESTSE1';
unix(sprintf('%sepi_reg --epi=%s --t1=%s --t1brain=%s_brain --out=%s',fslpath,functionalFN_RE,highresFN_RE,highresFN_RE,exfunc2highres_mat))
%exfunc2highres_mat='example_func2highres';
%highres2exfunc_mat='highres2example_func';
unix(sprintf('%sepi_reg --epi=%s --t1=%s --t1brain=%s_brain --out=%sNOFIELDMAP',fslpath,functionalFN_RE,highresFN_RE,highresFN_RE,exfunc2highres_mat))

%unix(sprintf('%sepi_reg --epi=%s.nii.gz --t1=%s.nii.gz --t1brain=%s_brain.nii.gz --out=%s3 --fmap=fieldmap_rads3 --fmapmag=magnitude3 --fmapmagbrain=magnitude3_brain --echospacing=0.00035 --pedir=y',fslpath,functionalFN_RE,highresFN_RE,highresFN_RE,exfunc2highres_mat))
unix(sprintf('%sepi_reg --epi=%s.nii.gz --t1=%s.nii.gz --t1brain=%s_brain.nii.gz --out=%s --fmap=fieldmap_rads --fmapmag=magnitude --fmapmagbrain=magnitude_brain --echospacing=0.00035 --pedir=y',fslpath,functionalFN_RE,highresFN_RE,highresFN_RE,exfunc2highres_mat))

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
% % test that it's the same thing just not rotated as before
% oldm = load('/Data1/code/motStudy03/data/3/reg/retrieval_anat_mask_orig.mat')
% oldm = oldm.mask_brain;
% figure;
% imagesc(oldm(:,:,10))
% MaskRot = zeros(size(mask_brain));
% for i = 1:size(mask_brain,3)
%     MaskRot(:,:,i) = rot90(mask_brain(:,:,i)); %rotates entire slice by 90 degrees
% end
% sum(find(MaskRot==oldm))
%% Now test with other data
scanNum = 11;
thisTR = 40;
[patterns.fileAvail filename] = GetSpecificFMRIFile(dicom_dir,scanNum,thisTR);
% now convert to nifti
niftiname = sprintf('nifti%3.3i', thisTR);
unix(sprintf('%sdicom2bxh %s%s %s.bxh',bxhpath,dicom_dir,filename,niftiname));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,niftiname,niftiname));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,niftiname,niftiname))
% now mcflirt!
unix(sprintf('%smcflirt -in %s.nii.gz -reffile exfunc_re',fslpath,niftiname))
