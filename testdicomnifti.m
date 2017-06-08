%testdicomnifti
dicomdir = '/Data1/subjects/20170120.0120171_motStudy03.0120171_motStudy03/';
scanNum =7;
exFuncScanNum = scanNum;
exFunc_scanstr = num2str(exFuncScanNum, '%2.2i');
exFunc_test_file = fullfile(dicomdir,['001_0000' exFunc_scanstr '_000008.dcm']);
brainmask = ones(64,64,36);
addpath(genpath('/Data1/code/motStudy03/code'))
roi = logical(brainmask);
image = dicomread(exFunc_test_file)
datafromdicom = ReadFile(exFunc_test_file,64,roi);

%% convert to nifti and then test
exfunc_str = exFunc_test_file; %general string for ALL mprage files**
t1 = GetSecs;
exffn = 'testingexfunc';
exfre = 'testingexfunc_reoriented';
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,exfunc_str,exffn));
%reorient bxh wrapper
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,exffn,exfre));
%convert the reoriented bxh wrapper to a nifti file
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,exfre,exfre))
t2 = GetSecs-t1;
addpath(genpath('/Data1/code/multibandutils/'))
niftidata = readnifti(exfre);
niftivec = niftidata(roi);

%% check with standard nifti
standard = '/Data1/code/motStudy03/data/MNI152_T1_2mm_brain.nii';
standardD = readnifti(standard);
%% now check with nifti after it's converted using dcm2niix
fancy = readnifti('test.nii');
figure;
imagesc(fancy(:,:,10))
dicomdir = '/Data1/subjects/20170120.0120171_motStudy03.0120171_motStudy03/'
fileOut = 'outputNifti'
process_dir = '/Data1/code/motStudy04/code';
dcm2path = '/opt/MRICROGL/2-2016/';
scanNum =7;
exFuncScanNum = scanNum;
exFunc_scanstr = num2str(exFuncScanNum, '%2.2i');
exFunc_test_file = fullfile(dicomdir,['001_0000' exFunc_scanstr '_000008.dcm']);

t1 = GetSecs;
unix(sprintf('%sdcm2niix %s -f %s -o %s -s y %s',dcm2path,dicomdir,fileOut,process_dir,exFunc_test_file))
t2 = GetSecs;
conv = t2-t1;
N1 = readnifti('outputNifti.nii');


% now test and make sure converts the same!
t1= GetSecs;
unix(sprintf('%sdicom2bxh %s %s.bxh',bxhpath,exFunc_test_file,exffn));
unix(sprintf('%sbxhreorient --orientation=LAS %s.bxh %s.bxh',bxhpath,exffn,exfre));
unix(sprintf('%sbxh2analyze --overwrite --analyzetypes --niigz --niftihdr -s %s.bxh %s',bxhpath,exfre,exfre))
t2 = GetSecs;
conv = t2-t1;

N2 = readnifti('testingexfunc_reoriented.nii.gz');
% same yay!
% next update the scripts to convert the dicom to nifti in this ways
% then test with mcflirt