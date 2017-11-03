function ProcessAnatMask_spock(s)
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
%roi_name = 'paraHCG';

fslpath = '/opt/pkg/FSL/fsl-5.0.9/';
functionalFN = 'reference';
highresFN = 'anat';
highres2exfunc_mat = 'highres2example_func';
roi_name = 'paraHCG';
proj_dir = '/jukebox/norman/amennen/motStudy05_transferred';
roi_dir = proj_dir;
svec = [1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40];
subjNum = svec(s);
subjectName = ['S' num2str(subjNum)];
subj_dir = [proj_dir '/' subjectName];
reg_dir = [subj_dir '/' 'analysis/firstlevel/localizer.feat/reg/'];
anat_dir = [subj_dir '/' 'anat/anat/nii/'];
ref_dir = [subj_dir '/' 'func/'];
%compute inverse transform (standard to highres)
unix(sprintf('%sconvert_xfm -inverse -omat %sstandard2highres.mat %shighres2standard.mat', fslpath,reg_dir,reg_dir));
unix(sprintf('%sinvwarp -w %shighres2standard_warp -o %sstandard2highres_warp -r %s%s_brain.nii.gz',fslpath,reg_dir,reg_dir,reg_dir,highresFN));

unix(sprintf('%sapplywarp -i %s/%s.nii.gz -r %s%s.nii -o %s%s_func.nii.gz -w %sstandard2highres_warp.nii.gz --postmat=%shighres2example_func.mat',fslpath,roi_dir,roi_name,ref_dir,functionalFN,ref_dir,roi_name,reg_dir,reg_dir));

%unzip, if necessary
if exist(sprintf('%s%s_func.nii.gz',ref_dir,roi_name),'file')
    unix(sprintf('gunzip %s%s_func.nii.gz',ref_dir,roi_name));
end


niftiFN = sprintf('%s%s_func.nii',ref_dir,roi_name);
maskData = readnifti(niftiFN);
unix(sprintf('%sbet %s%s.nii.gz %s%s_brain -R',fslpath,ref_dir,functionalFN,ref_dir,functionalFN)); % check that this is okay!

%%
% now unzip and convert to load into matlab
%unzip, if necessary
if exist(sprintf('%s%s.nii.gz',ref_dir,functionalFN),'file')
    unix(sprintf('gunzip %s%s.nii.gz',ref_dir,functionalFN));
end
if exist(sprintf('%s%s_brain.nii.gz',ref_dir,functionalFN),'file')
    unix(sprintf('gunzip %s%s_brain.nii.gz',ref_dir,functionalFN));
end

funcData = readnifti(sprintf('%s%s_brain.nii',ref_dir,functionalFN));
%%
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

%save anatomical mask
save(fullfile(ref_dir, [roi_name, '_mask']), 'mask_brain');
fprintf('Done with mask creation\n');
% if cd into the directory, cd out of it back to the general exp folder
end