projectName = 'motStudy05';
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
%roi_name = 'paraHCG';
code_dir = pwd; %change to wherever code is stored
addpath(genpath(code_dir));

setenv('FSLOUTPUTTYPE','NIFTI_GZ');
%add necessary package
fslpath = '/opt/pkg/FSL/fsl-5.0.9/';
functionalFN = 'reference';
highresFN = 'anat';
highres2exfunc_mat = 'highres2example_func';
roi_name = 'paraHCG.nii.gz';
proj_dir = '/jukebox/norman/amennen/motStudy05_transferred';
roi_dir = proj_dir;

svec = [1 3:6 8 9];

for s = 1:length(svec)
    subjNum = svec(s);
    subjectName = ['S' num2str(subjNum)];
    subj_dir = [proj_dir '/' subjectName];
    reg_dir = [subj_dir '/' 'analysis/firstlevel/localizer.feat/reg/'];
    anat_dir = [subj_dir '/' 'anat/anat/nii/'];
    func_dir = [subj_dir 
    %compute inverse transform (standard to highres)
    unix(sprintf('%sconvert_xfm -inverse -omat %sstandard2highres.mat %shighres2standard.mat', fslpath,reg_dir,reg_dir));
    unix(sprintf('%sinvwarp -w %shighres2standard_warp -o %sstandard2highres_warp -r %s%s_brain.nii.gz',fslpath,reg_dir,reg_dir,reg_dir,highresFN));
    
    unix(sprintf('%sapplywarp -i %s%s.nii.gz -r %s.nii -o %s_exfunc.nii.gz -w standard2highres_warp.nii.gz --postmat=%s.mat',fslpath,roi_dir,roi_name,functionalFN,roi_name,highres2exfunc_mat));
    
    %unzip, if necessary
    if exist(sprintf('%s_exfunc.nii.gz',roi_name),'file')
        unix(sprintf('gunzip %s_exfunc.nii.gz',roi_name));
    end
    
    
    niftiFN = sprintf('%s_exfunc.nii',roi_name);
    maskData = readnifti(niftiFN);
    unix(sprintf('%sbet %s.nii.gz %s_brain -R',fslpath,functionalFN_RE,functionalFN_RE)); % check that this is okay!
    
    %%
% now unzip and convert to load into matlab
%unzip, if necessary
if exist(sprintf('%s.nii.gz',functionalFN_RE),'file')
    unix(sprintf('gunzip %s.nii.gz',functionalFN_RE));
end
if exist(sprintf('%s_brain.nii.gz',functionalFN_RE),'file')
    unix(sprintf('gunzip %s_brain.nii.gz',functionalFN_RE));
end
    
    funcData = readnifti(sprintf('%s_brain.nii',functionalFN_RE));
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
    save(fullfile(process_dir, [roi_name, '_mask']), 'mask_brain');
    fprintf('Done with mask creation\n');
    % if cd into the directory, cd out of it back to the general exp folder
    cd (code_dir)
end