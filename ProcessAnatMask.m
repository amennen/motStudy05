projectName = 'motStudy03';
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
roi_name = 'paraHCG';
code_dir = ['/Data1/code/' projectName '/' 'code' '/']; %change to wherever code is stored
addpath(genpath(code_dir));
if IsLinux
    biac_dir = '/Data1/packages/BIAC_Matlab_R2014a/';
    bxhpath='/opt/BXH/1.11.1/bin/';
    fslpath='/opt/fsl/5.0.8/bin/';
end

%add necessary package
if ~exist('readmr','file')
    addpath(genpath(biac_dir));
    addpath([biac_dir '/mr/']);
    addpath([biac_dir '/general/'])
end
exfunc_reorient_fn = 'example_func_new_orientation';
highres2exfunc_mat = 'highres2example_func';

svec = [3 4 5 6 7 8 9 11 12 13];

for s = 1:length(svec)
    subjNum = svec(s);
    save_dir = ['/Data1/code/' projectName '/data/' num2str(subjNum) '/']; %this is where she sets the save directory!
    process_dir = [save_dir 'reg' '/'];
    roi_dir = ['/Data1/code/' projectName '/data/'];
    cd(process_dir)
    unix(sprintf('%sapplywarp -i %s%s.nii.gz -r %s.nii.gz -o %s_exfunc.nii.gz -w standard2highres_warp.nii.gz --postmat=%s.mat',fslpath,roi_dir,roi_name,exfunc_reorient_fn,roi_name,highres2exfunc_mat));
    
    %unzip, if necessary
    if exist(sprintf('%s_exfunc.nii.gz',roi_name),'file')
        unix(sprintf('gunzip %s_exfunc.nii.gz',roi_name));
    end
    
    %make bxh wrapper: do this to use for dicom files?
    unix(sprintf('%sbxhabsorb %s_exfunc.nii %s_exfunc.bxh',bxhpath,roi_name,roi_name));
    
    %load registered anatomical ROI
    maskStruct = readmr([roi_name '_exfunc.bxh'],'BXH',{[],[],[]});
    
    brainExtFunc = readmr([exfunc_reorient_fn '_mean_brain.bxh'], 'BXH',{[],[],[]});
    
    %load whole-brain mask-actual whole brain mask from example epi file:
    %see if masks are in the same space or not
    load(fullfile(process_dir,['mask_wholeBrain' '.mat']));
    
    %rotate anatomical ROI to be in the same space as the mask - check that this works for different scans/ROIs
    anatMaskRot = zeros(size(mask));
    brainExtRot = zeros(size(mask));
    for i = 1:size(maskStruct.data,3)
        anatMaskRot(:,:,i) = rot90(maskStruct.data(:,:,i)); %rotates entire slice by 90 degrees
        brainExtRot(:,:,i) = rot90(brainExtFunc.data(:,:,i));
    end
    
    %overwrite whole-brain mask
    mask = logical(anatMaskRot); %make it so it's just 1's and 0's
    brainExt = logical(brainExtRot);
    allinMask = find(anatMaskRot);
    allinBrainExt = find(brainExt);
    mask_indices = allinMask(find(ismember(allinMask,allinBrainExt))); %these are the good mask indices that are only brain
    [gX gY gZ] = ind2sub(size(mask),mask_indices);
    mask_brain = zeros(size(mask,1),size(mask,2),size(mask,3));
    for j=1:length(mask_indices)
        mask_brain(gX(j),gY(j),gZ(j)) = 1;
    end
    
    checkMask = 1;
    if checkMask
        plot3Dbrain(mask,[], 'mask')
        plot3Dbrain(mask_brain, [], 'mask_brain')
    end
    %save anatomical mask
    save(fullfile(process_dir, [roi_name, '_anat_mask']), 'mask_brain');
    fprintf('Done with mask creation\n');
    % if cd into the directory, cd out of it back to the general exp folder
    cd (code_dir)
end