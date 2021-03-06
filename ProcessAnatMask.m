projectName = 'motStudy05';
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
%roi_name = 'paraHCG';
roi_name = 'paraHCG.nii.gz';
code_dir = pwd; %change to wherever code is stored
addpath(genpath(code_dir));

setenv('FSLOUTPUTTYPE','NIFTI_GZ');
%add necessary package

functionalFN = 'reference';
highres2exfunc_mat = 'highres2example_func';



svec = [1 3:6 8 9];

for s = 1:length(svec)
    subjNum = svec(s);
    save_dir = ['/Data1/code/' projectName '/data/' num2str(subjNum) '/']; %this is where she sets the save directory!
    process_dir = [save_dir 'reg' '/'];
    roi_dir = ['/Data1/code/' projectName '/data/'];
    cd(process_dir)
    unix(sprintf('%sapplywarp -i %s%s.nii.gz -r %s.nii -o %s_exfunc.nii.gz -w standard2highres_warp.nii.gz --postmat=%s.mat',fslpath,roi_dir,roi_name,functionalFN,roi_name,highres2exfunc_mat));
    
    %unzip, if necessary
    if exist(sprintf('%s_exfunc.nii.gz',roi_name),'file')
        unix(sprintf('gunzip %s_exfunc.nii.gz',roi_name));
    end
    
    
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
   
    %save anatomical mask
    save(fullfile(process_dir, [roi_name, '_mask']), 'mask_brain');
    fprintf('Done with mask creation\n');
    % if cd into the directory, cd out of it back to the general exp folder
    cd (code_dir)
end