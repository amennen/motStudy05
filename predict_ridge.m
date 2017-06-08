function [labelsPredicted acts] = predict_ridge(data,model)

% data is the matrix for that single TR and nvoxels
% model is the scratchpad from classifier trainig
nConds = size(model.ridge.betas,2); %check this is right
for c = 1:nConds
    % get w vector 
    w = model.ridge.betas(:,c);
    acts(c) = w'*data';
end

%right now make threshold 0
[temp labelsPredicted] = max(acts);
%probability  = acts;

end