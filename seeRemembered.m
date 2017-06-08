svec = [8 12 14 15 16 18 20:22 26 27 28 29 30];
RT = [8 12 14 15 18 21 22];
YC = [16 20 26 27 28 29 30];
RT_m = [8 12 14 15 18 21 22];
YC_m = [16 28 20 26 27 29 30];
iRT = find(ismember(svec,RT));
iYC = find(ismember(svec,YC));
iRT_m = find(ismember(svec,RT_m));
for i = 1:length(YC_m)
    iYC_m(i) = find(svec==YC_m(i));
end


for i = 1:length(svec)
    n_rem(i) = length(findRememberedStim(svec(i)));
    remembered{i} = findRememberedStim(svec(i));
end

for i = 1:length(iYC_m)
    overlapping{i} = intersect(remembered{iRT_m(i)},remembered{iYC_m(i)});
end