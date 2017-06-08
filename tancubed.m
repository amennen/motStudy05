function ds = tancubed(categSep,Scale,OptimalForget,maxIncrement)
if categSep <-0.64
    categSep = -0.64;
elseif categSep > 0.84
    categSep = 0.84;
end
ds_temp = Scale * tan(((categSep- OptimalForget)*pi/2).^3);
ds = min([abs(maxIncrement) abs(ds_temp)]) * sign(ds_temp);
end