function ds = PID(categSep,Kp,Ki,Kd,OptimalForget,maxIncrement)

%categSep is vector from 1:t with all the available categSep points
errorC = categSep-OptimalForget;
if length(errorC) == 1
    errorC = [0 errorC];
end
% use trapezoidal integration for error control
ds_temp = Kp*errorC(end) + Ki*(.5*errorC(1) + sum(errorC(2:end-1)) + .5*errorC(end)) + Kd*(errorC(end)-errorC(end-1));
ds = min([abs(maxIncrement) abs(ds_temp)]) * sign(ds_temp);

end