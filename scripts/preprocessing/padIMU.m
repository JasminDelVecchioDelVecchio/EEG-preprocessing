function imu_pad = padIMU(imu, lags, val)

mi = min(lags);
ma = max(lags);

[l m] = size(imu);
lpad = l+(ma-mi);

imu_pad = val*ones(lpad, m);
inds = (-mi+1:lpad-ma)';

imu_pad(round(inds), :) = imu;

