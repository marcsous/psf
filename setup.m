%% initialize recon parameters (trajectory)
function g = setup(g)

% this code is slightly cut+paste from IDEA C++
% so apologies for the strangeness in places...

% extract some parameters from "g"
lRadialSpokes = g.lRadialSpokes;
lInterleaves = g.lInterleaves;
n = g.n;
fov = g.fov;
siz = g.siz;
lT0 = g.lT0;
osf = g.osf;
dGamma = g.dGamma;
dMaxGrad = g.dMaxGrad;
dNyquist = g.dNyquist;
dMaxSlew = g.dMaxSlew;
lDwellTime = g.lDwellTime;
lReadRewTime = g.lReadRewTime;
uReadoutType = g.uReadoutType;
lGRADdelay = g.lGRADdelay;
lAngularRate = g.lAngularRate;
lGRADsettle = g.lGRADsettle;
bRandomizePhase = g.bRandomizePhase;

% basic argument checks
if mod(lRadialSpokes,1); error('lRadialSpokes must be integer'); end
if mod(lInterleaves,1); error('lInterleaves must be integer'); end
if mod(lDwellTime,1); error('lDwellTime must be integer'); end

% anisotropic FOV requires more/less spokes than lRadialSpokes
dFOVAnisotropy = fov(1) / fov(end);
lLinesToMeasure = round(lRadialSpokes / dFOVAnisotropy); % estimate

% to help with cut and paste from C++ code
ny = lLinesToMeasure; M_PI = pi; smoothstep = @(x) x*x*(3-2*x);

% nominal increment along kz
dInc = 2 * dFOVAnisotropy / (lRadialSpokes-1);

% projection angles (golden angle)
for lI = 0:flintmax() % Inf
    dV = 1 - lI * dInc; % decrement from +1 to -1
    if dV < -1-dInc/2; break; end % allow a little overrun
    g.dPolarAngle(lI+1) = acos(max(dV,-1.0)); % stretch spacing along z
    g.dAzimuAngle(lI+1) = lI * M_PI * (3.0 - sqrt(5.0)); % golden angle
    
    if isequal(g.spokeDirections,'spiral')
        if lI == 0 || dV <= -1
            g.dAzimuAngle(lI+1) = 0;
        else
            g.dAzimuAngle(lI+1) = mod(g.dAzimuAngle(lI)+3.6/sqrt(lRadialSpokes*(1-dV*dV)),2*pi);
        end
    end
end
lLinesToMeasure = lI; ny = lLinesToMeasure; % this is now the correct no. spokes

% angles are saved in the MDH as uint16 (0-65535) so quantize
for lI = 0:lLinesToMeasure-1
    % put angles in 0-2pi range (no loss of precision)
    g.dPolarAngle(lI+1) = mod(g.dPolarAngle(lI+1), 2*M_PI);
    g.dAzimuAngle(lI+1) = mod(g.dAzimuAngle(lI+1), 2*M_PI);
    
    % quantize to "dconv" for perfect conversion into uint16
    dconv = (2*M_PI) / 65535;
    g.dPolarAngle(lI+1) = round(g.dPolarAngle(lI+1) / dconv) * dconv;
    g.dAzimuAngle(lI+1) = round(g.dAzimuAngle(lI+1) / dconv) * dconv;
end

% reorder spokes by interleaving (only useful for golden)
lI = 0;
for lInterleaf = 0:lInterleaves-1
    lJ = lInterleaf;
    while lJ < lLinesToMeasure
        lNewOrder(lI+1) = lJ+1; % MATLAB offsets
        lI = lI+1;
        lJ = lJ+lInterleaves;
    end
end
g.dAzimuAngle = g.dAzimuAngle(lNewOrder);
g.dPolarAngle = g.dPolarAngle(lNewOrder);

%% readout gradient

% target spatial frequency = dGamma (MHz/T) * integral{gradient (mT/m) * dt (s)} (1/mm)
kmax = dGamma * 0.5 * siz(1) * dMaxGrad * dNyquist * 1e-6; % 1/mm

% sanity check - scanner was crashing around here, probably due to nonsense values
if (kmax<=0 || dMaxGrad<=0 || dNyquist<=0 || dMaxSlew<=0)
    error('kmax=%f, dMaxGrad=%f, dNyquist=%f, dMaxSlew=%f',kmax,dMaxGrad,dNyquist,dMaxSlew);
end

% define the gradient shape: aGz=mT/m, dt=µs, slew=mT/m/ms
aGx = []; aGy = []; aGz = []; % arrays for Gx, Gy Gz waveforms
dt = 10;          % raster time (not tested except = 10) [µs]
k0 = kmax;        % spatial frequency to go density adapted
kz = 0.0;	      % current spatial frequency [phase rolls/mm]
Gz = 0.0;         % current value of the gradient [mT/m]
lGradSamples = 0; % no. points on the entire readout gradient
lRampTime = ceil(1e3 * dMaxGrad / dMaxSlew / dt) * dt; % us
lRampSamples = lRampTime/dt; % no. points on the gradient ramp

% readout rewinder = lGRADsettle on the flattop plus ramp up
if lReadRewTime
    kz = kz - dGamma * dMaxGrad * 1e-6 * (lGRADsettle + lRampTime / 2);
end

% ramp up to dMaxGrad
for lGradSamples = 1:lRampSamples
    
    aGz(lGradSamples+1) = Gz;
    kz = kz + dGamma * Gz * 1e-6 * dt;
    Gz = Gz + dMaxGrad / lRampSamples;
    
end

% density adapted readout
if uReadoutType~=0

    if lT0==0
        
        % optimal k0 
        k0 = sqrt(n*tan(M_PI/n)*lRadialSpokes/16/M_PI) / fov(1);
        
        % optimal lT0 based on k0 (us)
        lT0 = round(lRampTime + (k0-kz) / dGamma / dMaxGrad / 1e-6);
        
    end

    k0 = kz + max(lT0-lRampTime,0) * dGamma * Gz * 1e-6;

    % pass back to caller
    g.k0 = k0*fov(1);
    g.lT0 = lT0;

end

% FlatTop
while( kz < k0 )

    aGz(lGradSamples+1) = Gz;
    kz = kz + dGamma * Gz * 1e-6 * dt;
    
    lGradSamples = lGradSamples + 1;
    
end

% density adapted part
while( kz < kmax )
    
    % Nagel 2009 Eq 6-7
    tmp = 3 * dGamma * power(k0,2) * dMaxGrad  * (lGradSamples*dt-lT0) * 1e-6 + power(k0,3);
    
    if (tmp <= 0)
        error('something wrong with density adapted calculation');
    end
    Gz = power(k0,2) * dMaxGrad * power(tmp,-2.0/3.0);
    
    % slew rate protection as Gz drops from dMaxGrad
	if (aGz(lGradSamples)-Gz > dMaxSlew * dt * 1e-3)
        Gz = aGz(lGradSamples) - dMaxSlew * dt * 1e-3;
    end
    
    aGz(lGradSamples+1) = Gz;
    kz = kz + dGamma * Gz * 1e-6 * dt;
    
    lGradSamples = lGradSamples + 1;
    
end

% just repeat the last point - ramp down later
for lI = 0:lRampSamples
    aGz(lGradSamples+lI+1) = aGz(lGradSamples);
end

% check for any silliness
if lGradSamples > 50000
    warning('excessive no. points (%i) - probably an error.',lGradSamples)
    keyboard
end

%% non-projection options

% starting point in kz [1/mm]
if lReadRewTime==0
    kz = 0;  
else
    kz = -dGamma * dMaxGrad * 1e-6 * (lGRADsettle + lRampTime / 2);
end

% angle and radius of the oscillation
dAngle = 0; dRadius = 0;
 
% create corkscrew kx,ky on dt raster
for lI = 0:lGradSamples+lRampSamples
    
    % radial spacing along the spoke in time dt
    dkrad = dGamma * aGz(lI+1) * 1e-6 * dt; % 1/mm

    % increment kz
    kz = kz + dkrad; % 1/mm
  
    % corkscrew on the density adapted part
    if (kz > k0 && uReadoutType==2 && lAngularRate)
                
        % rotation increment in time dt
        dInc = lAngularRate * dt / dNyquist; % deg

        % increment the angle
        dAngle = dAngle + dInc * M_PI / 180; % rad
        
        % target radius (neglecting krad)
        dRadius = (kz/k0/2) / (1+sin(lAngularRate*M_PI/360));
        
        % target radius including krad (solve quadratic)
        %p(1) = 1 - sin(lAngularRate*M_PI/360)^2;
        %p(2) = -kz/k0;
        %p(3) = (kz/k0/2)^2 - (dkrad*fov(1)/2)^2;
        %dRadius = min(roots(p));
   
        % ease into the transition
        if dAngle < M_PI
            dRadius = dRadius * smoothstep(dAngle / M_PI);
        end

    end
    
    % create kx, ky (cycles / mm)
    aGx(lI+1) = dRadius * sin(dAngle) / fov(1);
    aGy(lI+1) = dRadius * cos(dAngle) / fov(1);
    
end

% convert trajectories to gradients
for lI = 1:lGradSamples+lRampSamples
    aGx(lI) = (aGx(lI+1)-aGx(lI)) / (dGamma * 1e-6 * dt); % mT/m
    aGy(lI) = (aGy(lI+1)-aGy(lI)) / (dGamma * 1e-6 * dt); % mT/m
end

% ramp down
for lI = 1:lRampSamples+1
    
    % first thing - stop increasing
    if (abs(aGx(lGradSamples+lI)) > abs(aGx(lGradSamples)))
        aGx(lGradSamples+lI) = aGx(lGradSamples);
    end
    if (abs(aGy(lGradSamples+lI)) > abs(aGy(lGradSamples)))
        aGy(lGradSamples+lI) = aGy(lGradSamples);
    end
    
    % now ramp down to zero
    ramp = double(lRampSamples-lI+1) / lRampSamples;
    
    aGz(lGradSamples+lI) = aGz(lGradSamples+lI) * ramp;
    aGx(lGradSamples+lI) = aGx(lGradSamples+lI) * smoothstep(ramp);
    aGy(lGradSamples+lI) = aGy(lGradSamples+lI) * smoothstep(ramp);
    
end

% check gradient and slew rate
dMaxGradXY = max(hypot(aGx,aGy));
dMaxSlewXYZ = vecnorm([diff(aGx);diff(aGy);diff(aGz)]) / dt / 1e-3;
dMaxSlewXYZ = max(dMaxSlewXYZ(lRampSamples+1:end-lRampSamples)); % skip the ramps
fprintf(' dMaxGradXY = %.2f dMaxSlewXYZ = %.2f (max = %.2f)\n',dMaxGradXY,dMaxSlewXYZ,dMaxSlew);

% save for display in calling function
g.G = [aGx;aGy;aGz];
g.dMaxGradXY = dMaxGradXY;
g.dMaxSlewXYZ = dMaxSlewXYZ;

%% create gradients for all projections

% estimate no. samples at 1 pt / dwell time
nx = ceil(lGradSamples * dt / lDwellTime / 1e-3);

% resample gradient to the sample times (pad ends to avoid out of bounds due to delay)
raster_time = dt * (-1:lGradSamples+lRampSamples+1); % us
sample_time = (lDwellTime * 1e-3) * (0:nx-1) - lGRADdelay; % us

aGx = interp1(raster_time,[0 aGx 0],sample_time);
aGy = interp1(raster_time,[0 aGy 0],sample_time);
aGz = interp1(raster_time,[0 aGz 0],sample_time);

% trajectory = dGamma (MHz/T) * integral[gradient (mT/m) * dt (ms)]
om = dGamma * cumsum([aGx;aGy;aGz] * lDwellTime * 1e-9, 2); % 1/mm

if lReadRewTime
    om(3,:) = om(3,:) - dGamma * dMaxGrad * (lGRADsettle + lRampTime / 2) * 1e-6; % 1/mm
end

% trim out of bounds points
nx = sum(om(3,:) <= kmax);
om = om(:,1:nx);

% rotate trajectory through all angles
g.om = zeros(3,nx,ny,'single');

for lLine = 1:ny
    
    % rotation matrix (same as on the scanner)
    aadRotationMatrix(1,1) = cos(g.dAzimuAngle(lLine))*cos(g.dPolarAngle(lLine));
    aadRotationMatrix(1,2) =-sin(g.dAzimuAngle(lLine));
    aadRotationMatrix(1,3) = cos(g.dAzimuAngle(lLine))*sin(g.dPolarAngle(lLine));
    aadRotationMatrix(2,1) = sin(g.dAzimuAngle(lLine))*cos(g.dPolarAngle(lLine));
    aadRotationMatrix(2,2) = cos(g.dAzimuAngle(lLine));
    aadRotationMatrix(2,3) = sin(g.dAzimuAngle(lLine))*sin(g.dPolarAngle(lLine));
    aadRotationMatrix(3,1) =-sin(g.dPolarAngle(lLine));
    aadRotationMatrix(3,2) = 0.0;
    aadRotationMatrix(3,3) = cos(g.dPolarAngle(lLine));
    
    % randomize starting phase
    if bRandomizePhase
        tmp = complex(om(1,:),om(2,:));
        tmp = tmp * exp(2*pi*i*rand(1));
        nom(1,:) = real(tmp);
        nom(2,:) = imag(tmp);
        nom(3,:) = om(3,:);
    else
        nom = om;
    end
    
    % apply rotation matrix
    g.om(:,:,lLine) = aadRotationMatrix * nom; % 1/mm

end

% change units for nufft (phase rolls/mm => phase rolls/fov)
g.om = g.om .* fov'; % phase rolls/fov (assuming fov in mm)

% check for any nasties - possibly interp1 out of bounds?
if any(~isfinite(g.om(:)))
    error('Non-finite values for om (probably from interp1). |GRADdelay| must be <=10µs')
end

%% radial density function for projection reconstruction

d = []; % automatic calculation

if isequal(g.density,'radial')
    if uReadoutType~=0
        error('density=''radial'' only works for uReadoutType=0');
    end
    k0 = sqrt(n*tan(M_PI/n)*lRadialSpokes/16/M_PI);
    d = min(k0^2,sum(g.om.^2));
    d = reshape(d,[],1) / k0^2;
end

%% create nufft object

g.recon = nufft_3d(g.om,g.siz,'J',g.J,'u',g.u,'radial',g.radial,'gpu',g.gpu,'d',d);
