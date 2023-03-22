% script to generate kspace trajectory and create psf for different acquisitions
%
% (1) reconstructs the psf (a decaying T2 signal)
% (2) reconstructs noise only
%
% tricky: code uses eval(expression) to change settings

% start fresh
clear g; 

%% important scanner settings (default values)
g.fov = [64 64 64]; % mm
g.siz = [64 64 64]; % matrix
g.lRadialSpokes = 40000; % no. spokes (full Nyquist ~ pi*g.siz(1)^2)
g.dGamma = 11.262; % gyromagnetic ratio (H=42.5775 NA=11.262)
g.uReadoutType = 0; % 0=radial 1=density adapted 2=corkscrew
g.lT0 = 0; % time to start density adapted (0=automatic)
g.bandwidth = 80; % Hz/pixel
g.osf = 1.5; % oversampling factor (1=Nyquist limit)
g.bRandomizePhase = 0; % randomize corkscrew phase
g.lAngularRate = 0; % corkscrew rotation rate (0=automatic)
g.n = 4; % tesselating polygon (3=triangle 4=square etc.)
g.spokeDirections = 'golden'; % 'golden' or 'spiral'

%% reconstruction settings (default values)
g.J = 4; % width of kaiser-bessel kernel
g.u = 2; % oversampling factor
g.radial = 0; % radial kernel (or separable)
g.maxit = 0; % no. of least squares iterations
g.damping = 0.0; % least squares Tikhonov damping
g.lambda = 0.0; % least squares phase constraint
g.density = ''; % '' for automatic or 'radial'
g.gpu = 1; % use GPU if available

%% evaluate "expression"

if exist('expression','var')
    eval(expression);
end

%% derived (or less important) settings
g.lInterleaves = 13; % reorder spokes into interleaves
g.dNyquist = 1e6 / g.bandwidth / g.siz(1); % us
g.lDwellTime = round(g.dNyquist / g.osf / 1e-3); % ns
g.lDwellTime = 64570
g.lGRADdelay = 0; % us gradient delay correction
g.lGRADsettle = 0; % us (settling time after ramp, 100us)
g.lReadRewTime = g.lGRADsettle; % us (time for rewinder, e.g. lGRADsettle) 

% heuristic angular rate calcuation (deg/dwell)
if g.lAngularRate == 0
    if g.lRadialSpokes/g.siz(1)^2 < 0.5
        g.lAngularRate = acsc(sqrt(pi/g.lRadialSpokes)*g.siz(1)-1) * 180 / pi;
    else
        g.lAngularRate = acsc(sqrt(2*pi)-1) * 180 / pi; % fix N/matrix^2 = 0.5
    end
    g.lAngularRate = round(g.lAngularRate); % in C++ code it is a long
end

% noise std dev (to give a ballpark SNR ~ 10)
sd = sqrt(g.bandwidth) / 20;

% allow flexibility for scalar fov and siz args
if numel(g.fov)==1; g.fov = repmat(g.fov,1,3); end
if numel(g.siz)==1; g.siz = repmat(g.siz,1,3); end

% reshape NA and T2 for vectorized code
NA = reshape(NA,[],1);
T2 = reshape(T2,[],1);

%% gradients (Prisma PERFORMANCE settings)

g.dMaxGrad = g.bandwidth * g.siz(1) / g.dGamma / g.fov(1); % mT/m
m_dMinRiseTime = 5.3; % from IDEA sys command (5.3 / 5.55 / 10 / 20)
g.dMaxSlew = 1e3 / m_dMinRiseTime;

if g.dMaxGrad > 37; warning('gradient overflow (%f)',g.dMaxGrad); end
if g.dMaxSlew > 188.7; warning('slew rate overflow (%f)',g.dMaxSlew); end

disp(g);

%% create the trajectory and do recon

g = setup(g);

% show end-points of the spokes (1 interleave only)
subplot(2,3,1)
k = 1:floor(g.lRadialSpokes/g.lInterleaves);
plot3(squeeze(g.om(1,end,k)),squeeze(g.om(2,end,k)),squeeze(g.om(3,end,k)),'.-');
mytitle = sprintf('1 interleave (%i spokes)',g.lRadialSpokes);
grid on; title(mytitle);
ylabel('kx'); xlabel('ky'); zlabel('kz');

% show a few spokes at the north pole
subplot(2,3,2);
k=3; while numel(k)<7; k(end+1) = k(end)+ceil(g.lRadialSpokes/g.lInterleaves); end
k = [1 k]; k(k>size(g.om,3)) = []; % include north pole and discard out of bounds
plot3(squeeze(g.om(1,1:1:end,k)),squeeze(g.om(2,1:1:end,k)),squeeze(g.om(3,1:1:end,k)),'LineWidth',1);
mytitle = sprintf('%i spokes',numel(k));
axis tight; grid on; title(mytitle); 
ylabel('kx','FontSize',12); xlabel('ky','FontSize',12); zlabel('kz','FontSize',12);

% show the readout gradient
subplot(2,3,3);
t = (0:size(g.G,2)-1) / 100; % time (ms)
h = plot(t,g.G);
hold on;
plot(t,vecnorm(g.G),'--black');
td = (0:size(g.om,2)-1) * g.lDwellTime * 1e-6; % ms
plot(td,interp1(t,g.G(1,:),td),'.','Color',h(1).Color);
plot(td,interp1(t,g.G(2,:),td),'.','Color',h(2).Color);
plot(td,interp1(t,g.G(3,:),td),'.','Color',h(3).Color);
hold off
mytitle = sprintf('Gradients (T_0 = %i us)',g.lT0);
title(mytitle); xlabel('time (ms)'); grid on;
xlim([0 max(t)]); ylabel('G (mT/m)');
legend({'Gx','Gy','Gz','||G||'});

%% reconstruction

t = te0 + (0:size(g.om,2)-1) * g.lDwellTime * 1e-6; % ms
data = sum(NA.*exp(-t./T2),1).'; % sum the components
data = repmat(data,1,size(g.om,3)); % all spokes

% reconstruct delta function
img = g.recon.iNUFT(data,g.maxit,g.damping,[],'phase-constraint',g.lambda);

% reconstruct noise
noise = complex(randn(size(data)),randn(size(data))) * sd;
noise = g.recon.iNUFT(noise,g.maxit,g.damping,[],'phase-constraint',g.lambda);

%% display the psf and artefact-free fov

psf = real(img);

% middle slice
mid = 1+g.siz(3)/2; 

% artefact-free fov (radial)
artefact_free_fov = sqrt(g.lRadialSpokes/pi);

% artefact_free_fov is the diameter but simulations show it is
% underestimated by a factor of ~2 so we use it as the radius
% below (similar to the 2x overestimate of Nyquist no. spokes)
radius = artefact_free_fov; 

subplot(2,3,4);
imagesc(log(abs(psf(:,:,mid))),[-12 0]);
rectangle('Position',[mid-radius,mid-radius,2*radius,2*radius],...
          'Curvature',[1 1],'EdgeColor','#0072BD','LineWidth',2.0);
title(sprintf('FOV %.1f (%i spokes)',artefact_free_fov,g.lRadialSpokes));

%% interpolate to give a smoother psf

osf = 10; % oversampling (interpolation)
npixels = 5; % +/- no. pixels to retain

psf = interpft(psf,osf*g.siz(1),1);
range = (1+g.siz(1)*osf/2)+(-npixels*osf:npixels*osf);
psf = psf(range,:,:);

psf = interpft(psf,osf*g.siz(2),2);
range = (1+g.siz(2)*osf/2)+(-npixels*osf:npixels*osf);
psf = psf(:,range,:);

psf = interpft(psf,osf*g.siz(3),3);
range = (1+g.siz(3)*osf/2)+(-npixels*osf:npixels*osf);
psf = psf(:,:,range);

%% remove scaling effects

% peak amplitude of the psf
peak = max(psf(:));

% normalize psf to 1 for display (FWHM unchanged)
psf = psf / peak;

% factor out ifft scaling 
noise = noise * sqrt(prod(g.siz));

% volume of the psf (sum all pixels)
vol = sum(psf(:)) / osf^3;

%% display the center pixels of the psf

% middle slice
mid = npixels*osf+1;

subplot(2,3,5);
try;surfl(psf(:,:,mid));end;view(45,30);shading interp;
xlabel('x (mm)'); ylabel('y (mm)'); axis tight;
title(sprintf('PSF: [volume %.3f mm^3]',vol));

% beautify axes
range = linspace(-npixels,npixels,2*npixels+1);
set(gca,'XTick',1:osf:2*osf*npixels+1)
set(gca,'XTickLabel',range)
set(gca,'YTick',1:osf:2*osf*npixels+1)
set(gca,'YTickLabel',range)
%axis([1.5-eps numel(range)-0.5+eps 1.5 numel(range)-0.5 -0.1 1.03])

%% estimate the 1/2 max location along each axis

range = linspace(-npixels,npixels,2*osf*npixels+1);

xvalues = squeeze(psf(:,mid,mid));
yvalues = squeeze(psf(mid,:,mid));
zvalues = squeeze(psf(mid,mid,:));

 % duration of the readout
g.Taq = max(t)-min(t);

% extend range for valid convolution
kmaxr = range; 
while numel(kmaxr)<2*numel(range)
    dkmaxr = abs(kmaxr(2)-kmaxr(1)); % same spacing
    kmaxr = [kmaxr(1)-dkmaxr kmaxr kmaxr(end)+dkmaxr];
end
kmaxr = pi * kmaxr; % pi * kmax * r from Bornert paper 

% theory (Bornert MRM 2006;55:1075-1082)
Psv = 3*(sin(kmaxr)-kmaxr.*cos(kmaxr))./(kmaxr).^3;
Pdec = 6*(T2./g.Taq).^3./(1+(kmaxr.*T2./g.Taq).^2).^2;

% correct 1/0 problems
Psv(kmaxr==0) = 1;
k = ~isfinite(T2);
Pdec(k,:) = 0; Pdec(k,ceil(numel(kmaxr)/2)) = 1; % delta function

% convolve (sum the components? doesn't exactly agree with simulation)
Ptot = sum(NA.*conv2(Pdec,Psv,'same'),1);

% crop the invalid range
while numel(Ptot)>numel(range)
    Ptot = Ptot(2:end-1);
end
tvalues = Ptot / max(Ptot);

% fzero finds the 1/2 width (x2 for full width) 
myfunc = @(pt,values)interp1(range,values,pt,'spline')-max(values)/2;

fwhm3d(1) = fzero(@(pt)myfunc(pt,xvalues),0.6) * 2;
fwhm3d(2) = fzero(@(pt)myfunc(pt,yvalues),0.6) * 2;
fwhm3d(3) = fzero(@(pt)myfunc(pt,zvalues),0.6) * 2;
fwhm3d(4) = fzero(@(pt)myfunc(pt,tvalues),0.6) * 2;

% plot the (normalized) psf in each dimension
subplot(2,3,6);
plot(range,[xvalues(:) yvalues(:)]/max(psf(:)));
hold on; plot(range,zvalues(:)/max(psf(:)),'--'); hold off;
mytitle = sprintf('\nFWHM: [%.2f  %.2f  %.2f  %.2f]',fwhm3d);
hold on; plot(range,tvalues,':'); hold off; grid on; axis tight;
title(mytitle); legend({'x','y','z','theory'}); xlabel('pixel (mm)')

%% look at noise at center of the volume
range = -npixels:npixels;
tmp = noise(1+g.siz(1)/2+range,1+g.siz(2)/2+range,1+g.siz(3)/2+range);
stddev = std(tmp(:));

% local noise std map
% noisemap = stdlocal(noise,[1 1 1]*numel(range));
% ims(log10(noisemap(:,:,1+g.siz(3)/2))); title('log noise std');

%% update the plot (do once only for speed)
drawnow;
