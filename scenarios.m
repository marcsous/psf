% test various parameters in radial, density adapted and corkscrew sampling.
%
% the first set of simulations neglect T2, the second set use biexponential.
%
% the code keeps everything in the workspace to pass between scripts
%  - scenarios.m: sets up a few key situations to examine
%  - setup.m: creates trajectories and technical settings
%  - simulate.m: executes the reconstruction and displays
%
% tricky: simulate.m uses eval(expression) to override default settings.
%
% the total runtime is ~2 hours on a recent 8Gb GPU or ~24 hours on CPU.
%
clear all
close all

% need nufft_3d (https://github.com/marcsous/nufft_3d)
if ~exist('nufft_3d.m','file'); error('nufft_3d.m required'); end

% echo time (start of data acquisition)
te0 = 0; % ms

% reconstruction (1=regridding 10=least squares)
maxit = 10; 

% open all windows to prevent them popping up randomly
for k = 1:9; h=figure(k); h.Position=[8 604 424 363]; end; drawnow;

%% 1: single component long T2 - "ideal" case

NA = 1; T2 = Inf;

%% voxel fwhm vs nspokes
h=figure(1);h.Name='fwhm vs. nspokes';drawnow

myfov = 128;

% range of spokes to simulate
mynyquist = pi*myfov^2;
myspokes = round(linspace(mynyquist*0.06,mynyquist*0.96,25));

% 1=radial 2=density, 3=corksrew
for myl = 1:3
    
    for myk = 1:size(myspokes,2)

        expression = ['g.lRadialSpokes = myspokes(myk); g.fov = myfov; g.siz = g.fov;'];
        expression = [expression 'g.uReadoutType = myl-1; g.maxit = maxit;'];

        simulate;
        
        myvol(myk,myl) = vol;
        mystd(myk,myl) = stddev;
        myfwhm(myk,myl) = fwhm3d(1);
        
    end
    
end

figure(1); subplot(1,1,1);
h = plot(myspokes'./myfov.^2,myfwhm');
xlabel('N / matrix^2'); ylabel('FWHM (Δx)'); grid on;
h(1).Marker = 'o'; h(2).Marker = 'sq'; h(3).Marker = 'd';
h(1).MarkerIndices = 1:2:numel(myspokes);
h(2).MarkerIndices = 2:2:numel(myspokes);
legend({'radial', 'density adpt','corkscrew'}); xlim([0.1 3.1]);
drawnow;

% save for Figure 6
save_myfov_no_T2 = myfov;
save_myfwhm_no_T2 = myfwhm; 
save_myspokes_no_T2 = myspokes;

clear my*

%% fwhm vs corkscrew angular rate
h=figure(2);h.Name='fwhm vs corkscrew angular rate';drawnow

myfov = 128;
myangles = 0.5:0.5:95;
myspokes = round([0.1 0.25 0.5] * myfov^2);

for myj = 1:numel(myspokes)
    
    for myk = 1:numel(myangles)
        
        expression = ['g.lRadialSpokes = myspokes(myj); g.fov = myfov; g.siz = g.fov;'];
        expression = [expression 'g.uReadoutType = 2; g.lAngularRate = myangles(myk);'];
        expression = [expression 'g.maxit = maxit; g.osf = 1;']; % osf=1 shows peaks better 
        
        simulate;
        
        myvol(myj,myk) = vol;
        mystd(myj,myk) = stddev;
        myfwhm(myj,myk) = fwhm3d(1);
        mymaxslew(myj,myk) = g.dMaxSlewXYZ;
    
        subplot(2,3,3);
        plot(myangles(1:myk),myfwhm(1:myj,1:myk)');
        title('Effect of ϕ versus no. spokes');
        xlabel('Angular Rate (deg/dwelltime)'); ylabel('FWHM (Δx)');
        legend(num2str(myspokes(1:myj)'/myfov^2,'%.2f matrix^2'));
        xlim([0 max(myangles)]); drawnow;
        
    end
    
end

subplot(1,1,1);
[h,ax1,ax2] = plotyy(myangles,myfwhm',myangles,min(mean(mymaxslew),100));
set(ax2,'linestyle','--','color','black');
xlabel('Angular Rate (deg/dwelltime)');
ylabel(h(1),'FWHM (Δx)');
ylabel(h(2),'Slew Rate (% of max)','color','black'); 
set(h(2),{'ycolor'},{'black'})
legend(num2str(myspokes'/myfov^2,'%.2f matrix^2'),'Location','NorthWest');
xlim(h(1),[0 95]);
xlim(h(2),[0 95]);
xticks(h(1),0:10:90)
xticks(h(2),0:10:90)
ylim(h(1),[1.5 2.2]);
ylim(h(2),[0 100.01]);
yticks(h(1),1.5:0.1:3.3);
yticks(h(2),[0 100]);
grid on; drawnow;

% target phi for given N
h=figure(3);h.Name='target phi versus N';drawnow

myfov = 128;
myangles = 0:0.5:90;
myspokes = round((sin(myangles*pi/360)./(1+sin(myangles*pi/360))).^2*pi*myfov^2);
mytarget = 2*acsc(sqrt(pi./myspokes)*myfov-1)*180/pi; % should == myspokes

plot(myspokes/myfov^2,mytarget); xlabel('N / matrix^2'); ylabel('Angular Rate (deg/dwelltime)');
hold on; plot(myspokes/myfov^2,mytarget/2,'black--'); hold off; grid on; axis tight;
legend({'Target ϕ for a given N','Heuristic (divided by 2)'},'location','NorthWest');
drawnow

clear my*

%% randomize phase versus FWHM
h=figure(4);h.Name='randomize phase versus FWHM';drawnow

myfov = 128;
mynyquist = pi*myfov^2;
myspokes = round(linspace(mynyquist*0.06,mynyquist*0.96,25));

for myj = 1:2
    
    for myk = 1:numel(myspokes)
        
        % corkscew
        expression = ['g.lRadialSpokes = myspokes(myk); g.fov = myfov; g.siz = g.fov; g.maxit = maxit;'];
        expression = [expression 'g.uReadoutType = 2; g.bRandomizePhase = (myj==2)'];

        simulate;
        
        myvol(myj,myk) = vol;
        mystd(myj,myk) = stddev;
        myfwhm(myj,myk) = fwhm3d(1);
        
        % display
        subplot(1,1,1);
        plot(myspokes(1:myk),myfwhm(1:myj,1:myk),'o-');
        xlabel('N'); ylabel('FWHM (Δx)'); title('randomize phase');
        legend({'Constant','Random'}); grid on; ylim([1.5 1.8]); drawnow;
        
    end
    
end

clear my*

%% calculate point density in a shell
h=figure(5);h.Name='density in a shell';drawnow

myfov = 128;
mynyquist = pi*myfov^2;
myspokes = ceil(mynyquist*0.1);

mybandwidth = 80; % lT0 has a dependence on bw... ugly
mylT0 = [1135 1600 2250] * (100/mybandwidth); % matrix=128

for myj = 1:numel(mylT0)
    
    % get the same result with osf=1 but "noisier" so prefer osf=2
    expression = ['g.lRadialSpokes = myspokes; g.fov = myfov; g.siz = g.fov; g.maxit = maxit;'];
    expression = [expression 'g.uReadoutType = 1; g.osf = 2; g.lT0 = mylT0(myj);'];
    expression = [expression '; g.bandwidth = mybandwidth'];
    
    simulate;
    
    % calculate no. points in the shell (width=1)
    for myk = 0.5:0.5:g.siz(1)/2
        myl = (vecnorm(g.om)>=myk-1) & (vecnorm(g.om)<myk);
        mypts(2*myk) = nnz(myl);
        myvol(2*myk) = (4/3)*pi*myk^3 - (4/3)*pi*(myk-1)^3;
    end
    mydensity(:,myj) = mypts./myvol./g.osf; % divide out osf
    myk0(myj) = g.k0;
    
    % fudging to get the first plots to display on top
    if myj==3 && numel(mylT0)==3
        myk=1; while mydensity(myk+1,3)>=mydensity(myk+1,2); mydensity(myk,3)=NaN; myk=myk+1; end;
        myk=1; while mydensity(myk+1,2)>=mydensity(myk+1,1); mydensity(myk,2)=NaN; myk=myk+1; end;
    end
    
    subplot(1,1,1)
    plot(0:0.5:g.siz(1)/2-0.5,mydensity); ylim([0 3.5]); xlim([10 63]); xlabel('k-space radius (cycles/fov)');
    grid on; legend(num2str(myk0','k_0 = %.1f')); ylabel('shell density / osf'); drawnow;

end

clear my*

%% regrid(iterative), regrid(analytical), least squares
h=figure(6);h.Name='reconstruction';drawnow

myfov = 128;
mynyquist = pi*myfov^2;
myspokes = round(linspace(mynyquist*0.06,mynyquist*0.96,25));

mymaxit(1) = 1; % regrid (iterative weighting)
mymaxit(2) = 1; % regrid (analytical weighting)
mymaxit(3) = 10; % least squares (iterative weighting)
mymaxit(4) = 10; % least squares (analytical weighting)

for myj = 1:numel(mymaxit)

    for myk = 1:size(myspokes,2)
        
        % analytical only works for radial (uReadoutType=0)
        expression = ['g.lRadialSpokes = myspokes(myk); g.fov = myfov; g.siz = g.fov;'];
        expression = [expression 'g.uReadoutType = 0; g.maxit = mymaxit(myj);'];
        expression = [expression '; g.dMaxSlew = 100000;']; % minimize ramp (for radial density formula)
        if myj==2 || myj==4; expression = [expression 'g.density = ''radial'';']; end
        
        simulate;
        
        myvol(myj,myk) = vol;
        mystd(myj,myk) = stddev;
        myfwhm(myj,myk) = fwhm3d(1);
        
    end
    
end

subplot(1,1,1);
plot(myspokes./myfov.^2,myfwhm(1,:,1)','o-');
hold on;
plot(myspokes./myfov.^2,myfwhm(2,:,1)','sq-');
plot(myspokes./myfov.^2,myfwhm(3,:,1)','d-');
plot(myspokes./myfov.^2,myfwhm(4,:,1)','--');
hold off;
legend({'regridding (iterative)','regridding (analytic)','least squares (iterative)','least squares (analytic)'}); xlim([0.1 3.1]);
xlabel('N / matrix^2'); ylabel('FWHM (Δx)'); grid on; drawnow;

clear my*

%% 2: biexponential T2 - "sodium" case

NA = [0.4 0.6]; % typical signal fractions
T2 = [ 23   3]; % typical decay constants (ms)

%% voxel fwhm vs. nspokes + T2
h=figure(7);h.Name='voxel fwhm vs. nspokes + T2';drawnow

myfov = 128;
mynyquist = pi*myfov^2;
myspokes = round(linspace(mynyquist*0.06,mynyquist*0.96,25));
        
% 1=radial, 2=density, 3=corkscrew
for myl = 1:3
    
    for myk = 1:numel(myspokes)
        
        expression = ['g.lRadialSpokes = myspokes(myk); g.fov = myfov; g.siz = g.fov;'];
        expression = [expression 'g.uReadoutType = myl-1; g.maxit = maxit;'];
        
        simulate;
        
        myvol(myk,myl) = vol;
        mystd(myk,myl) = stddev;
        mypeak(myk,myl) = peak;
        myfwhm(myk,myl) = fwhm3d(1);
        myreadoutduration(myk,myl) = size(g.om,2) * g.lDwellTime * 1e-6; % ms
        
    end
    
end

subplot(1,1,1);
h = plot(myspokes'./myfov.^2,myfwhm(:,1),'o-');
hold on
h = plot(myspokes'./myfov.^2,myfwhm(:,2)','sq-');
h = plot(myspokes'./myfov.^2,myfwhm(:,3)','d-');
hold off
xlabel('N / matrix^2'); ylabel('FWHM (Δx)');
xlim([0.1 3.1]); ylim([1.5 3.1]); grid on;

% overlay the no T2 case for radial
hold on
plot(save_myspokes_no_T2'./save_myfov_no_T2.^2,save_myfwhm_no_T2(:,1),'--black')
hold off

legend({'radial', 'density adpt','corkscrew','radial (no T2)'});

clear my*

%% voxel fwhm vs. bandwidth + T2
h=figure(8);h.Name='voxel fwhm vs. bandwidth + T2';drawnow

myfov = 128;
myspokes = 10000;
mybandwidth = linspace(5,200,19);

% 1=radial, 2=density, 3=corkscrew
for myl = 1:3
    
    for myj = 1:numel(mybandwidth)
        
        expression = ['g.lRadialSpokes = myspokes; g.fov = myfov; g.siz = g.fov;'];
        expression = [expression 'g.uReadoutType = myl-1; g.maxit = maxit;'];
        expression = [expression 'g.bandwidth = mybandwidth(myj);'];
        
        simulate;
        
        myvol(myj,myl) = vol;
        mystd(myj,myl) = stddev;
        mypeak(myj,myl) = peak;
        myfwhm(myj,myl) = fwhm3d(1);
        
        myreadoutduration(myj,myl) = size(g.om,2) * g.lDwellTime * 1e-6; % ms
        
    end
end

subplot(1,1,1);
h = plot(mybandwidth,myfwhm,'o-');
h(2).Marker = 'sq';
h(3).Marker = 'diamond';
xlabel('Bandwidth (Hz/pixel)'); ylabel('FWHM (Δx)');
xlim([0 200]); ylim([1.5 3.1]); grid on; 
legend({'radial', 'density adpt','corkscrew'}); drawnow;

h=figure(9);h.Name='snr vs. bandwidth + T2';drawnow
snr = (mypeak.*sqrt(myreadoutduration))';
h = plot(mybandwidth,snr,'o-');
h(2).Marker = 'sq';
h(3).Marker = 'diamond';
xlabel('Bandwidth (Hz/pixel)'); ylabel('single voxel SNR');
xlim([0 200]); ylim([0 1]); grid on;
legend({'radial', 'density adpt','corkscrew'}); drawnow;

clear my*
