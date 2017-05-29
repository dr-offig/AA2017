% Run the piston model for a set of test frequencies and a particular piston radius and listener radius
a = 0.150;                % 2mm radius driver
r = 1;                    % 1m listening distance
freqs = [250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000];

for f=freqs
  [SPL,SPL_FF,VR,VT,VRL,VTL,THETA,PRESSURE] = piston(f,a,1,1024);
  
%   scrsz = get(groot,'ScreenSize');
%   figure('Position',[1 1 scrsz(3)/2 scrsz(3)/2])
%   polarplot(THETA,SPL./2);
%   title(sprintf('Piston SPL: F=%dHz r=%dm a=%dmm',round(f),round(r),round(a*1000))); 
%   
   save(sprintf('Piston_F%dHz_r%dm_a%dmm.mat',round(f),round(r),round(a*1000)),'SPL','VR','VT','THETA');
%   animate_polar_log(real(VT)./real(VR),10e-9,THETA,sprintf('Piston VT / VR: F=%dHz r=%dm a=%dmm',round(f),round(r),round(a*1000)));

  
end
