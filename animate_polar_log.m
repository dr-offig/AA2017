function mov = animate_polar_log(Z,q,THETA,description)

X = max(log10(abs(real(Z))./q)/log10(1/q),zeros(size(Z)));

[MX,I] = max(X);
[mx,j] = max(MX);

scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(3)/2])
polarplot(THETA,X(:,j))
axis tight manual
ax = gca;
ax.RLim = [0,mx];
ax.NextPlot = 'replaceChildren';
title(description);

sz = size(X);
T = sz(2);
mov(T) = struct('cdata',[],'colormap',[]);
for t = 1:T
    polarplot(THETA,X(:,t));
    drawnow
    mov(t) = getframe;
end