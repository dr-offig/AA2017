function mov = animate_polar_linear(Z,THETA)

X = abs(real(Z));

[MX,I] = max(X);
[mx,j] = max(MX);

scrsz = get(groot,'ScreenSize');
figure('Position',[1 1 scrsz(4) scrsz(4)])
polarplot(THETA,X(:,j))
axis tight manual
ax = gca;
ax.RLim = [0,mx];
ax.NextPlot = 'replaceChildren';

sz = size(X);
T = sz(2);
mov(T) = struct('cdata',[],'colormap',[]);
for t = 1:T
    polarplot(THETA,X(:,t));
    drawnow
    mov(t) = getframe;
end