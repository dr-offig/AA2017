function mov = animate_bipolar_log(W,Z,THETA)

X = abs(real(W));
Y = abs(real(Z));


[MX,I] = max(X);
[mx,a] = max(MX);

[MY,J] = max(Y);
[my,b] = max(MY);

[MM,C] = max([mx,my]);

figure
if (C == 1)
  polarplot(THETA,[X(:,a),Y(:,a)]);
else
  polarplot(THETA,[X(:,b),Y(:,b)]);
end
axis tight manual
ax = gca;
ax.RLim = [0,mx];
ax.NextPlot = 'replaceChildren';

sz = size(X);
T = sz(2);
mov(T) = struct('cdata',[],'colormap',[]);
for t = 1:T
    polarplot(THETA,[X(:,t),Y(:,t)]);
    drawnow
    mov(t) = getframe;
end