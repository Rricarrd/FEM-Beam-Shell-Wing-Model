function plotModes(file,Phi,freq,imodes)

% Load data
switch file
    case 'shell'
        load('shell.mat','xn','Tn','Tc');
    case 'wing'
        load('wing.mat','xn','Tn_wb','Tn_rb','Tn_sk','Tc');
        Tn = [Tn_wb;Tn_rb;Tn_sk];
end

% Precompute
scale = 1;
x0 = xn(:,1);
y0 = xn(:,2);
z0 = xn(:,3);

figure
for i = 1:length(imodes)
    I = imodes(i);
    x = x0+scale*Phi(1:6:end,I)/Phi(find(abs(Phi(:,I))==max(abs(Phi(:,I))),1),I);
    y = y0+scale*Phi(2:6:end,I)/Phi(find(abs(Phi(:,I))==max(abs(Phi(:,I))),1),I);
    z = z0+scale*Phi(3:6:end,I)/Phi(find(abs(Phi(:,I))==max(abs(Phi(:,I))),1),I);
    subplot(3,3,i)
    hold on
    patch(x0(Tc)',y0(Tc)',z0(Tc)',ones(size(Tc))','facecolor','none','edgecolor',0.5*[1,1,1]);
    patch(x(Tc)',y(Tc)',z(Tc)',ones(size(Tc))','facecolor','none','edgecolor','k');
    patch(x(Tn)',y(Tn)',z(Tn)',ones(size(Tn))','EdgeColor','none','FaceColor',0.5*[1,1,1]);
    view(40,20);
    set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none');
    title(sprintf('f_{%i} = %.3f Hz',I,freq(I)));
    axis equal;
    axis vis3d;
    h = gca;    % Get the handle to the current axes
    set(h, 'FontSize', 20); % Set font size
end

end