function plotDeformed(file,xn,Tn,u,scale,sigVM)

% Load data
switch file
    case 'shell'
        load('shell.mat','Tc');
    case 'wing'
        load('wing.mat','Tc');
end

% Precompute
x0 = xn(:,1);
y0 = xn(:,2);
z0 = xn(:,3);
x = x0+scale*u(1:6:end);
y = y0+scale*u(2:6:end);
z = z0+scale*u(3:6:end);

% Initalize figure
figure
hold on

% Plot undeformed contour
patch(x0(Tc)',y0(Tc)',z0(Tc)',ones(size(Tc))','facecolor','none','edgecolor',0.5*[1,1,1]);

% Plot deformed structure
patch(x(Tc)',y(Tc)',z(Tc)',ones(size(Tc))','facecolor','none','edgecolor','k');

if ~exist('sigVM','var')
    patch(x(Tn)',y(Tn)',z(Tn)',ones(size(Tn))','facecolor',0.5*[1,1,1],'edgecolor','none');
else
    SigVM = InterpFunction(sigVM,Tn);
    patch(x(Tn)',y(Tn)',z(Tn)',SigVM,'facecolor','interp','edgecolor','none');
    % Color axis
    caxis([0,max([SigVM(:);1])]);
    cb = colorbar;
    set(cb,'ticks',[0,max([SigVM(:);1])]);
end

% Additional options
title(sprintf('Deformation scale = %g',scale));
set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none');
axis equal
axis vis3d
view(40,20);

end

function sign = InterpFunction(sige,Tn)
    a = [-1,1,1,-1];
    b = [-1,-1,1,1];
    sign = zeros(size(Tn'));
    for e = 1:size(Tn,1)
        N4 = zeros(4,4);
        for k = 1:4
            for i = 1:4
                N4(k,i) = (1+a(i)*a(k)/sqrt(3))*(1+b(i)*b(k)/sqrt(3))/4;
            end
        end
        sign(:,e) = N4\(sige(e,:)');
    end
end