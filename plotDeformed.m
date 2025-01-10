function fig = plotDeformed(file,xn,Tn,u,type,element,scale,sigVM)

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
fig = figure
hold on

% Plot undeformed contour
patch(x0(Tc)',y0(Tc)',z0(Tc)',ones(size(Tc))','facecolor','none','edgecolor',0.5*[1,1,1]);

% Plot deformed structure
patch(x(Tc)',y(Tc)',z(Tc)',ones(size(Tc))','facecolor','none','edgecolor','k');

if ~exist('sigVM','var')
    patch(x(Tn)',y(Tn)',z(Tn)',ones(size(Tn))','facecolor',0.5*[1,1,1],'edgecolor','none');
else
    SigVM = InterpFunction(sigVM,Tn,element);
    if strcmp(element, 'shell')
        patch(x(Tn)',y(Tn)',z(Tn)',SigVM,'facecolor','interp','edgecolor','none');

    elseif strcmp(element, 'beam')
        % Normalize SigVM to the range [0, 1] for colormap
        SigVM_normalized = (SigVM(1,:) - min(SigVM(1,:))) / (max(SigVM(1,:)) - min(SigVM(1,:)));
        
        % Create a colormap (e.g., 'jet' or any other colormap)
        colors = colormap; % Get the RGB values of the colormap
        nColors = size(colors, 1);
        
        % Map normalized SigVM values to colormap indices
        colorIndices = round(SigVM_normalized * (nColors - 1)) + 1;
        pointColors = colors(colorIndices, :); % RGB colors for each point

        % Plot line segments between specified pairs
        for i = 1:size(Tn, 1)
            % Get the indices of the pair
            p1 = Tn(i, 1);
            p2 = Tn(i, 2);
            
            % Extract the coordinates for the pair
            x_pair = [x(p1), x(p2)];
            y_pair = [y(p1), y(p2)];
            z_pair = [z(p1), z(p2)];
            
            % Interpolate the color for the line based on the pair's SigVM values
            lineColor = mean([pointColors(i, :); pointColors(i, :)], 1);
            
            % Plot the line segment
            plot3(x_pair, y_pair, z_pair, '-', 'Color', lineColor, 'LineWidth', 1.5);
        end

    end

    % Color axis
    clim([0,max([SigVM(:);1])]);
    cb = colorbar;
    set(cb,'ticks',[0,max([SigVM(:);1])]);
    yl = ylabel(cb,'Von Mises Stress [MPa]');
    yl.Position(1) = min(xlim(cb));
    yl.VerticalAlignment = 'bottom';
end

% Additional options
title(sprintf('Deformation for %s elements. Scale = %g',type,scale));
set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none');
axis equal
axis vis3d
axis equal;  % Ensures equal scaling of x, y, and z axes
view(40,20);



end

function sign = InterpFunction(sige,Tn,element)
    a = [-1,1,1,-1];
    b = [-1,-1,1,1];
    sign = zeros(size(Tn'));
    for e = 1:size(Tn,1)
        
        if strcmp(element, 'shell')
            N4 = zeros(4,4);
            for k = 1:4
                for i = 1:4
                    N4(k,i) = (1+a(i)*a(k)/sqrt(3))*(1+b(i)*b(k)/sqrt(3))/4;
                end
            end
            sign(:,e) = N4\(sige(e,:)');

        elseif strcmp(element, 'beam')
            sign(:,e) = sige(e,:)';
        end
    end
end

