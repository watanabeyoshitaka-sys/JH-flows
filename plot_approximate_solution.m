function plot_approximate_solution(N,h,R,a,uh,eq)

% Plot of the figure at (-1,1)
%hold on
X = (0:N)'*mid(h)-1;
u = [0 uh' 0]';
plot(X,u,'k')
% title('uh')
axis auto
fname = sprintf("uh_R%.3f_a%.3f_eq%.d.pdf", R, a, eq);
exportgraphics(gca, fname, 'ContentType','vector');
hold off

% parameter setting
% gamma = pi/4;       % option
n = ceil(a)*50;       % Number of segments in an arc

% load("app.mat")     % reading approximate solution (option)
np = size(uh,1)+1; xp = -1+2/np:2/np:1-2/np; x = -1+2/n:2/n:1-2/n;
vh = interp1(xp,uh,x,'pchip'); vh=[0,vh,0]; % interpolate values based on the number of segments
theta = linspace(-a, a, 100); % Angle Range of an Arc
r = ones(size(theta)); % An arc with a radius of 1

% Angle Division
angles = linspace(-a, a, n+1);          
positions = [cos(angles); sin(angles)]; 

s = vh;
switch eq % Convert to the function (-1,1)->(-a,a) according to the equation
  case 1
    s = R*s/a;
  case 2
    s = 2*s/a;
  otherwise
    s = s;
end
velocityVectors = [cos(angles); sin(angles)] .* s;

% Normalize by the maximum value
maxNorm = max(vecnorm(velocityVectors));
velocityVectors = velocityVectors / maxNorm;

% Coordinates of an Arc
x_arc = r .* cos(theta);
y_arc = r .* sin(theta);

% coordinates of a ray
x1 = 2*[0, cos(-a)];
y1 = 2*[0, sin(-a)];
x2 = 2*[0, cos(a)];
y2 = 2*[0, sin(a)];

figure;
hold on;
plot(x_arc, y_arc, 'b:', 'LineWidth', 2); 
plot(x1, y1, 'b-', 'LineWidth', 2.);    
plot(x2, y2, 'b-', 'LineWidth', 2.);    
if a<pi/2, pl=-0.085; else pl=0.025; end
text(pl,0,'$0$','Interpreter', 'latex','FontSize',15,'BackgroundColor', '[1 1 1 0.85]','Margin',1);

if a<=pi/2, pl=1.5; else pl=1.75; end
text(pl*cos(-a),pl*sin(-a),'$\theta\!=\!-\alpha$','Interpreter', 'latex','FontSize', 15,'BackgroundColor', '[1 1 1 0.85]','Margin',1);
text(pl*cos(a),pl*sin(a),'$\theta\!=\!\alpha$','Interpreter', 'latex','FontSize', 15,'BackgroundColor', '[1 1 1 0.85]','Margin',1);

%axis equal;
% grid on;
%title(['r = 1, \theta \in [-\alpha, \alpha], a = ', num2str(a)]);
%title(['maximum of vector field: ', num2str(norm(s,inf))]);
%legend('Arc on r=1', 'Line at -a', 'Line at a');

x_start = positions(1,:);
y_start = positions(2,:);

u = velocityVectors(1,:);
v = velocityVectors(2,:);

x_end = x_start + u/4.5;
y_end = y_start + v/4.5;

%quiver(x_start, y_start, u, v, 0.5, 'LineWidth', 1.5);
%scale = pi/(10*a);
quiver(x_start(1:5:n), y_start(1:5:n), u(1:5:n)/5, v(1:5:n)/5, 'off', 'LineWidth', 1.5);

plot(x_end, y_end, 'r-', 'LineWidth', 2);        

%axis equal;
ax = gca;
ax.YTickLabel = cell(size(ax.YTickLabel)); 
ax.XTickLabel = cell(size(ax.XTickLabel)); 
xticks([]);
yticks([]);
ax.XColor = 'none';   
ax.YColor = 'none';   
% axis off;
% grid on;
%title('velocity vectors on unit circle');

allX = [x_arc, x1, x2, positions(1,:), x_start, x_end];
allY = [y_arc, y1, y2, positions(2,:), y_start, y_end];

xmin = min(allX); xmax = max(allX);
ymin = min(allY); ymax = max(allY);
%xlim([-0.15, xmax*1.25]);
xlim([min(xmin,-0.1) xmax*1.25]);
if a<=pi/2, pl=1; else pl=1.2; end
ylim([ymin*pl, ymax*pl]);
daspect([1 1 1]);   

%stval = num2str(norm(s,inf),'%.4f');
%label = sprintf('$\\|u_h\\|_\\infty\\approx%s$',stval);
%text(min(xlim), min(ylim), label,'Interpreter', 'latex','FontSize', 15,'VerticalAlignment','bottom','HorizontalAlignment','left');
%num2str(norm(s,inf))

fname = sprintf("STLine_R%.3f_a%.3f_eq%.d.pdf", R, a, eq);
exportgraphics(gca, fname, 'ContentType','vector');
hold off
