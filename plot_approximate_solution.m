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

% パラメータ設定
% gamma = pi/4;         % 傾き角 [今後のオプション]
n = ceil(a)*50;       % 円弧の分割数

% load("app.mat")       % 近似解の読込み
np = size(uh,1)+1; xp = -1+2/np:2/np:1-2/np; x = -1+2/n:2/n:1-2/n;
vh = interp1(xp,uh,x,'pchip'); vh=[0,vh,0]; % 分割数に応じて値を補間

theta = linspace(-a, a, 100); % 円弧の角度範囲
r = ones(size(theta)); % 半径1の円弧

% 角度の分割
angles = linspace(-a, a, n+1);          % [-a,a]の内点を分割
positions = [cos(angles); sin(angles)]; % 単位円上の座標 (2×n 行列)

% 速度ベクトルの計算
%xvals = angles / a;                    % 正規化された x 値
%dw = R * f(xvals);      % f(x) に R を掛けた値
%velocityVectors = [dw; tan(gamma) * dw]; % 速度ベクトル (2×n 行列)
%s = flip(vh); % 関数を原点中心に反転 
s = vh;
switch eq % 方程式に応じて関数 (-1,1)->(-a,a) に変換
  case 1
    s = R*s/a;
  case 2
    s = 2*s/a;
  otherwise
    s = s;
end
velocityVectors = [cos(angles); sin(angles)] .* s;

% 無限ノルムで正規化（最大値で割る）
maxNorm = max(vecnorm(velocityVectors));
velocityVectors = velocityVectors / maxNorm;

% 円弧の座標
x_arc = r .* cos(theta);
y_arc = r .* sin(theta);

% 放射線の座標（原点から角度-aとaへ）
x1 = 2*[0, cos(-a)];
y1 = 2*[0, sin(-a)];
x2 = 2*[0, cos(a)];
y2 = 2*[0, sin(a)];

% プロット
figure;
hold on;
plot(x_arc, y_arc, 'b:', 'LineWidth', 2); % 円弧
plot(x1, y1, 'b-', 'LineWidth', 2.);    % 線1（-a方向）
plot(x2, y2, 'b-', 'LineWidth', 2.);    % 線2（a方向）
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

% 始点座標
x_start = positions(1,:);
y_start = positions(2,:);

% ベクトル成分
u = velocityVectors(1,:);
v = velocityVectors(2,:);

% 放射状ベクトルの先端座標
x_end = x_start + u/4.5;
y_end = y_start + v/4.5;

% 矢印を描画
%quiver(x_start, y_start, u, v, 0.5, 'LineWidth', 1.5);
%scale = pi/(10*a);
quiver(x_start(1:5:n), y_start(1:5:n), u(1:5:n)/5, v(1:5:n)/5, 'off', 'LineWidth', 1.5);

% 先端を結ぶ線を追加
plot(x_end, y_end, 'r-', 'LineWidth', 2);        % 先端を結んだ曲線

%axis equal;
ax = gca;
ax.YTickLabel = cell(size(ax.YTickLabel)); % Y軸の目盛り上の数値を削除
ax.XTickLabel = cell(size(ax.XTickLabel)); % X軸の目盛り上の数値を削除
xticks([]);
yticks([]);
ax.XColor = 'none';   % x軸の線を消す
ax.YColor = 'none';   % y軸の線を消す
% axis off;
% grid on;
%title('velocity vectors on unit circle');

% すべてのX/Yデータを結合
allX = [x_arc, x1, x2, positions(1,:), x_start, x_end];
allY = [y_arc, y1, y2, positions(2,:), y_start, y_end];

xmin = min(allX); xmax = max(allX);
ymin = min(allY); ymax = max(allY);
%xlim([-0.15, xmax*1.25]);
xlim([min(xmin,-0.1) xmax*1.25]);
if a<=pi/2, pl=1; else pl=1.2; end
ylim([ymin*pl, ymax*pl]);
daspect([1 1 1]);   % XとYの比率だけを固定

%stval = num2str(norm(s,inf),'%.4f');
%label = sprintf('$\\|u_h\\|_\\infty\\approx%s$',stval);
%text(min(xlim), min(ylim), label,'Interpreter', 'latex','FontSize', 15,'VerticalAlignment','bottom','HorizontalAlignment','left');
%num2str(norm(s,inf))

% ファイル名の生成と保存
fname = sprintf("STLine_R%.3f_a%.3f_eq%.d.pdf", R, a, eq);
exportgraphics(gca, fname, 'ContentType','vector');
hold off
