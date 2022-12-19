%% Plot boxes as patches
%
% DOC_boxes is 1x7 or 7x1 vector for the colors of the boxes.

function [cHandle] = plot_boxes(DOC_boxes, varargin)

if ~any(strcmp(varargin, 'nonew'))
    %figure('color', 'white', 'position', [715,599,432,261])
    figure('color', 'white', 'position', [5,5,432,261])
end

if any(strcmp(varargin, 'clabel'))
    idx = find(strcmp(varargin, 'clabel'));
    cl  = varargin{idx+1}; 
else
    cl = 'DOC [mmolC/m^3]';
end

if any(strcmp(varargin, 'caxlim'))
    idx = find(strcmp(varargin, 'caxlim'));
    clim  = varargin{idx+1}; 
else
    clim = [];
end

PE = setup_model;
PS = setup_run(PE);
PO = setup_ocean(PE,PS);

% smallest box is square, respectively set width of all boxes
W = repmat(sqrt(min(PO.A)), 1, 7); % width of all boxes [m]

% calculate length of all boxes from their width
L = PO.A./W; % length of all boxes [km]
H = PO.H; % height of all boxes [m]
assert(all(W.*L - PO.A==0), 'wrong width or length')

% NADW
V_NADW = PO.V(6);
h = (H(3)+H(5))-H(2);
hx = (V_NADW-L(2)*W(2)*h-L(4)*W(4)*h)/((L(2)+L(3)+L(4))*W(2));
NADW_bottom = 1000+hx; % NADW bottom depth

% AABW
V_AABW = PO.V(7);
hy = (V_AABW-L(1)*W(1)*(750+hx))/(L(1)*W(1)+(L(2)+L(3)+L(4))*W(2));
AABW_bottom = 1000+hx+hy;

% Set up coordinates of boxes 
AAx = [0 0 L(1) L(1) L(1) L(1) L(1) L(1)];
AAy = [H(1) 0 0 H(1) H(1) H(1) H(1) H(1)];
SAx = [L(1) L(1) L(1)+L(2) L(1)+L(2) L(1)+L(2) L(1)+L(2) L(1)+L(2)...
    L(1)+L(2)];
SAy = [H(2) 0 0 H(2) H(2) H(2) H(2) H(2)];
LLx = [L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3)) sum(L(1:3)) ...
    sum(L(1:3)) sum(L(1:3)) sum(L(1:3))];
LLy = [H(3) 0 0 H(3) H(3) H(3) H(3) H(3)];
NAx = [sum(L(1:3)) sum(L(1:3)) sum(L(1:4)) sum(L(1:4)) sum(L(1:4)) ...
    sum(L(1:4)) sum(L(1:4)) sum(L(1:4))];
NAy = [H(4) 0 0 H(4) H(4) H(4) H(4) H(4)];
TCx = [L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3)) sum(L(1:3)) ...
    sum(L(1:3)) sum(L(1:3)) sum(L(1:3))];
TCy = [H(3)+H(5) H(3) H(3) H(3)+H(5) H(3)+H(5) H(3)+H(5) H(3)+H(5)...
    H(3)+H(5)];
NADWx = [L(1) L(1) L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3))...
    sum(L(1:4)) sum(L(1:4))];
NADWy = [NADW_bottom H(2) H(2) H(3)+H(5) H(3)+H(5) H(4) H(4) NADW_bottom];
AABWx = [0 0 L(1) L(1) sum(L(1:4)) sum(L(1:4)) sum(L(1:4)) sum(L(1:4))];
AABWy = [AABW_bottom H(1) H(1) NADW_bottom NADW_bottom AABW_bottom...
     AABW_bottom AABW_bottom];
cx = [AAx; SAx; LLx; NAx; TCx; NADWx; AABWx]';
cy = [AAy; SAy; LLy; NAy; TCy; NADWy; AABWy]';

patch(cx, cy, DOC_boxes)

ylabel('Depth [m]'), axis tight
cHandle = colorbar; ylabel(cHandle, cl)
%colormap(flipud(cbrewer('div', 'RdYlBu', 100, 'PCHIP')))
%colormap(parula)
%%%%% change here for colormap
colm=[141,36,0;...
    218,89,5;...
    254,152,40;...
    255,219,121;...
    255,249,191]/255;
colmi=flipud(interp1(1:1:length(colm),colm,1:0.00001:length(colm)));
colormap(colmi)
%%%%%%end colormap
legend('off')
set(gca, 'YDir', 'reverse')

%bluemap = [flipud(cbrewer('seq', 'Blues', 30)); 1 1 1];
%redmap = [1 1 1; cbrewer('seq', 'Reds', 30)];
%diffmap = [bluemap; redmap];
%diffmap=parula;
diffmap=colmi;
% for comparison plot
if any(strcmp(varargin, 'diff'))
    if all(DOC_boxes>0)
%         diffmap = flipud(hot(30));
        colormap(redmap)
%         ylabel(cHandle, sprintf('Difference Model-Data \n[mmolC/m³]'))
        caxis([0 max(10, max(DOC_boxes))])
    elseif all(DOC_boxes<=0)      
        colormap(bluemap)
%         ylabel(cHandle, sprintf('Difference Model-Data \n[mmolC/m³]'))
        caxis([min(-10, min(DOC_boxes)) 0])
    else
%         load('diffmap'); %cbrewer('div', 'PuOr', 30); % diffmap = colormap(gca)
        colormap(diffmap)
%         ylabel(cHandle, sprintf('Difference Model-Data \n[mmolC/m³]'))
        caxis([min(-10, -max(abs(DOC_boxes))) max(10, max(abs(DOC_boxes)))])
    end  
end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

ax = gca;
ax.XTick = ax.XLim;
ax.XTickLabel = {'South', 'North'};
ax.YTick = [0, 1000, 2000, 3000];
caxis([clim(1) clim(2)])

%set(cHandle, 'YTickLabel', ...
%    cellstr(num2str(reshape(get(cHandle, 'YTick'),[],1),'%0.0f')) )
end
