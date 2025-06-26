n = 66;
b = bar3(tran_hist_deep(:,locs(n)-15:locs(n)+15)./max(tran_hist_deep(:,locs(n)-15:locs(n)+15),[],2));

% Set white background
set(gcf,'Color','white');  % Figure background to white
set(gca,'Color','white');  % Axes background to white

% Set white grid lines and ticks
set(gca,'GridColor','white');  % Grid lines to white
set(gca,'XColor','white');     % X-axis ticks to white
set(gca,'YColor','white');     % Y-axis ticks to white
set(gca,'ZColor','white');     % Z-axis ticks to white
set(gca, 'LineWidth', 200);
axis off;

%% (Commented out code - translated but not executed)
% % Set fonts
% set(gca, 'Fontname','Arial','FontSize',32,'FontWeight','bold','Color','white');
% 
% ylim([-3.5 3.5]);
% 
% % Set white axis labels
% hx = xlabel('Time(s)', 'FontSize', 32, 'FontName', 'Arial','Rotation',43,'Color','white');
% hy = ylabel(['Anterior-Posterior ' ...
%     'bins'], 'FontSize', 32, 'FontName', 'Arial','Rotation', -45,'Color','white');
% set(hx, 'Position', [5, -5, -2]);  % X-axis label position
% set(hy, 'Position', [-5, 5, -2]);  % Y-axis label position
% 
% % Set axis line width and ticks
% set(gca, 'LineWidth', 3);
% 
% set(gca, 'XTick', linspace(1, 17, 3)); 
% 
% colormap sky;
% daspect([2,1,0.15]);