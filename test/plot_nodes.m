function plot_nodes(nodefile, facefile)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

%plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.');
hold on;
%nn = size(nodes,1);
%plot3(nodes(83:nn,1), nodes(83:nn,2), nodes(83:nn,3), 'r.');
gix=[]+1;
for i=1:length(gix)
    plot3(nodes(gix(i),1), nodes(gix(i),2), nodes(gix(i),3), 'bo');
end
hold off;

grid;

% Plot faces
faces = load(facefile);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
znodes = nodes(:,3);

xfaces = xnodes(faces(:,:));
yfaces = ynodes(faces(:,:));
zfaces = znodes(faces(:,:));

if(size(faces,1) == 1)
    xfaces = xfaces.';
    yfaces = yfaces.';
    zfaces = zfaces.';
end

hold on;
axis off;
for i = 1:size(xfaces,1)
    patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w','EdgeColor','k','LineWidth',2);
end
patch(xfaces(1,:), yfaces(1,:), zfaces(1,:),'g','EdgeColor','k','LineWidth',2);
patch(xfaces(690,:), yfaces(690,:), zfaces(690,:),'b','EdgeColor','k','LineWidth',2);
patch(xfaces(1336,:), yfaces(1336,:), zfaces(1336,:),'r','EdgeColor','k','LineWidth',2);
%patch(xfaces(11816,:), yfaces(11816,:), zfaces(11816,:),'r');
hold off;

view([-45 45]);
camzoom(2.0);

set(gca, 'FontSize', 14);