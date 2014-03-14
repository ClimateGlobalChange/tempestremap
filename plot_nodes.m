function plot_nodes(nodefile, facefile)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

plot3(nodes(1:82,1), nodes(1:82,2), nodes(1:82,3), 'k.');
hold on;
nn = size(nodes,1);
plot3(nodes(83:nn,1), nodes(83:nn,2), nodes(83:nn,3), 'r.');
hold off;
grid;

% Plot faces
faces = load(facefile);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
znodes = nodes(:,3);

xfaces = xnodes(faces(:,1:4));
yfaces = ynodes(faces(:,1:4));
zfaces = znodes(faces(:,1:4));

if(size(faces,1) == 1)
    xfaces = xfaces.';
    yfaces = yfaces.';
    zfaces = zfaces.';
end

hold on;
for i = 1:size(xfaces,1)
    patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w');
end
hold off;

set(gca, 'FontSize', 14);