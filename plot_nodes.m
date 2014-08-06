function plot_nodes(nodefile, facefile)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.');
hold on;
%nn = size(nodes,1);
%plot3(nodes(83:nn,1), nodes(83:nn,2), nodes(83:nn,3), 'r.');
gix=[];
for i=1:length(gix)
    plot3(nodes(gix(i),1), nodes(gix(i),2), nodes(gix(i),3), 'go');
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
for i = 1:size(xfaces,1)
    patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w');
end
%patch(xfaces(11816,:), yfaces(11816,:), zfaces(11816,:),'r');
hold off;

set(gca, 'FontSize', 14);