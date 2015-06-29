function plot_meshes(nodefile1, facefile1, facelist1, nodefile2, facefile2, facelist2)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes from first mesh
nodes1 = load(nodefile1);
faces1 = load(facefile1);

if (facelist1 == 0)
    facelist1 = [1:size(faces1,1)];
end

nodeset = unique(faces1(facelist1,:));

plot3(nodes1(nodeset,1), nodes1(nodeset,2), nodes1(nodeset,3), 'k.');
hold on;

nodeix = 2403+1;
plot3(nodes1(nodeix,1), nodes1(nodeix,2), nodes1(nodeix,3), 'mo','MarkerSize',12);
nodeix = 2405+1;
plot3(nodes1(nodeix,1), nodes1(nodeix,2), nodes1(nodeix,3), 'ro','MarkerSize',12);
nodeix = 2404+1;
plot3(nodes1(nodeix,1), nodes1(nodeix,2), nodes1(nodeix,3), 'go','MarkerSize',12);

xnodes = nodes1(:,1);
ynodes = nodes1(:,2);
znodes = nodes1(:,3);

xfaces = xnodes(faces1(facelist1,:));
yfaces = ynodes(faces1(facelist1,:));
zfaces = znodes(faces1(facelist1,:));

if(length(facelist1) == 1)
    xfaces = xfaces.';
    yfaces = yfaces.';
    zfaces = zfaces.';
end

for i = 1:length(facelist1)
    patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w');
end

% Plot nodes from second mesh
nodes2 = load(nodefile2);
faces2 = load(facefile2);

nodes2 = nodes2 * 1.000001;

if (facelist2 == 0)
    facelist2 = [1:size(faces2,1)];
end

nodeset = unique(faces2(facelist2,:));

plot3(nodes2(nodeset,1), nodes2(nodeset,2), nodes2(nodeset,3), 'bo');

nodeix = 94965-40335+1; %132478 + 1 - size(nodes1,1);
%plot3(nodes2(nodeix,1), nodes2(nodeix,2), nodes2(nodeix,3), 'go','MarkerSize',12);
nodeix = 54271+1;
%plot3(nodes2(nodeix,1), nodes2(nodeix,2), nodes2(nodeix,3), 'bo','MarkerSize',12);

xnodes = nodes2(:,1);
ynodes = nodes2(:,2);
znodes = nodes2(:,3);

xfaces = xnodes(faces2(facelist2,:));
yfaces = ynodes(faces2(facelist2,:));
zfaces = znodes(faces2(facelist2,:));

if(length(facelist2) == 1)
    xfaces = xfaces.';
    yfaces = yfaces.';
    zfaces = zfaces.';
end

xfaces(:,size(xfaces,2)+1) = xfaces(:,1);
yfaces(:,size(yfaces,2)+1) = yfaces(:,1);
zfaces(:,size(zfaces,2)+1) = zfaces(:,1);

for i = 1:length(facelist2)
    plot3(xfaces(i,:), yfaces(i,:), zfaces(i,:),'b-');
end

xlabel('x');
ylabel('y');
zlabel('z');

hold off;

return;
