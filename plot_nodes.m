function plot_nodes(nodefile, facefile)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.');

grid;

% Plot faces
faces = load(facefile);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
znodes = nodes(:,3);

xfaces = xnodes(faces(:,1:4));
yfaces = ynodes(faces(:,1:4));
zfaces = znodes(faces(:,1:4));

for i = 1:size(xfaces,1)
    patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w');
end

hold on;
plot3(6.6531356996e-01, -5.1429550210e-01, -5.1429550210e-01, 'r.');
hold off;

% bnd = [529 25;528 529;527 528;526 527;29 526;30 29;31 30;32 31;523 32;522 523;521 522;35 521;36 35;360 36;361 360;370 361;371 370;380 371;381 380;390 381;399 390;408 399;417 408;416 417;425 416;424 425;433 424;432 433;108 432;107 108;512 107;513 512;514 513;104 514;103 104;102 103;516 102;517 516;518 517;519 518;520 519;97 520;96 97;276 96;275 276;274 275;265 274;264 265;255 264;254 255;245 254;236 245;227 236;218 227;219 218;210 219;211 210;202 211;203 202;204 203;24 204;25 24] + 1;
% 
% hold on;
% for i = 1:size(bnd,1)
%    %arrow3d([xnodes(bnd(i,1)) ynodes(bnd(i,1)) znodes(bnd(i,1))], [xnodes(bnd(i,2)) ynodes(bnd(i,2)) znodes(bnd(i,2))]);
%    plot3([xnodes(bnd(i,1)) xnodes(bnd(i,2))], [ynodes(bnd(i,1)) ynodes(bnd(i,2))], [znodes(bnd(i,1)) znodes(bnd(i,2))], 'b-', 'LineWidth', 4);
% end
% hold off;

set(gca, 'FontSize', 14);