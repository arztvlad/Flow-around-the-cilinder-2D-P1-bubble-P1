function [left, right, upper, below, circ] = make_edges(p,id)

%% left %%%
points=[p(id.left,1) p(id.left,2)];
angles = points(:,2);
[~, sorted_indices] = sort(angles);
ID=id.left(sorted_indices);
left=[ID(1:end-1) ID(2:end)];

%% right %%%
points=[p(id.right,1) p(id.right,2)];
angles = points(:,2);
[~, sorted_indices] = sort(angles);
ID=id.right(sorted_indices);
right=[ID(1:end-1) ID(2:end)];

%% upper %%%
points=[p(id.upper,1) p(id.upper,2)];
angles = points(:,1);
[~, sorted_indices] = sort(angles);
ID=id.upper(sorted_indices);
upper=[ID(1:end-1) ID(2:end)];

%% upper %%%
points=[p(id.below,1) p(id.below,2)];
angles = points(:,1);
[~, sorted_indices] = sort(angles);
ID=id.below(sorted_indices);
below=[ID(1:end-1) ID(2:end)];

%% circ %%%
points=[p(id.circ,1) p(id.circ,2)];
center = mean(points, 1);
angles = atan2(points(:,2) - center(2), points(:,1) - center(1));
[~, sorted_indices] = sort(angles);
ID=id.circ(sorted_indices);
circ=[ID(1:end-1) ID(2:end)];
end