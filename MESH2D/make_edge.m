function [left, right, top, bottom,circ,circ_mask] = make_edge(p,ID)

    x_min = min(p(ID(:),1)); x_max = max(p(ID(:),1));
    y_min = min(p(ID(:),2)); y_max = max(p(ID(:),2));
    
    points=[p(ID,1) p(ID,2)];

    % Logické masky pro jednotlivé hrany
    left_mask = (points(:,1) == x_min); 
    right_mask = (points(:,1) == x_max);
    bottom_mask = (points(:,2) == y_min);
    top_mask = (points(:,2) == y_max);



    left=ID(left_mask);
    right=ID(right_mask);
    top=ID(top_mask);
    bottom=ID(bottom_mask);

      bottom=setdiff(bottom,left);
    top=setdiff(top,left);

    %right=setdiff(setdiff(right,bottom),top);

    IDo = [left;right;top;bottom];
    IDi = setdiff(ID,IDo);

    points=[p(IDi,1) p(IDi,2)];
    eps=1e-8;
    circ_mask = and( ((points(:,1)-0.2).^2 +(points(:,2)-0.0).^2) > 0.05^2-eps,...
        ( (points(:,1)-0.2).^2 + (points(:,2)-0.0).^2) < 0.05^2+eps);
    circ=IDi(circ_mask);
  
end