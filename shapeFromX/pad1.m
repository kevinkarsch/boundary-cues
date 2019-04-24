function Zp = pad1(Z)

  Zp = [Z(1,1), Z(1,:), Z(1,end); 
        Z(:,1), Z,      Z(:,end);
        Z(end,1), Z(end,:), Z(end,end)];

% Zp = [[Z(1,1); Z(:,1); Z(end,1)], [Z(1,:); Z; Z(end,:)], [Z(1,end); Z(:,end); Z(end,end)]];


