function S = concatFlow( F1, F2 )
% Implementation of 'o' concatenation operator in Suwajanakron paper
% Computes S = F o F' :
%       S_x = F'_x + W_F'(F_x)
%       S_y = F'_y + W_F'(F_y)

    Sx = F2.x + flowWarp(F2, F1.x);
    Sy = F2.y + flowWarp(F2, F1.y);
    
    S = struct('x',Sx,'y',Sy);   

end

