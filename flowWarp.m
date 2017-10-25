function result = flowWarp( flow, img )
% Warps input image using given flow field
%  --> Flows should only be 1-channel

    % construct displacement field from flow
    D = cat(3, flow.x(:,:), flow.y(:,:));
      
    result = zeros(size(img));
    for c = 1:size(img,3)
        
        % do warping on each channel
        result(:,:,c) = imwarp(img(:,:,c), D, 'Interp', 'cubic');
        
    end


%     result = zeros(size(img));
%     for c = 1:size(img,3)
%                 
%         % do warping on each channel
%         [wImg,~] = iat_pixel_warping(img(:,:,c), flow.x, flow.y);
%         result(:,:,c) = wImg;
%         
%     end    


end

