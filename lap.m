function L = lap(phi,dx,type)
    [M,N] = size(phi);
    L = zeros(M,N);
    l1 = [0,1,0;1,-4,1;0,1,0];
    l2 = [0.25,0.5,0.25;0.5,-3,0.5;0.25,0.5,0.25];
    l3 = [1,1,1;1,-8,1;1,1,1];

    for i = 2:(M-1)
        for j = 2:(N-1)
            sub_mat = phi((i-1):(i+1),(j-1):(j+1));
            if(type==1)
                L(i,j) = sum(sum(sub_mat.*l1));
            end
            if(type==2)
                L(i,j) = sum(sum(sub_mat.*l2));
            end    
            if(type==3)
                L(i,j) = sum(sum(sub_mat.*l3));
            end
        end
    end      
    
    % reflective bc
    L(1,:) = L(3,:);
    L(end,:) = L(M-2,:);
    L(:,1) = L(:,3);
    L(:,end) = L(:,N-2);
    L = L/dx^2;
end   