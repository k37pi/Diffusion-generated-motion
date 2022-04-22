function G = grad(phi,Dx,Dy,dx,dt) % Dx = pd of D wrt x
    [M,N] = size(phi);
    G = phi;

    for i = 2:(M-1)
        for j = 2:(N-1)
           alpha_x = dt*Dx(i,j)/(2*dx); alpha_y = dt*Dy(i,j)/(2*dx);
           G(i,j) =  phi(i,j) + alpha_x*(phi(i+1,j)- phi(i-1,j)) + alpha_y*(phi(i,j+1)- phi(i,j-1));
        end
    end      
    
    % reflective bc
    G(1,:) = G(3,:);
    G(end,:) = G(M-2,:);
    G(:,1) = G(:,3);
    G(:,end) = G(:,N-2);
end   