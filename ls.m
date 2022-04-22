close all;
clear vars-global

%% Level set 
I=imread('circle.jpg');%imread('overlap_sq.png');%imread('overlap_circ.jpg');%imread('terminal.gif');
%I = I(:,:,1);
I = rgb2gray(I);
I = I/max(max(I));
I = double(I);
I = imresize(I,[256,256]);
I = I - 0.5; % for circle.jpg
[M,N] = size(I);
%I=im2double(I);
circle_xi = I;
circle_xi(abs(circle_xi - 0.5) < 0.00001 & circle_xi > 0) = 0; % for the circle image
circle_xi(circle_xi ~= 0) = 1/2;
counter = 0; old_counter = 0;
for r = 1:M
    cind = find(circle_xi(r,:)~=0);
    if ~isempty(cind)
        old_counter = counter;
        counter = counter + 1;
    end
    if counter == 1 || counter == old_counter
        1;
    else
        circle_xi(r,(min(cind)+1):(max(cind)-1)) = 1;
    end    
end

counter = 0; old_counter = 0;
for c = 1:N
    rind = find(circle_xi(:,c)~=0);
    if ~isempty(rind)
        old_counter = counter;
        counter = counter + 1;
    end
    if counter == 1 || counter == old_counter
        1;
    else
        circle_xi((min(rind)),c) = 1/2;
        circle_xi((max(rind)),c) = 1/2;
    end    
end
cbdry = zeros(M,N);
for r = 1:M
    c = find(circle_xi(r,:)>0);
    cbdry(r,min(c)) = 1/2;
    cbdry(r,max(c)) = 1/2;
end
for c = 1:N
    r = find(circle_xi(:,c)>0);
    cbdry(min(r),c) = 1/2;
    cbdry(max(r),c) = 1/2;
end
%%{

square = zeros(256);
square(44:223,43) = 1/2;
square(43,44:223) = 1/2;
square(223,44:223) = 1/2;
square(44:223,223) = 1/2;
square(43,43) = 1/2;

phi0 = bwdist(cbdry)+ bwdist(1-cbdry)-70; % circle boundary
%phi0 = bwdist(square)-60; % square boundary
[M,N] = size(phi0);

% Values at t = 0
phi = phi0; its = 100; mu = 100; dt = 0.5; name = 'circle';

figure();contour(phi, [0 0], 'r');
title(['0 /' num2str(its) ' Iterations']);
hold off;
drawnow;

idx = 0; % plot index
xi_news = {1:4};

for n=1:its
    % Curvature
    phi_dub = double(phi);
    phi_pad = padarray(phi_dub,[1,1],1,'both'); % getting the 'ghost' points

    % central difference
    fy = (phi_pad(3:end,2:N+1)-phi_pad(1:M,2:N+1));
    fx = (phi_pad(2:M+1,3:end)-phi_pad(2:M+1,1:N));
    fyy = phi_pad(3:end,2:N+1)+phi_pad(1:M,2:N+1)-2*phi_dub;
    fxx = phi_pad(2:M+1,3:end)+phi_pad(2:M+1,1:N)-2*phi_dub;
    fxy = (1/4).*(phi_pad(3:end,3:end)-phi_pad(1:M,3:end)+phi_pad(3:end,1:N)-phi_pad(1:M,1:N));
    K = ((fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5))).*(fx.^2+fy.^2).^(0.5);

    K(1,:) = eps;
    K(end,:) = eps;
    K(:,1) = eps;
    K(:,end) = eps;
    K = K./max(max(abs(K)));

    % Evolution
    F = mu*K; % 
    F = F./max(max(abs(F)));

    phi = phi+dt*F; 

    contour(phi, [0 0], 'r');
    title([num2str(n) '/' num2str(its) ' Iterations']);
    hold off;
    drawnow;
    
    
    if(mod(n,its/4)==0)
        idx = idx + 1;
        phi_news{idx} = phi;
    end
    close;
    
end
%-- End

nImages = length(phi_news);
figure(4);
for idx = 1:nImages
    h{idx} = subplot(2,2,idx);
    hLine{idx} = contour(phi_news{idx}, [0 0], 'r');
    title([num2str(its/4*idx) '/' num2str(its) ' Iterations']);
end
imwrite(getframe(gcf).cdata, name+"/dt"+dt+"_ls.jpg", "Quality", 100)

%}
