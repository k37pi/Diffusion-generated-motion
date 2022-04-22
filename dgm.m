close all;
clear vars-global; 

%% Inputs
trip_pt = imread('triple_point.PNG');
trip_pt = rgb2gray(trip_pt);
trip_pt = imresize(trip_pt,[400,400]);
trip_pt = double(trip_pt);
trip_pt(trip_pt~=0) = 1/2;
trip_pt = trip_pt - 1/2;
trip_pt = trip_pt*(-1);

square1 = zeros(size(trip_pt));
square2 = zeros(size(trip_pt));
square3 = zeros(size(trip_pt));
square1(1:227,:) = trip_pt(1:227,:);
square2(:,1:206) = trip_pt(:,1:206); square2(224,206)=0;square2(225:233,206)=0;square2(228:231,204)=1/2;square2(224,205)=0;
square3(:,207:end) = trip_pt(:,207:end); square3(227:end,:) = trip_pt(227:end,:); square3(226,206)=1/2; square3(227,204)=0; square3(225,206)=1/2;

[r1,c1] = find(square1 > 0);
for c = 1:length(c1)
    rows = 1:r1(c);
    square1(rows,c1(c)) = 1/2;
end       
%figure();imagesc(square1);

[r2,c2] = find(square2 > 0);
for r = 1:length(r2)
    cols = 1:c2(r);
    square2(r2(r),cols) = 1/2;
end       
%figure();imagesc(square2);

[r3,c3] = find(square3 > 0);
for r = 1:length(r3)
    cols = c3(r):size(trip_pt,2);
    square3(r3(r),cols) = 1/2;
end       
%figure();imagesc(square3);

%{
trip_pt(trip_pt==0) = -1;
trip_pt(trip_pt==1/2) = 1;

square1(square1==0) = -1;
square1(square1==1/2) = 1;
square2(square2==0) = -1;
square2(square2==1/2) = 1;
square3(square3==0) = -1;
square3(square3==1/2) = 1;
%}


trip_pt_xi0 = {square1,square2,square3};
figure(1)
for k1 = 1:4
    h{k1} = subplot(2,2,k1);
    if k1 == 1
        hLine{k1} = imagesc(trip_pt);
    else
        hLine{k1} = imagesc(trip_pt_xi0{k1-1});
    end    
end
clf(1);

%%{
self_int = imread('self_intersection.jpg');
self_int = rgb2gray(self_int);
self_int = imresize(self_int,[256,256]);
self_int = double(self_int);
self_int = self_int/max(max(self_int));
self_int(abs(self_int-0) < 0.07) = 0;
self_int(self_int~=0) = 1;
self_int = self_int - 1;
self_int = self_int*(-1);
[M,N] = size(self_int);

imagesc(self_int);

square1 = zeros(size(self_int));
square1(:,1:150) = self_int(:,1:150); 
square1(182:195,142:150) = 0;
square1(196:end,100:150) = 0;
[r1,c1] = find(square1 > 0);

square2 = square1;
square2(:,131:end) = 0;
square2(1:70,:) = 0; 
square2(170:end,:) = 0; 
square2(70:75,124:130) = 0; 
[r2,c2] = find(square2 > 0);
inds = min(unique(r2)):max(unique(r2));

temp_square2 = square2;
square1(temp_square2~=0) = 0;

for r = inds
    rinds = find(r2==r);
    cinds = c2(rinds);
    square2(r,min(cinds):max(cinds)) = 1;
end

for r = 1:M
    cind = find(square1(r,:)~=0,1,'last');
    square1(r,1:cind) = 1;
end    
square1(square2~=0) = 0;

square3 = zeros(M,N);
square3(1:83,:) = self_int(1:83,:);
square3(70:90,35:110) = 0;
for r = 1:83
    cind = find(square3(r,:)~=0);
    square3(r,min(cind):max(cind)) = 1;
end    

square4 = zeros(M,N);
square4(1:100,151:end) = self_int(1:100,151:end);
for r = 1:100
    cind = find(square4(r,:)~=0);
    if length(cind)<10
         square4(r,min(cind):end) = 1;
    else
        square4(r,min(cind):max(cind)) = 1;
    end    
end    

square5 = zeros(M,N);
square5(84:end,131:end) = self_int(84:end,131:end);
square5(175:185,1:145) = 0;
square5(179:end,1:148) = 0;
%square5(84,184) = 0;
for c = 131:N
    rind = find(square5(:,c)~=0);
    square5(min(rind):max(rind),c) = 1;
end

square67 = zeros(M,N);
square67(179:end,:) = self_int(179:end,:);

square6 = square67;
square6(:,1:120) = 0;
square6(175:190,1:145) = 0;
for c = 120:N
    rind = find(square6(:,c)~=0);
    square6(min(rind):end,c) = 1;
end

square7 = square67;
square7(:,149:end) = 0;
for r = 175:M
    cind = find(square7(r,:)~=0);
    square7(r,min(cind):max(cind)) = 1;
end

self_int_xi0 = {square1,square2,square3,square4,square5,square6,square7};
figure(2)
for k1 = 1:8
    h{k1} = subplot(2,4,k1);
    if k1 == 1
        hLine{k1} = imagesc(self_int);
    else
        hLine{k1} = imagesc(self_int_xi0{k1-1});
    end    
end
clf(2);

dif = imread('dif.PNG');
dif = rgb2gray(dif);
dif = imresize(dif,[256,256]);
dif = double(dif);
dif = dif/max(max(dif));
dif(dif<1) = 0;
kd = kappa(dif);

dx = 10;
dif_dif = 4*del2(dif,dx);
for d = 1:10
    dif_dif = 4*del2(dif_dif,dx);
    dif_dif = dif_dif/max(max(dif_dif));
end
%imagesc(dif_dif);

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
circle_xi(circle_xi ~= 0) = 1;
bdry = zeros(M,N);
for r = 1:M
    c = find(circle_xi(r,:)>0);
    bdry(r,min(c)) = 1/2;
    bdry(r,max(c)) = 1/2;
end
for c = 1:N
    r = find(circle_xi(:,c)>0);
    bdry(min(r),c) = 1/2;
    bdry(max(r),c) = 1/2;
end

circle_xi0 = {circle_xi};
circle_xi0{2} = 1-circle_xi0{1};

square = zeros(256);
square(44:223,43) = 1;
square(43,44:223) = 1;
square(223,44:223) = 1;
square(44:223,223) = 1;
square(43,43) = 1;

square_xi = square;
[M,N] = size(square_xi);
for r = 1:256
    c = find(square(r,:)>0);
    square_xi(r,min(c):max(c)) = 1;
end
square_xi0 = {square_xi};
square_xi0{2} = 1-square_xi0{1};

cross = zeros(256,256);
cross(128,:) = 1;
cross(:,128) = 1;
[M,N] = size(cross);

cross1 = zeros(M,N);
cross1(1:128,1:128) = 1;
cross2 = zeros(M,N);
cross2(1:128,129:256) = 1;
cross3 = zeros(M,N);
cross3(129:256,129:256) = 1;
cross4 = zeros(M,N);
cross4(129:256,1:128) = 1;
cross_xi0 = {cross1,cross2,cross3,cross4};


%-- End
%% DGM
name = "trip_pt";
I = trip_pt;
[M,N] = size(I);
xi0 = trip_pt_xi0;
l = length(xi0);
its = 16;
dx = 0.5; dt = 0.001; D = 0.01; stencil = 2;

Dx = zeros(M,N); Dy = zeros(M,N); cx = round(M/2); cy = round(N/2);
for i = 1:M
    for j = 1:N
        r = sqrt((i)^2+(j)^2); %%r = sqrt((cx-i)^2+(cy-j)^2);
        %rterm = (-sin(r))/r; Dx(i,j) = rterm*(i-cx); Dy(i,j) = rterm*(j-cy);
        rterm = (-sin(r))/r; Dx(i,j) = rterm*(i); Dy(i,j) = rterm*(j);
        %rterm = (sec(r))^2/r; Dx(i,j) = rterm*(i); Dy(i,j) = rterm*(j);
        %rterm = (2*r)/r; Dx(i,j) = rterm*(i); Dy(i,j) = rterm*(j);
        %Dx(i,j) = 1/i; Dy(i,j) = 0; 
        %Dx(i,j) = log(i+j); Dy(i,j) = log(i+j); 
    end
end    

dtype = "d1";

type = "nc"; % "c" for constant, any other string for variable
xi = xi0; %xi = cellfun(@(x) imgaussfilt(x),xi,'UniformOutput',false);

figure(1);
for k1 = 1:(l+1)
    h{k1} = subplot(2,round((l+1)/2),k1);
    if k1 == 1
        hLine{k1} = imagesc(I);
    else
        hLine{k1} = imagesc(xi0{k1-1});
    end    
end

K = xi0; tmin = zeros(length(xi),1); tmax = zeros(length(xi),2); %dt = zeros(length(xi),1); 
%{
for i = 1:length(xi)
    K{i} = kappa(xi{i});
    K_ind = find(K{i}>0.001);
    tmin(i) = dx/(D*( max(abs(K{i}(K_ind)))));
    tmax(i) = 1/(D*(min(min(abs(K{i}(K_ind))))));
    %dt(i) = (tmax(i)+tmin(i))/2;
end
%}
%dt = max(tmax);

%xi = cellfun(@(x) x/max(max(abs(x))), xi, 'UniformOutput',false);

change = xi0; xi_Lap = xi0; 
idx = 0; % plot index
xi_news = {1:4}; im = {1;4};
%xi0norm = zeros(its,1);

for n=1:its
    % Discrete Laplacian
    %xi_Lap = del2(xi,dx);
    
        if lower(type) == "c" % Constant D 
         for d = 1:25   
           for len = 1:l
                %xi_Lap = cellfun(@(x) del2(x,dx),xi,'UniformOutput',false);
                %change = cellfun(@(x) dt*D*x,xi_Lap,'UniformOutput',false);
                %xi_Lap{len} = 4*del2(xi{len},dx);
           
             xi_Lap{len} = lap(xi{len},dx,stencil);
             change{len} = dt*D*xi_Lap{len};
             xi{len} = xi{len}+change{len};
             xi{len} = xi{len}/max(max(abs(xi{len})));
           end 
         end 
       else  
        for d = 1:100   
            for len = 1:l
                xi{len} = grad(xi{len},Dx,Dy,dx,dt);
                xi{len} = xi{len}/max(max(abs(xi{len})));
            end
        end
       end
        %xi = cellfun(@plus, xi, change, 'UniformOutput',false);%cellfun(@(x) x+change,xi,'UniformOutput',false); 
        %xi = cellfun(@(x) x/max(max(abs(x))), xi, 'UniformOutput',false);
        %K = cellfun(@(x) kappa(x),xi,'UniformOutput',false);
 

    figure(2);
    for k1 = 1:(l)
        h{k1} = subplot(2,round((l)/2),k1);
        hLine{k1} = imagesc(xi{k1});
    end
    sgtitle([num2str(n) '/' num2str(its) ' Iterations']);
    hold off;
    drawnow;

    %%{
    
    thr = 10^-5;
    %{
    for reg = 1:(length(xi)-1)
        for reg1 = (reg+1):length(xi)
            xi_comp = xi{reg} ~= 0 & xi{reg1} ~= 0;
            indcomp = find(xi_comp~=0);
            for inds = indcomp
                if (xi{reg}(inds)) >= (xi{reg1}(inds))
                    xi{reg1}(inds) = 0;
                else
                    xi{reg}(inds) = 0; 
                end    
            end
        end    
    end
    %}

    for reg = 1:(length(xi)-1)
        for reg1 = (reg+1):length(xi)
            xi_comp = abs(xi{reg}) >= thr & abs(xi{reg1}) >= thr;
            indcomp = find(xi_comp~=0);
            for inds = indcomp
                if (xi{reg}(inds)) >= (xi{reg1}(inds))
                    xi{reg1}(inds) = 0;
                else
                    xi{reg}(inds) = 0; 
                end    
            end
        end    
    end


    figure(3);
    for k1 = 1:(l)
        h{k1} = subplot(2,round((l)/2),k1);
        hLine{k1} = imagesc(xi{k1});
    end
    sgtitle([num2str(n) '/' num2str(its) ' Iterations']);
    hold off;
    drawnow;

    xi_new = zeros(M,N);
    for l = 1:length(xi)
        %xi{l}(xi{l} ~= 0) = 1;
        xi{l}(abs(xi{l}) > thr) = 1;
        %xi{l}(xi{l} > 1/2) = 1;
        xi{l}(xi{l} ~= 1) = 0;
        
        % row search 
        for r = 1:M
            cind = find(xi{l}(r,:)>0);
            diff_cind = diff(cind);
            diff_cind_ind = find(diff_cind > 1);

            if isempty(diff_cind_ind) % no gaps
                
                if length(cind)==1
                   xi_new(r,cind) = 1;
                end         
    
                if length(cind)>1
                    if (max(cind) == N && min(cind) ~= 1) ||  (max(cind) ~= N && min(cind) == 1)
                        xi_new(r,min(cind)) = 1;
                    end    
                    if (max(cind) ~= N && min(cind) ~= 1) 
                        xi_new(r,min(cind)) = 1;
                        xi_new(r,max(cind)) = 1;
                    end    
                end                
            end   
        end

        % col search 
        for c = 1:N
            rind = find(xi{l}(:,c)>0);
            diff_rind = diff(rind);
            diff_rind_ind = find(diff_rind > 1);

            if isempty(diff_rind_ind) % no gaps
                
                if length(rind)==1
                   xi_new(rind,c) = 1;
                end         
    
                if length(rind)>1
                    if (max(rind) == M && min(rind) ~= 1) ||  (max(rind) ~= M && min(rind) == 1)
                        xi_new(min(rind),c) = 1;
                    end    
                    if (max(rind) ~= M && min(rind) ~= 1) 
                        xi_new(min(rind),c) = 1;
                        xi_new(max(rind),c) = 1;
                    end    
                end                
            end   
        end

    end
    
    %xi = cellfun(@(x) imgaussfilt(x),xi,'UniformOutput',false);

    %imagesc(xi);
    figure(1);
    for k1 = 1:(l+1)
        h{k1} = subplot(2,round((l+1)/2),k1);
        if k1 == 1
            hLine{k1} = imagesc(xi_new);
        else
            hLine{k1} = imagesc(xi{k1-1});
        end    
    end
    sgtitle([num2str(n) '/' num2str(its) ' Iterations']);
    hold off;
    drawnow;
    %}
    
    fig = figure;
    if(mod(n,its/4)==0)
        idx = idx + 1;
        imagesc(xi_new)
        title([num2str(n) '/' num2str(its) ' Iterations'])
        drawnow
        frame = getframe(fig);
        im{idx} = frame2im(frame);
        xi_news{idx} = xi_new;
    end
    close;
    
    %if mod(n,25)==0
    %    imwrite(xi_new, name+"/dx"+dx+"_dt"+dt+"_D"+dtype+"its_"+n+".jpg", "Quality", 100) 
    %end
     
    %{
    for i = 1:length(xi)
        %K{i} = kappa(xi{i});
        %tmin(i) = dx/(D*(min(min(abs(K{i}+eps)))));
        %tmax(i) = 1/(D*(max(max(abs(K{i}.^2+eps)))));
        %dt(i) = tmax(i);

        K{i} = kappa(xi{i});
        K_ind = find(K{i}>0.001);
        if ~isempty(K{i}(K_ind))
            tmin(i) = dx/(D*( max(abs(K{i}(K_ind)))));
            tmax(i) = 1/(D*(min(min(abs(K{i}(K_ind))))));
            %dt(i) = (tmax(i)+tmin(i))/2;
            %dt(i) = tmin(i)/2;
        end    
    end
    %}
end

nImages = length(xi_news);
figure(4);
for idx = 1:nImages
    h{idx} = subplot(2,2,idx);
    hLine{idx} = imagesc(xi_news{idx});
    title([num2str(its/4*idx) '/' num2str(its) ' Iterations']);
end
imwrite(getframe(gcf).cdata, name+"/dx"+dx+"_dt"+dt+"_D"+dtype+"_stencil"+stencil+".jpg", "Quality", 100)

filename = name+"/dx"+dx+"_dt"+dt+"_D"+dtype+"_stencil"+stencil+".gif"; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

   

%figure();plot(xi0norm);      

%-- End

function K = kappa(phi)
    phi_dub = double(phi);
    [M,N] = size(phi_dub);
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
end   


%figure(5);
%subplot(2,2,1); imagesc(bdry);
%subplot(2,2,2); imagesc(square);
%subplot(2,2,3); imagesc(trip_pt);
%subplot(2,2,4); imagesc(self_int);

