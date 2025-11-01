% Domain length
lx=1; 
ly=0.5;

% no. of elements
n_x=input('nx'); 
n_y=input('ny');

% This specifies the the stretch direction 
% 1 for positive, -1 for negative and 0 for no stretch
% e=input("x direction");
% w=input('y direction');

x_size= lx/n_x;
y_size= ly/n_y;

% For equidistant mesh
dx= zeros(1,n_x);
dy= zeros(1,n_y);

% Corner points 
z_x=zeros(1,n_x);
z_y=zeros(1,n_y);

z_x2=ones(1,n_x);
z_y2=ones(1,n_y);

% Defining grid array
grid_x=zeros(1,n_x+1);
grid_y=zeros(1,n_y+1);

%Defining nodes array
node_x=zeros(n_x,n_y);
node_y=zeros(n_x,n_y);

% Node value 
Temp = zeros(n_x+2,n_y+2);
T_new = zeros(n_x+2,n_y+2);

% Coefficient matrix 
A_e=zeros(n_x,n_y);
A_w=zeros(n_x,n_y);
A_n=zeros(n_x,n_y);
A_s=zeros(n_x,n_y);

a_e=zeros(n_x+2,n_y+2);
a_w=zeros(n_x+2,n_y+2);
a_n=zeros(n_x+2,n_y+2);
a_s=zeros(n_x+2,n_y+2);

Mesh = input('0 for no refinement; 1 for Refinement- ');

% MESHING

% Equidistant Mesh
if Mesh == 0

    for i=1:n_x
   
    grid_x(i+1)=grid_x(i)+x_size;
 
    end
    

    for i=1:n_y
  
    grid_y(i+1)=grid_y(i)+y_size;
 
    end
 
end  

% Stretched mesh 
if Mesh == 1

    %Stretch Factor for refinement
    st_x=0.9;
    st_y=0.9;
    
    %Length of first element using GP
    dx(1)=((lx)*(1-st_x)/(1-(st_x^(n_x))));

    dy(1)=(ly*(1-st_y)/(1-(st_y^n_y)));
    
    for i=1:n_x
    
    dx(i+1)=dx(i)*st_x;
    
    end

    for i=1:n_y

    dy(i+1)=dy(i)*st_y;

    end

    % X refinement in 1 direction

    for i=1:n_x

        grid_x(i+1)=grid_x(i)+dx(i);

    end

    % Negative refinement
    %      for i=1:n_x
    % 
    %         grid_x(i+1)=grid_x(i)+dx((n_x+1)-i);
    %      end
    

    
    % Refinement in x in both sides
    
    % dx(1)=((lx/2)*(1-st_x)/(1-(st_x^(n_x/2))));
    % for i=1:n_x/2
    % 
    % dx(i+1)=dx(i)*st_x;
    % 
    % end

    % for i=1:n_x/2
    % 
    % grid_x(i+1)=grid_x(i)+dx(((n_x/2)+1)-i);
    % 
    % end
    % 
    % 
    % for j= 1:n_x/2
    % 
    % p = n_x/2;
    % grid_x(p+1+j)=grid_x(p+j)+dx(j);
    % 
    % end

   

    % Y direction ( NO refinement)

    % for i=1:n_y
    % 
    % grid_y(i+1)=grid_y(i)+y_size;
    % 
    % end
     
    
    % Positive refinement

    for i=1:n_y

        grid_y(i+1)=grid_y(i)+dy(i);

    end
    
    % Negative refinement
    %      for i=1:n_y
    % 
    %         grid_y(i+1)=grid_y(i)+dy((n_y+1)-i);
    %      end
    


end 

cell_x=zeros(1,n_x);
cell_y=zeros(1,n_y);

for l=1:n_x
cell_x(l)= grid_x(l+1)-grid_x(l);
end

for l=1:n_y
cell_y(l)= grid_y(l+1)-grid_y(l);
end
   

%Loop for node location
for j=1:n_x
    for k=1:n_y
       node_x(j,k)=(grid_x(j)+grid_x(j+1))/2;
       node_y(j,k)=(grid_y(k)+grid_y(k+1))/2;

    end
end



% Coefficient update

for j=1:n_x
    for k=1:n_y
     % East
     if j == n_x
         A_e(j,k)= (32*((node_y(n_x,k)/0.5) + 1)*cell_y(k))/(grid_x(n_x+1)-grid_x(n_x));

    else
        A_e(j,k)= (16*((node_y(j,k)/0.5) + 1)*cell_y(k))/(grid_x(j+1)-grid_x(j));

    end 
    % West
    if j == 1
       A_w(j,k)= (32*((node_y(j,k)/0.5) + 1)*cell_y(k))/(grid_x(2));    
     
    else
       A_w(j,k)= (16*((node_y(j,k)/0.5) + 1)*cell_y(k))/(grid_x(j)-grid_x(j-1));

    end 

   % North
    if k == n_y

        A_n(j,k)=(32*((node_y(j,n_y)/0.5) + 1)*cell_x(j))/(grid_y(n_y + 1)-grid_y(n_y));
         
    else
        A_n(j,k)=(16*((node_y(j,k)/0.5) + 1)*cell_x(j))/(grid_y(k+1)-grid_y(k));
       
    end 

    % South 
    if k == 1

         A_s(j,k)=(32*((node_y(j,1)/0.5) + 1)*cell_x(j))/(grid_y(2));
   
    else
         A_s(j,k)=(16*((node_y(j,k)/0.5) + 1)*cell_x(j))/(grid_y(k)-grid_y(k-1));
        
    end 



    end
end


% Boundary condition 
Temp(:,1)=15;
Temp(:,n_y+2)=10;
% Temp(1,:)=12;

N_x=node_x(1,:);
N_y=node_y(1,:);

for k = 2:n_y+1
Temp(n_x+2,k)= 5*(1-(N_y(k-1)/ly)) + (15*sin((pi*N_y(k-1))/ly));
end

a_w(2:n_x+1,2:n_y+1)=A_w;
a_e(2:n_x+1,2:n_y+1)=A_e;
a_n(2:n_x+1,2:n_y+1)=A_n;
a_s(2:n_x+1,2:n_y+1)=A_s;
a_p=zeros(n_x+2,n_y+2);

Dx=zeros(1,n_x+2);
Dy=zeros(1,n_y+2);
Dx(2:n_x+1)=cell_x;
Dy(2:n_y+1)=cell_y;


% disp(a_w);
for m=2:n_x+1
  for n=2:n_y+1
    a_p(m,n) = a_w(m,n) + a_e(m,n) + a_s(m,n) + a_n(m,n)+(30000*Dx(m)*Dy(n));
  end
end
% disp(a_p);

residual = 1e5;
iter=0;
residual_list = []; 

while residual > 1e-6


% Guass scidel
    for m = 2:n_x+1
        for n = 2:n_y+1
    
        Temp(m,n)= (a_e(m,n)*Temp(m+1,n)+ a_w(m,n)*Temp(m-1,n) + a_n(m,n)*Temp(m,n+1) + a_s(m,n)*Temp(m,n-1)+(500000*Dx(m)*Dy(n)))/a_p(m,n);
       
        end
    end 
    Temp(1,:) = Temp(2,:);
    % Temp(1,:)= Temp(2,:)+(-16*((node_y(a-1,b-1)/0.5) + 1)* Dx(2)*5000);
    residual = 0;
    iter=iter+1;


% Residual Calculation
    for a = 2:n_x+1
        for b = 2:n_y+1

        RHS = a_e(a,b)*Temp(a+1,b)+ a_w(a,b)*Temp(a-1,b) + a_n(a,b)*Temp(a,b+1) + a_s(a,b)*Temp(a,b-1)+(500000*Dx(a)*Dy(b));
        LHS =Temp(a,b)*a_p(a,b);

        r= RHS-LHS;
        residual= residual + r;

        end
    end 
residual_list(end+1) = residual;
% disp(residual);
end

disp(iter);


% Heat flux
Qx=zeros(n_x,n_y);
Qy=zeros(n_x,n_y);

 for a = 2:n_x+1
    for b = 2:n_y+1

        Qx(a-1,b-1)=-16*((node_y(a-1,b-1)/0.5) + 1)* (Temp(a+1,b) - Temp(a-1,b)) / (2*Dx(a));
        Qy(a-1,b-1)=-16*((node_y(a-1,b-1)/0.5) + 1)* (Temp(a,b+1) - Temp(a,b-1)) / (2*Dy(b));

    end
 end

% Temperature Contour
[X, Y] = meshgrid(grid_x,grid_y); 
figure;
contourf(X, Y, transpose(Temp(1:end-1,1:end-1)),20,'LineStyle','none');     
colorbar;
xlabel('L'); ylabel('H');
title('Temperature Contours');

% Heat Flux Vector
[Xp, Yp] = meshgrid(node_x(:,1),node_y(1,:));
figure;
quiver(Xp, Yp, transpose(Qx),transpose(Qy), 'AutoScale','on','AutoScaleFactor', 2);  
xlabel('x'); ylabel('y');
title('Heat Flux Vectors');
axis equal;




% Mesh Plotting
 
% [X,Y]=meshgrid(grid_x,grid_y);
% figure
% plot(X,grid_y); hold on
% plot(transpose(grid_x),transpose(Y));
% 
% % Plotting node points
% plot(node_x, node_y, 'k.');
% 
% % Plotting corner points
% plot(node_x,z_y,'k.');
% plot(z_x,node_y,'k.');
% plot(node_x,z_y2,'k.');
% plot(z_x2,node_y,'k.');
% xlim([0,1]);
% ylim([0,0.5]);

% Residual vs Iteration plotting

figure;
semilogy(1:iter, residual_list, 'b-*', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Residual');
title('Residual vs No. of iterations');
grid("on");
