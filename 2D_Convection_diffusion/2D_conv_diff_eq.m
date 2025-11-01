close all
clear
clc

% read xc
load xc.dat
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
% For a 10x10 mesh, ni=no of nodes=12,nim1=no of grid lines=11
ni=nim1+1;

% read yc
load yc.dat
njm1=length(yc);
nj=njm1+1;

% read u
load u.dat
u2d=reshape(u,ni,nj);
% read v
load v.dat
v2d=reshape(v,ni,nj);

% compute the x-coordinates of the cell centres
for i=2:nim1
   xp(i)=0.5*(xc(i)+xc(i-1));
end
xp(1)=xc(1);
xp(ni)=xc(nim1);
%
% take the transpose of x
xp=xp';


% compute the y-coordinates of the cell centres
for j=2:njm1
   yp(j)=0.5*(yc(j)+yc(j-1));
end
yp(1)=yc(1);
yp(nj)=yc(njm1);
%
% take the transpose of y
yp=yp';


%To compute the delta_xe:dist b/w node 
for i=1:nim1
d_x(i)= xp(i+1)-xp(i);
end

for i=1:njm1
d_y(i)= yp(i+1)-yp(i);
end

 
for i=1:ni-2
    Dx(i)= xc(i+1)-xc(i);
end

for i=1:nj-2
    Dy(i)= yc(i+1)-yc(i);
end



Ae=zeros(ni,nj);
Aw=zeros(ni,nj);
An=zeros(ni,nj);
As=zeros(ni,nj);

De=zeros(ni,nj);
Dw=zeros(ni,nj);
Dn=zeros(ni,nj);
Ds=zeros(ni,nj);

Fe=zeros(ni,nj);
Fw=zeros(ni,nj);
Fn=zeros(ni,nj);
Fs=zeros(ni,nj);

rho = 1;
Cp = 1/50;

% Dx  = zeros(1,ni-2);
% Dy  = zeros(1,nj-2);

for i= 2:nim1
    for j = 2:njm1
Fxe(i)= 0.5 * Dx(i-1) / d_x(i);
Fxw(i) = 0.5 * Dx(i-1) / d_x(i-1);

% Different indexing for velocity as extracting from data file

uw(i,j) = Fxw(i) * u2d(j,i-1) + (1-Fxw(i)) * u2d(j,i);
ue(i,j) = Fxe(i) * u2d(j,i+1) + (1-Fxe(i)) * u2d(j,i);

Fyn(j)= 0.5 * Dy(j-1) / d_y(j);
Fys(j)= 0.5 * Dy(j-1) / d_y(j-1);

vn(i,j) = Fyn(j) * v2d(j,i+1) + (1-Fyn(j)) * v2d(j,i);
vs(i,j) = Fys(j) * v2d(j,i-1) + (1-Fys(j)) * v2d(j,i);


    end
end


for i =2:nim1
 for j=2:njm1


Fe(i,j)=rho*ue(i,j)*Dy(j-1);
Fw(i,j)=rho*uw(i,j)*Dy(j-1);
Fn(i,j)=rho*vn(i,j)*Dx(i-1);
Fs(i,j)=rho*vs(i,j)*Dx(i-1);


De(i,j)= ((1/50)*d_y(j-1)) / d_x(i);
Dw(i,j)= ((1/50)*d_y(j-1)) / d_x(i-1);
Dn(i,j)=((1/50)*d_x(i-1)/ d_y(j));
Ds(i,j)=((1/50)*d_x(i-1)/ d_y(j-1));

Aw(i,j)=max([Fw(i,j),Dw(i,j)+(Fw(i,j)/2),0]);
Ae(i,j)=max([-Fe(i,j),De(i,j)-(Fw(i,j)/2),0]);
An(i,j)=max([-Fn(i,j),Dn(i,j)-(Fn(i,j)/2),0]);
As(i,j)=max([Fs(i,j),Ds(i,j)+(Fs(i,j)/2),0]);

 end
end

Ap = Aw + Ae +As + An ;


Temp=25*ones(ni,nj);

Temp(ni,7:end)=50;
Temp(1,22:27)=20;
% Temp(2:nim1,1)= 15 ;

residual = 1e5;
iter=0;
residual_list = []; 

f=rho*1*0.136*30;

method = input(" Enter 1 For Gauss seidel , 2 for TDMA -");


if method == 1
    % Guass scidel
    
    while residual>1e-2
        for m = 2:nim1
            for n = 2:njm1
    
            Temp(m,n)= (Ae(m,n)*Temp(m+1,n)+ Aw(m,n)*Temp(m-1,n) + An(m,n)*Temp(m,n+1) + As(m,n)*Temp(m,n-1))/Ap(m,n);
    
            end
        end 
    
          % Left
          Temp(1,2:21) = Temp(2,2:21);
    
          % % Bottom
          Temp(2:nim1,1)=Temp(2:nim1,2);
    
          % Top
          Temp(2:nim1,nj)=Temp(2:nim1,nj-1);
    
          % Right
          Temp(ni,2:6)= Temp(ni-1,2:6); 
    
          residual = 0;
          iter=iter+1;
    
    
    % Residual Calculation
        for a = 2:nim1
            for b = 2:njm1
    
            RHS = Ae(a,b)*Temp(a+1,b)+ Aw(a,b)*Temp(a-1,b) + An(a,b)*Temp(a,b+1) + As(a,b)*Temp(a,b-1);
            LHS =Temp(a,b)*Ap(a,b);
    
            r= abs(RHS-LHS)/f;
            residual= residual + r;
    
            end
        end 
    residual_list(end+1) = residual;
    % disp(residual);
    end

end


if method == 2

        P=zeros(ni,1);
    Q=zeros(ni,1);
    d=zeros(ni,1);
    
    % for o = 1:15
    while residual>1e-3
    for j = 2:njm1
    
        Temp(ni,7:end)=50;
        Temp(1,22:27)=20;
    
    
        for h = 2:nim1
            d(h) = An(h,j)*Temp(h,j+1) + As(h,j)*Temp(h,j-1);
        end
    
        if (j >= 22 && j <= 27)       
        Tw = 20;
        d(2) = An(2,j)*Temp(2,j+1) + As(2,j)*Temp(2,j-1) + Aw(2,j)*Tw;
        P(2) = Ae(2,j) / Ap(2,j);
        Q(2) = d(2) / Ap(2,j);
        
        else
    
        d(2) = An(2,j)*Temp(2,j+1) + As(2,j)*Temp(2,j-1);
        P(2) = Ae(2,j) / (Ap(2,j) - Aw(2,j)); 
        Q(2) = d(2) / (Ap(2,j) - Aw(2,j));
        end
    
        % P(2)= Ae(2,j)/(Ap(2,j));
        % Q(2)= d(2) / (Ap(2,j));
    
        for k = 3:nim1
    
        P(k) = Ae(k,j) / (Ap(k,j)-Aw(k,j)*P(k-1));
        Q(k) = (Aw(k,j)*Q(k-1) + d(k) ) / ( Ap(k,j) -Aw(k,j)*P(k-1));
        end
    
        for n = nim1:-1:2
    
            Temp(n,j)= P(n)*Temp(n+1,j) + Q(n);
        end 
    
    
    end
    
          % Left
          Temp(1,2:21) = Temp(2,2:21);
    
          % Bottom
          Temp(2:nim1,1)=Temp(2:nim1,2);
    
          % Top
          Temp(2:nim1,nj)=Temp(2:nim1,nj-1);
    
          % Right
          Temp(ni,2:6)= Temp(ni-1,2:6); 
    
        residual = 0;
        iter=iter+1;
    
    
    % Residual Calculation
        for a = 2:nim1
            for b = 2:njm1
    
            RHS = Ae(a,b)*Temp(a+1,b)+ Aw(a,b)*Temp(a-1,b) + An(a,b)*Temp(a,b+1) + As(a,b)*Temp(a,b-1);
            LHS =Temp(a,b)*Ap(a,b);
    
            r = (abs(RHS-LHS))/f;
            residual= residual + r;
    
            end
        end 
    
    residual_list(end+1) = residual;
    
    % disp(residual);
    
    end

end 

% Heat Flux 

% Top and Bottom 

Qb = sum ((1/50)*(Temp(2:end-1,1)-Temp(2:end-1,2))/d_y(1));
Qt = sum ((1/50)*(Temp(:,nj)-Temp(:,nj-1))/d_y(end));

UA = 1;
UC = 1;
Ql=0;
Qr=0;

for x = 22:nj-1

 Qlcond= (1/50)* (Temp(2,x) - Temp(1,x)) / d_x(1);
Qlconv = (0.136 * rho *d_y(x-1)* UA * Temp(1,x));
Ql = Ql + Qlconv + Qlcond;
end

for x = 2:6
    
Qrcond= (1/50) * (Temp(ni,x) - Temp(nim1,x)) / d_x(end);
Qrconv = (0.136 * d_y(x-1) *rho * UC * Temp(ni,x));
Qr =Qr+ Qrconv + Qrcond;
end



disp("Bottom heat flux = " + num2str(Qb));
disp("Top heat flux = " + num2str(Qt));
disp("Right heat flux = " + num2str(Qr));
disp("Left heat flux = " + num2str(Ql));


[X, Y] = meshgrid(xc,yc); 
figure;
contourf(X, Y, transpose(Temp(1:end-1,1:end-1)),25,'LineStyle','none');     
colorbar;
xlabel('L'); ylabel('H');
title('Temperature Contours');




figure;
semilogy(1:iter, residual_list, 'b-*', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Residual');
title('Residual vs No. of iterations');
grid("on");

% vec= 5;
% figure
% quiver(xp,yp,u2d,v2d,vec)
% axis('equal');
% xlabel('x'); ylabel('y')
% %print vectxy.ps -deps
