%ADI scheme is used to solve the diffussion equation

clear all;
% discretization of teh given domain lamda values are kept similar to the
% previous case while solving using explicit method 
Ax=0;
Bx=2*pi;
By=2*pi;
Ay=0;
N=3;
T=100;
dx=Bx/(N+1);
dy=By/(N+1);
del_x=Ax:dx:Bx;
del_y=Ay:dy:By;
del_t=0.25*(dx)^2;
lamda=del_t/dx^2;
alpha=T/round(del_t,1);

u_bc1=zeros(size(del_y));
u_bc2=zeros(size(del_y));
u_bc3=zeros(size(del_x));
u_sol(N+2,N+2)=zeros(size(N+2:N+2));
u_step_update(N+2,N+2)=zeros(size(N+2:N+2));
time_step=zeros(size(alpha));
temp=0;
error=zeros(N^3);
matrix_A=zeros(N+2,N+2);
%initializing time stepping scheme which update the temperature keeping
% central difference discretization explicit in y and implicit in x 
 u_new(N+2,N+2)=zeros(size(N+2:N+2));
for j=2:alpha+1
    time_step(j)=time_step(j-1)+del_t;
    
    
end

% adding boundary condition to the solution
for i=1:N+2
    u_bc1(i)=(del_x(i)^2*(cos(del_x(i)))); %Bc_horizontal_bottom Ay
    u_bc2(i)=del_x(i)^3  ; %Bc_horizontal top_By
    u_bc3(i)=((2*pi)^2*cos(2*pi))+((del_y(i)/(2*pi))*((2*pi)^3-((2*pi)^2*cos(2*pi))));            %Bc vertical_right_Bx
    
end
u_sol(1,:)=u_bc1;
u_sol(N+2,:)=u_bc2;
u_sol(:,N+2)=u_bc3;
tolerance=0;
iteration=0;
increment=0;
topper=1;
matrix_A=zeros(size(N+1,N+1));
matrix_B=zeros(size(N,N));

%creating tridiagnal matrix which is used to calcuate all the point s while
%transversing through one row at a time and matrix_b while calculating
%values while transversing through y direction 
matrix_A= full(gallery('tridiag',N+1,-lamda/2,1+lamda,-lamda/2));
matrix_b=full(gallery('tridiag',N,-lamda/2,1+lamda,-lamda/2));
matrix_A(1,2)=-lamda;
%while topper> tolerance
for inc=1:2*alpha
    u_new=u_sol;
    increment=increment+del_t;
    % calculating rhs of the equation consisting explicit in y 
    for q=2:N+1   
       for w=1:N+1
                
        f(q,w)=(lamda/2)*(u_new(q-1,w)+u_new(q,w))+(1-lamda)*u_sol(q,w);
        
        
         
        end
            
    end 
    %calculating temp values in x direction using inverse of the
    %tridiagonal matrix
     for p=2:N+1
         for o=1:N+1
             u_new(p,o)=(matrix_A(p,:)')\(f(:,o));
         end 
     end
   %creating rhs while transversing throgh y direction or columns in the
   %matrix case 
    for z=2:N+1   
       for c=1:N+1
           if c== 1 
                
        r(z,c)=(lamda)*(u_new(z,c+1))+(1-lamda)*u_sol(z,c);
           else
        r(z,c)=(lamda/2)*(u_new(z,c-1)+u_new(z,c+1))+(1-lamda)*u_sol(z,c);
           end
           
           
        
         
        end
            
    end 
    r1=r(2:N+1,2:N+1);
    %Calculating for the values in y direction
     for n=2:N+1
         for m=1:N+1
             u_new(n,m)=(matrix_A(n,:)')\(r(:,m));
         end 
     end
     %using error check to iterate untill problem becomes independednt of
     %time 
    error=u_new-u_sol;
    check=max(max(abs(error)));
     u_sol=u_new;
     iteration=iteration+1;
     

    if check==tolerance
      
        
        break;

    end
    u(:,:,inc)=u_new;
        
end
[x,y]=meshgrid(del_x);
figure(1)
contourf(x,y,u_new)

figure(2)
surfc(x,y,u_new)

figure(4)
hold on
plot(u(:,:,10))
plot(u(:,:,20))
plot(u(:,:,30))
plot(u(:,:,40))
plot(u(:,:,50))
plot(u(:,:,60))



