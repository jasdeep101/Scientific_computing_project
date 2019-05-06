clear all;
%discretization of the problem
Ax=0;
Bx=2*pi;
By=2*pi;
Ay=0;
N=3;
T=10;
dx=Bx/(N+1);
dy=By/(N+1);
del_x=Ax:dx:Bx;
del_y=Ay:dy:By;
del_t=dx^2/4;
lamda=del_t/dx^2;
alpha=T/round(del_t,1);


u_bc1=zeros(size(del_y));
u_bc2=zeros(size(del_y));
u_bc3=zeros(size(del_x));
u_sol(N+2,N+2)=zeros(size(N+2:N+2));
time_step=zeros(size(alpha));


error=zeros(N^3);
 u_new(N+2,N+2)=zeros(size(N+2:N+2));
 u=zeros(size(u_new));
 [x,y]=meshgrid(del_x);

for j=2:alpha+1
    time_step(j)=time_step(j-1)+del_t;
    
end


% assigning BC to the problem
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
 u_new=u_sol;
 %using time stepping scheme using explicit method
for time=1:4*alpha

    increment=increment+del_t;
    
    %condition to calculate for ther fist time step to check if values are
    %updated right 
    if time==1
    for row=2:N+1 
        for col=1:N+1 
            if col==1
                u_new(row,col)=lamda*(2*u_sol(row,col+1)+u_sol(row-1,col)+u_sol(row+1,col))+ (1-4*lamda)*(u_sol(row,col));
            else

                u_new(row,col)=lamda*(u_sol(row,col+1)+u_sol(row-1,col)+u_sol(row+1,col)+u_sol(row,col-1))+ (1-4*lamda)*(u_sol(row,col)) ;
            end
             


        end 
    end
    
    else
        
        
        for row=2:N+1 
        for col=1:N+1 
            if col==1
                u_new(row,col)=lamda*(2*u_sol(row,col+1)+u_sol(row-1,col)+u_sol(row+1,col))+ (1-4*lamda)*(u_sol(row,col));
            else

                u_new(row,col)=lamda*(u_sol(row,col+1)+u_sol(row-1,col)+u_sol(row+1,col)+u_sol(row,col-1))+ (1-4*lamda)*(u_sol(row,col)) ;
            end
           
        end 
        
        end
    end
   u(:,:,time)=u_new;
   %calculating error
error=u_new(row,col)-u_sol(row,col);
check=max(abs(error));


iteration=iteration+1;
if check==tolerance

    break;


end
u_sol=u_new;



    
end




%vizulization

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
plot(u(:,:,max(time)))

    