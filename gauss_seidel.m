clear all;

%given domain
Ax=0;
Bx=2*pi;
By=2*pi;
Ay=0;
%Number of discretized points
N=3;

%itilization of del_x and del_y values
dx=Bx/(N+1);
dy=By/(N+1);
del_x=Ax:dx:Bx;
del_y=Ay:dy:By;
u_loop_sol=zeros(size(N));
f_sol=zeros(size(N));
tolerance=10^-06;
cal_u_vec=zeros(size(N+1));
cal_values=zeros(size(N+1,N+1));
error=zeros(size(N+1));
cal_values_norm=zeros(size(N+1,N+1));

%Creating matrix A to check the eroor
matrix_A=zeros(size(N+1,N+1));
for g=1:N+1
    for h =1:N+1
        
        if g==h
            matrix_A(g,h)=-4;
        elseif g-h==1 || h-g==1
             matrix_A(g,h)=1;
             if g==N+1 && h== N
                 matrix_A(g,h)=2;
             end
                 
        end
    end
end
K=diag(matrix_A);

U=-triu(matrix_A,1);
L=tril(matrix_A);
L_iv=inv(L);
N_ins_p=inv(L);


 spectral_radius=max(abs((eig(N_ins_p*U)))); % value of spectral radius calculated of teh triadiagonal matrix formed at each iteration
 


D1=matrix_A-diag(K);
D=inv(diag(K));
B=norm(D1*D); %value of the matrix of the norm which is equivalent to less then 1, 0.572 for gauss seidel.



%Boundary conditions being initiated along with  u_initila_solution which
%consists zeros apart from BC 
u_bc1=zeros(size(del_y));
u_bc2=zeros(size(del_y));
u_bc3=zeros(size(del_x));
u_sol(N+2,N+2)=zeros(size(N+2:N+2));
f(N+2,N+2)=zeros(size(N+2:N+2));
%Extra values storage to create the error check
u_step_update(N+2,N+2)=zeros(size(N+2:N+2));


[x,y]=meshgrid(del_x);



 %Implying Bc to u_initial which for this code and the remaining code will
 %be u_sol also calculating Rhs the forcing function
for i=1:N+2
    u_bc1(i)=(2*pi-del_y(i))^2*cos(del_y(i)/2); %Bc_vertical_left_Ax
    u_bc2(i)=(del_y(i))*((2*pi)-del_y(i))^2;    %Bc_vertical_right_Bx
    u_bc3(i)=4*pi^2-2*pi*(del_x(i));            %Bc horizontal_bottom_Ay
    for j =1:N+2
        f(i,j)=dy^2*(cos(pi/2*((del_x(j)/pi)+1))*sin(del_y(i)));
        %f(i,j)=0;
    end

end


%Imposing Bc  for the solution at the gridpoints
u_sol(1,:)=(u_bc3(1,:));
u_sol(:,1)=(u_bc1(1,:));
u_sol(:,N+2)=(u_bc2(1,:));

 %Implementing Guass seidel
 
 iteration=0;
 u_new=u_sol;
 %u_step_update is used the check the calculated error for every r
 %increment and u_new is the new values created every increment
for e=1:(N+2)^3
    
    for t =2 :N+1
        for r = 2:N+2

            if r ==2 
                u_step_update(r,t)=0.25*(f(r,t)+(u_sol(r,t-1)+u_sol(r,t+1)+u_sol(r-1,t)+u_sol(r+1,t)));
                u_new(r,t)=0.25*(f(r,t)+(u_sol(r,t-1)+u_sol(r,t+1)+u_sol(r-1,t)+u_sol(r+1,t)));
              
            elseif  r>2 && r<N+2
                u_step_update(r,t)=0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+u_new(r-1,t)+u_sol(r+1,t)));
                u_new(r,t)=0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+u_new(r-1,t)+u_sol(r+1,t)));
                
            elseif  r==N+2    
                 u_step_update(r,t)=0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+2*u_new(r-1,t)));
                 u_new(r,t)=0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+2*u_new(r-1,t)));
                 
            end

            u_loop_sol(r)=u_step_update(r,t);
            f_sol(r)=f(r,t)';
            u_loop_sol_trans=u_loop_sol(2:r)';
            f_sol_trans=f_sol(2:r)';
            

        end
        
        cal_u_vec=matrix_A*u_loop_sol_trans;
       
        for z=1:N+1

             for q =1:N+1
                 cal_values(z,q)=matrix_A(z,q)*u_loop_sol_trans(q);
                 cal_values_norm(z,q)=norm(cal_values(z,q));

             end   
         end

     %after calculating u_new it will check for error and tolerace provided
     %which for this code is 10e-08
        for o=1:N+1
        error(o)= norm(cal_u_vec(o)-f_sol_trans(o))/(sum(cal_values_norm(o,:)) + norm(f_sol_trans(o)));
        error2(e)=max(max(abs(u_new-u_sol)));
        end
    
    end
    if error2(e)-tolerance<0
        break;
       

    else
       
        u_sol=u_new;
        
    end
        
    

    


 iteration=iteration+1;  
end
    
 %plots and vizulization.
 
 %grid convergence
 spectral_radius
 K=6/-log10(spectral_radius)
 Tol_iteration=iteration
 u_new
 %here K fetches a value of 12 which is smaller then total number of
 %iteration and iterations are greater then value of K which prooves that
 %grid converges.
figure(1)


contourf(x,y,u_new)

xlabel('x');
ylabel('y');
zlabel('u');
grid;
  