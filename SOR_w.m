clear all;

Ax=0;
Bx=2*pi;
By=2*pi;
Ay=0;
N=3;

%similar to gauss seidel aprt form adding the SOR-factor to the equation s
%which is calculated using the equation below
dx=Bx/(N+1);
dy=By/(N+1);
del_x=Ax:dx:Bx;
del_y=Ay:dy:By;
temp(N+2,N+2)=zeros(size((N+2),(N+2)));
temp_x(N+2,N+2)=zeros(size((N+2),(N+2)));
temp_y(N+2,N+2)=zeros(size((N+2),(N+2)));

u_loop_sol=zeros(size(N));
f_sol=zeros(size(N));
tolerance=10^-06;
cal_u_vec=zeros(size(N+1));
cal_values=zeros(size(N+1,N+1));
error=zeros(size(N+1));
cal_values_norm=zeros(size(N+1,N+1));


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

Sor_fac=2/1+sin(pi*dx);
K=diag(matrix_A);

U=triu(matrix_A,1)*1;
L=tril(matrix_A)*1/Sor_fac;
L_iv=inv(L);
N_ins_p=inv(L);
D1= full(gallery('tridiag',N+1,0,K(1),0));
D2=D1*(1/Sor_fac)*(1-Sor_fac);
U1=D2-U;



 spectral_radius=max(abs((eig(N_ins_p*U1)))); % value of spectral radius calculated of the triadiagonal matrix formed at each iteration i s0.165 for the sor case 


u_bc1=zeros(size(del_y));
u_bc2=zeros(size(del_y));
u_bc3=zeros(size(del_x));
u_sol(N+2,N+2)=zeros(size(N+2:N+2));
f(N+2,N+2)=zeros(size(N+2:N+2));

u_step_update(N+2,N+2)=zeros(size(N+2:N+2));

for k=1:(N+2)
    for l=N+2:-1:1
        temp(k,l)=(del_x(k));
        temp_x=(temp);
        temp_y=(temp);
        
    end 
end
[x,y]=meshgrid(del_x);



 
for i=1:N+2
    u_bc1(i)=(2*pi-del_y(i))^2*cos(del_y(i)/2); %Bc_vertical_left_Ax
    u_bc2(i)=(del_y(i))*((2*pi)-del_y(i))^2;    %Bc_vertical_right_Bx
    u_bc3(i)=4*pi^2-2*pi*(del_x(i));            %Bc horizontal_bottom_Ay
    for j =1:N+2
        f(i,j)=dy^2*(cos(pi/2*((del_x(j)/pi)+1))*sin(del_y(i)));
    end

end


%Imposing Bc  for the solution at the gridpoints
u_sol(1,:)=(u_bc3(1,:));
u_sol(:,1)=(u_bc1(1,:));
u_sol(:,N+2)=(u_bc2(1,:));

 w=0;
 %Implementing Sor using the sor_fac calculated using research paper found
 %online which say w is aprox equivalent 2/(1+sinpi*del_x)
 
 iteration=0;
 u_new=u_sol;
for e=1:(N+2)^2
    
   
    for t =2 :N+1
        for r = 2:N+2

            if r ==2 
                u_step_update(r,t)=Sor_fac*0.25*(f(r,t)+(u_sol(r,t-1)+u_sol(r,t+1)+u_sol(r-1,t)+u_sol(r+1,t))) +(1-Sor_fac)*(u_sol(r,t));
                u_new(r,t)=Sor_fac*0.25*(f(r,t)+(u_sol(r,t-1)+u_sol(r,t+1)+u_sol(r-1,t)+u_sol(r+1,t)))+(1-Sor_fac)*(u_sol(r,t));
                %u_new=u_sol+u_step_update;
            elseif  r>2 && r<N+2
                u_step_update(r,t)=Sor_fac*0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+u_new(r-1,t)+u_sol(r+1,t)))+(1-Sor_fac)*(u_sol(r,t));
                u_new(r,t)=Sor_fac*0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+u_new(r-1,t)+u_sol(r+1,t)))+(1-Sor_fac)*(u_sol(r,t));
                %u_new=u_sol+u_step_update;
            elseif  r==N+2    
                 u_step_update(r,t)=Sor_fac*0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+2*u_new(r-1,t)))+(1-Sor_fac)*(u_sol(r,t));
                 u_new(r,t)=Sor_fac*0.25*(f(r,t)+(u_new(r,t-1)+u_sol(r,t+1)+2*u_new(r-1,t)))+(1-Sor_fac)*(u_sol(r,t));
                 %u_new(r,t)=u_sol(r,t)+u_step_update(r,t);
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

        
        for o=1:N+1
        error(o)= norm(cal_u_vec(o)-f_sol_trans(o))/(sum(cal_values_norm(o,:)) + norm(f_sol_trans(o)));
        error2(e)=max(max(abs(u_new-u_sol)));
        end
        
    end
     if error2(e)-tolerance<0
        break;


    else
        iteration=iteration+1;
        u_sol=u_new;
    end


    


   
end
%grid_convergence
spectral_radius
 K=8/-log10(spectral_radius)
 Tol_iteration=iteration
 u_new
    
        


figure(1)

contourf(x,y,u_new)

xlabel('x');
ylabel('y');
zlabel('u');
grid;
  