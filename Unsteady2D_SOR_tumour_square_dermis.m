format long
%dermis layer properties
Rod=1200;
cd=3.3;
kd=0.445;
%fat layer properties
Rof=1000;
cf=2.674;
kf=0.185;
%muscle layer properties
Rom=1085;
cm=3.8;
km=0.51;
wb=0.0027;
%blood properties(for perfusion)
Rob=1060;
cb=3.77;
%tumour properties
kt=0.56;
Rot=1050;
ct=3.77;

%properties
Tinf=25;
Tb=37;
h=10;
L=0.060; %Length in vertical direction
H=0.0116; %Length in horizontal direction
M=300; %No. of nodes in the vertical direction
N=72; %No of nodes in the horizontal direction
dx=H/(M-1);
dy=L/(N-1);
q=500; %heat flux
dt=1;
K=wb*Rob*cb*Tb;
%position of the boundaries between layers
p=floor((1.6/11.6)*M);
g=floor((3.6/11.6)*M);
%tumour in muscle
a=4; %towards i
b=11; %towards j
p1x=round((1.5/11.6)*M)+a;
p1y=round((1.5/60)*N)+b;

%constants to make the discretized equation easier to represent
a1x=(kd*dt)/(Rod*cd*dx^2);
a1y=(kd*dt)/(Rod*cd*dy^2);
a2x=(kf*dt)/(Rof*cf*dx^2);
a2y=(kf*dt)/(Rof*cf*dy^2);
a3x=(km*dt)/(Rom*cm*dx^2);
a3y=(km*dt)/(Rom*cm*dy^2);
a4x=((kd+kf/2)*dt)/(((Rod+Rof)/2)*((cd+cf)/2)*dx^2);
a4y=((kd+kf/2)*dt)/(((Rod+Rof)/2)*((cd+cf)/2)*dy^2);
a5x=((kf+km/2)*dt)/(((Rof+Rom)/2)*((cf+cm)/2)*dx^2);
a5y=((kf+km/2)*dt)/(((Rof+Rom)/2)*((cf+cm)/2)*dy^2);
atx=(kt*dt)/(Rot*ct*dx^2);
aty=(kt*dt)/(Rot*ct*dy^2);
%initialization
T_SOR=37*ones(M,N);
T=T_SOR;
T_new=T;
T_old = T;
T_previous_nt = T_old ;
error=1000;
tol=10e-10;
alpha=1.2; %SOR coefficient
time=10000;
timeT_Surface=zeros(time,1); %matrix to store temperature at each timestep, at the skin surface
timeT_muscle=zeros(time,1); %matrix to store temperature at each timestep at the muscle surface

k=1;
while(k<time+1)
     sor_iter = 1 ;
     %the while loop below performs SOR for each timestep, until solution
     %converges at that timestep.
     while(error > tol) 
         %heat flux boundary condition
         for j=2:floor(N/3)
             T(1,j)=((T(2,j)*(kd*dy/dx))+((T(1,j+1)+T(1,j-1))*(kd*dx/(2*dy)))+(q*dy))/((kd*dy/dx)+(kd*dx/dy));
             %T(1,j)=T(2,j)+(q*dx/kd);-for personal refernce
         end
         %convective boundary condition
         for j=floor(N/3)+1:N-1 
             T(1,j)=((T(2,j)*(kd*dy/dx))+((T(1,j+1)+T(1,j-1))*(kd*dx/(2*dy)))+(h*dy*Tinf))/((kd*dy/dx)+(kd*dx/dy)+(h*dy));
             %T(1,j)=(1/((1/dx)-(h/kd)))*((T(2,j)/dx)-(Tinf*h/kd));-ignore
         end
         %SOR through interior nodes
         %for tumour
         for i = a+1:p1x-1
             for j = b+1:p1y-1
                 T(i,j)=(1/(1+2*atx+2*aty))*(T_previous_nt(i,j)+atx*T_old(i+1,j)+atx*T(i-1,j)+aty*T_old(i,j+1)+aty*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
         end
         %other parts of the dermis
         for i = 2:p-1
             for j = 2:b-1
                 T(i,j)=(1/(1+2*a1x+2*a1y))*(T_previous_nt(i,j)+a1x*T_old(i+1,j)+a1x*T(i-1,j)+a1y*T_old(i,j+1)+a1y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
             for j = p1y+1:N-1
                 T(i,j)=(1/(1+2*a1x+2*a1y))*(T_previous_nt(i,j)+a1x*T_old(i+1,j)+a1x*T(i-1,j)+a1y*T_old(i,j+1)+a1y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
         end
         for j = b:p1y
             for i = 2:a-1
                 T(i,j)=(1/(1+2*a1x+2*a1y))*(T_previous_nt(i,j)+a1x*T_old(i+1,j)+a1x*T(i-1,j)+a1y*T_old(i,j+1)+a1y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
             for i = p1x+1:p-1
                 T(i,j)=(1/(1+2*a1x+2*a1y))*(T_previous_nt(i,j)+a1x*T_old(i+1,j)+a1x*T(i-1,j)+a1y*T_old(i,j+1)+a1y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
         end
         %for fat
         for i = p+1:g-1
             for j = 2:N-1
                 T(i,j)=(1/(1+2*a2x+2*a2y))*(T_previous_nt(i,j)+a2x*T_old(i+1,j)+a2x*T(i-1,j)+a2y*T_old(i,j+1)+a2y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
         end
         %For muscle layer
         for i = g+1:M-1
             for j = 2:N-1
                 T(i,j)=(1/(1+2*a3x+2*a3y+((wb*cb*Rob*dt)/(Rom*cm))))*((K*dt/(Rom*cm))+T_previous_nt(i,j)+a3x*T_old(i+1,j)+a3x*T(i-1,j)+a3y*T_old(i,j+1)+a3y*T(i,j-1));
                 T(i,j)=T_old(i,j)*(1-alpha)+alpha*T(i,j);
             end
         end
         
         error = max(max(abs(T_old - T)));
         T_old=T;
         sor_iter = sor_iter + 1;
     end
     %after the above loop, temperature values updated forthenext timestep
     %continuity equations
     %-at the tumour boundaries
     for j=b:p1y %horizontal walls
         T(a,j)=(1/(km+kt))*(kt*T(a+1,j)+km*T(a-1,j));
         T(p1x,j)=(1/(km+kt))*(km*T(p1x+1,j)+kt*T(p1x-1,j));
     end
     for i=a+1:p1x-1 %vertical walls
         T(i,b)=(1/(km+kt))*(kt*T(i,b+1)+km*T(i,b-1));
         T(i,p1y)=(1/(km+kt))*(km*T(i,p1y+1)+kt*T(i,p1y-1));
     end
     %-at the skin layer interfaces
     for j=2:N-1
            T(p,j)=(1/(kd+kf))*(kf*T(p+1,j)+kd*T(p-1,j));
            T(g,j)=(1/(km+kf))*(km*T(g+1,j)+kf*T(g-1,j));
     end
     %boundary condition at the bottom
     T(M,:)=37;
     %adiabatic conditions
     for i=1:M-1  
        T(i,1)=T(i,2);
        T(i,N)=T(i,N-1);
     end
    %putting temperature into a time matrix
    timeT_Surface(k,1)=T(2,floor(N/3)/2);
    timeT_muscle(k,1)=T(g,floor(N/3/2));
    error=1000;
    if T(5,12)>43.5
        break
    end
    %calculating the residue
    error1 = max(max(abs(T - T_previous_nt)));
    %updating the tempertaure value for next time steps
    T_previous_nt = T;
    %to track the number of iterations
    k=k+1;
 end