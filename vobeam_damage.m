% v_tt+\lambda v_t +\lambda v +\lambda D^{\alpha(t)} v =g
% v(0)=v_0, v^{\prime}(0)=\check{v}_0
% function [u]=vobeam(T,M,numberElements,damping,f,condition,m1,bd,alpha,h,omega)
%% parameter of beam
function u_data=vobeam_damage(numberElements,E_element,x_points,T,M0,space_point,...
    time_gap,omega,condition,appro,location)
format long
r=1;
bd='CT';
m1=1;
height=0.01;
T_init=0;
E_0=1.99948e11; 
points(1).x=0;
points(1).y=0;
points(2).x=1;
points(2).y=0;
amplitude=1;
L=points(2).x-points(1).x;
% E_element=E_0*E_element;
scale=sqrt(4*5e-6);
Ef=errf(x_points,scale,E_element,L,condition,appro);  
Ef=@(x)Ef(x).*E_0;
damp=[0.8,1];
numberNodes=numberElements+1;
h=L/numberElements;
x_nodes=0:h:L;
GDof=2*numberNodes;
damping=0;
z=0;
Area=0.1*height; I=0.1*height^3/12;
rho=8193.2518;
%% SS
if bd=='SS'
    u0=@(x)sin(m1*pi*x/L);%need to be chosen accordingly
    u0x=@(x)m1*pi*cos(m1*pi*x/L)/L;
    beta=m1*pi/L;
end
%% CC
if bd=='CC'
    if m1<=5
        Beta=[4.73004074,7.85320462,10.9956079,14.1371655 17.2787597];
        beta=Beta(m1);
    end
    if m1>5
        beta=(2*m1+1)*pi/2;
    end
    Sigma=(cosh(beta*L)-cos(beta*L))./(sinh(beta*L)-sin(beta*L));
    u0=@(x)cosh(beta*x)-cos(beta*x)-Sigma*(sinh(beta*x)-sin(beta*x));
    u0x=@(x)beta*(sinh(beta*x)+sin(beta*x)-Sigma*(cosh(beta*x)-cos(beta*x)));
    u0xxxx=@(x)beta^4*(cosh(beta*x)-cos(beta*x)-Sigma*(sinh(beta*x)-sin(beta*x)));
end

if bd=='CT'
    if m1<=5
        Beta=[1.87510407 4.69409113 7.85475744 10.99554073 14.13716839];
        Sigma=[0.7341 1.0185 0.9992 1.0000 1.0000];
        beta=Beta(m1);
        sigma=Sigma(m1);
    else
        beta=(2*m1-1)*pi/2;
        sigma=1;
    end
    u0=@(x)cosh(beta*x)-cos(beta*x)-sigma*(sinh(beta*x)-sin(beta*x));
    u0x=@(x)beta*(sinh(beta*x)+sin(beta*x)-sigma*(cosh(beta*x)-cos(beta*x)));
    u0xxx=@(x)beta^3*(sinh(beta*x)-sin(beta*x)-sigma*(cosh(beta*x)+cos(beta*x)));
end
step=2*pi/M0/omega;
M=ceil((T-T_init)/step);
tau=step;
%% right hand side
q=@(t)amplitude*cos(omega*t);
%%  BC
if bd=='SS'
    %%pined-pinned SSBC
    fixedNodeU =[1 2*numberElements+1]' ;
    fixedNodeV=[];
end
if bd=='CC'
    %%fixed-fixed'
    fixedNodeU =[1 2*numberElements+1]' ;
    fixedNodeV=[2 2*numberElements+2]';
end

if bd=='CT'
    %%pined-pinned SSBC
    fixedNodeU =[1]' ;
    fixedNodeV=[2]';
end
%
prescribedDof=[fixedNodeU;fixedNodeV];
activeDof=setdiff([1:GDof]',[prescribedDof]);
% c=zeros(M+1,1);
%% FEM matrix
K2= matrix1(points(1).x,points(2).x,numberElements);%transverse mass
K21=rho*Area*K2;%transverse mass
K221=damping*rho*Area*K2;
K22= healthy_matrix2(points(1).x,points(2).x,numberElements,Ef,I);%tranverse stiffness
K21=K21(activeDof,activeDof);
K221=K221(activeDof,activeDof);
K22=K22(activeDof,activeDof);
U=zeros(length(activeDof),M+1);
%% SS
if bd=='SS'
    d4=(m1*pi)^4/L^4;
end
if bd=='CC'
    d4=beta^4;
end
if bd=='CT'
    d4=beta^4;
end
n0=2;
rhsb=vecload(@(x)0.5,0,L,numberElements,'uniform');

for n=n0:M
    %tn=time(n);
    tn=T_init+(n-1)*step;
    %% central difference
    Q=(q(tn))*rhsb;
    Q=Q(activeDof,1);
    rhs=tau^2*Q+2*(K21)*U(:,n)-(K21)*U(:,n-1)+ tau*K221*U(:,n)-tau^2*K22/2*U(:,n)-tau^2*K22/4*U(:,n-1);
    A=K21+tau*K221+tau^2*K22/4;
     U(:,n+1)=A\rhs;
end

u=zeros(2*numberElements+2,M+1);
u(activeDof,:)=U;
u_data=zeros(numberElements/space_point,floor(M/time_gap)+1);
if space_point<numberElements
 u_data(1:end,:)=u(1+2*space_point:2*space_point:2*numberElements+1,1:time_gap:end);
elseif space_point==numberElements && strcmp(location,'left')
u_data(1:end,:)=u(3,1:time_gap:end);
elseif space_point==numberElements && strcmp(location,'right')
    u_data(1:end,:)=u(2*numberElements+1,1:time_gap:end);
end
end