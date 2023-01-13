function [x,error]=dich3(numberElements,T,M0,space_point,time_gap,condition,appro,location)
% fid=fopen('result_inverse_new.txt','at');
% numberElements=60;
% E_0=3.4e9;
E_0=1.99948e11;
points(1).x=0;
points(1).y=0;
points(2).x=1;
points(2).y=0;
amplitude=1;
L=points(2).x-points(1).x;
numberNodes=numberElements+1;
h=L/numberElements;
height=0.01;
Area=0.1*height; I=0.1*height^3/12;
% rho=1220;
rho=8193.2518;
damage_ratio=0.8;
scale=sqrt(4*5e-6);
%damage vector
 E_element=ones(numberElements,1);
% % E_element(10:13)=damage_ratio; %N=16
%   E_element(5)=damage_ratio; %N=16
   E_element(6)=damage_ratio; %N=8

% E_element=ones(3,1);
% E_element(2)=damage_ratio;
%known damage location
% x_true=[1/4,sqrt(2)/4]';
% x_true=[5/8,3/4]';
%  x_true(5)=sqrt(2)/4;

x_true=(h:h:L-h)';
%  x_true(5)=sqrt(2)/4;

% x_true(5)=x_true(5)-1/16;
% x_true(6)=x_true(6)+1/16;

% x_true(5)=5/8; %length=1/16
% x_true(6)=11/16;
num=2;
omega=1.87510407^2/L^2*sqrt(E_0*I/rho/Area)/num;
% omega=4.69409113^2/L^2*sqrt(E_0*I/rho/Area)/num;
x_points=(h:h:L-h)';


if strcmp(condition,'free')
    E_element=[E_element;x_true];
%     x=[ones(numberElements,1);(h:h:L-h)'];%initial guess
%  x=[ones(3,1);[3/16,6/16]'];%initial guess
x=[ones(3,1);[5/8,11/16]'];%initial guess
end
if strcmp(condition,'fix')
    x=ones(numberElements,1);%initial guess
% x=[1,1,1]';
end
u_damaged=vobeam_damage(numberElements,E_element,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
u_healthy=vobeam_damage(numberElements,x,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
step=2*pi/M0/omega;
E_true=errf(x_true,scale,E_element,L,condition,appro);
% u_sol=vobeam_damage(numberElements,E_element,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
% u_sol_vector=u_sol(:);
rho=0.75; sg=0.25;
epi=10^(-6);
Tol=1e-14;
N=length(x);
I=eye(N);
x_results=[];
for simulate=1:50
u_sol=vobeam_damage(numberElements,E_element,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
u_sol_vector=u_sol(:);
u_sol_vector=u_sol_vector.*(ones(length(u_sol_vector),1)...
    +0.001*randn(length(u_sol_vector),1));
max=100;
J=zeros(length(u_sol_vector),N);
u=zeros(length(u_sol_vector),N+1);
error_functional=zeros(max,1);
t1=clock;

if strcmp(condition,'fix')
    x=ones(numberElements,1);%initial guess
% x=[1,1,1]';
end
if size(x_results,1)==30
    break;
end
pp=1;
for i=1:max
    flag=0;
    u_appro=vobeam_damage(numberElements,x,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
    u(:,1)=u_appro(:);
    for mmm=1:N
        u_appro1=vobeam_damage(numberElements,x+epi*I(:,mmm),x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
        u(:,mmm+1)=u_appro1(:);
    end
    for mmm=2:N+1
        J(:,mmm-1)=(u(:,mmm)-u(:,1))/epi;
    end
    r=u(:,1)-u_sol_vector;
    rhs=J'*r;
    Joc=J'*J;
    d_l=-(Joc+pp*I)\rhs;   
    for m=1:50
        [u_appro1]=vobeam_damage(numberElements,x+rho^m*d_l,x_points,T,M0,space_point,time_gap,omega,condition,appro,location);
        errd=1/2*norm( u_appro1(:)-u_sol_vector)^2 - 1/2*norm(u(:,1)-u_sol_vector)^2 -sg*rho^m*rhs'*d_l;
        if errd <=0
            break;
        end
    end % for
    if norm(rho^m*d_l)<Tol%%norm(rhs)
        break
    end
    rd=rho^m*d_l;
    x=x+rd;
    pp=pp/2;
    error_functional(i)=step*time_gap/2*norm(r)^2;
    t2=clock;
    CPU_time=etime(t2,t1);    
end
if x(end)<0 || x(end)>1.5 || x(end)<0.5 
    x=[];
end

num_iter=i-1;
x_results=[x_results;x'];
end
%% plots
TT=5.*(T==0.5);
eval(['ad=''C:/Users/liyiq/OneDrive/document/USC/work/beam/lab/HW1''']);
figname1 = [ad,'/',appro,condition,'_x_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location];
figname2 = [ad,'/',appro,condition,'_errfunc_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location];
figname3 = [ad,'/',appro,condition,'_time_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location];
figure(1)
plot(0:step:T,u_healthy(1:length(0:step:T)),'-',0:step:T,u_damaged(1:length(0:step:T)),'o','LineWidth',1.5,'MarkerSize',1.5)
% hold on
% plot(0:step:T,u_damaged(1:length(0:step:T)))
k=legend('healthy','delaminated','location','northwest');
set([k],'Interpreter','latex') 
k1=xlabel('$t$');
k2=ylabel('$w$','Rotation',0);
set([k1 k2 ],'Interpreter','latex') 

print(figname3,'-depsc2', '-r600') %eps -depsc2
% x_nodes=0:L/numberElements:L;


sensor_x=[1];
sensor_y=[1];
% sensor_x=[0.25 0.5 0.75 1];
% sensor_y=[1 1 1 1  ];
% x_location=x_points;
 x_location=x([4:end]);
%  x_location=x([20,22]);
% x_location=x(numberElements+1:end);

y_points=0.72*ones(length(x_points));
Y_points=0.73*ones(length(x_location));

%   sensor_y=[1 1 1 1 1 1 1 1 ];

Ef=errf(x_points,scale,x,L,condition,appro);
x_nodes=0:1/500:L;
y_healthy=ones(length(x_nodes),1);
figure(2)
%  X_nodes=[0:1/20:0.57,0.5701:1/400:0.65,0.6501:1/20:0.71,0.72:1/400:0.8,0.801:1/20:1];
 X_nodes=0:1/32:L;
%   X_nodes=[0:1/20:0.5,0.501:1/400:0.6,0.601:1/20:0.78,0.7801:1/400:0.85,0.8501:1/20:1];
% X_nodes=[0:1/20:0.1,0.101:1/400:0.4,0.401:1/20:0.5,0.501:1/400:0.6,0.601:1/20:0.78,0.7801:1/400:0.9,0.901:1/20:1];
% X_nodes=[0:1/100:0.5,0.501:1/400:0.6,0.601:1/100:0.78,0.7801:1/200:1];
% X_nodes=[0:1/20:0.5,0.501:1/400:0.85,0.8501:1/20:1];

   %   X_nodes=[0:1/20:0.57,0.5701:1/400:0.73,0.7301:1/20:1];
%  X_nodes=[0:1/20:0.1,0.101:1/60:0.5,0.501:1/150:0.57,0.5701:1/400:0.73,0.7301:1/20:1];

plot(x_nodes(1:end),E_true(x_nodes(1:end)),'b-',x_nodes,y_healthy',':','LineWidth',2)
hold on
plot(X_nodes(1:end),Ef(X_nodes(1:end)),'r.','LineWidth',2,'MarkerSize',10)
hold on
plot(sensor_x,sensor_y,'k.','MarkerSize',30)
hold on
plot([x_points(length(x_points)),x_points(length(x_points))],[0.7,y_points(length(x_points))],...
    'ms:','LineWidth',2,'MarkerSize',4)
hold on
plot([x_location(length(x_location)),x_location(length(x_location))],[0.7,Y_points(length(x_location))],'g-','LineWidth',2)
hold on
for i=1:length(x_points)-1
plot([x_points(i),x_points(i)],[0.7,y_points(i)],'ms:','LineWidth',2,'MarkerSize',4)
hold on
end

hold on
for i=1:length(x_location)-1
plot([x_location(i),x_location(i)],[0.7,Y_points(i)],'g-','LineWidth',2)
hold on
end
hold off

k1=xlabel('$x$');
k2=ylabel('$\kappa$','Rotation',0);
k3=legend('$\kappa_{true}(x) $','$\kappa_{healthy}(x) $','$\kappa_{pred}(x)$','sensor','FEM mesh','$\kappa_{pred}$ mesh','location','northeast');
set([k1 k2 k3],'Interpreter','latex') 
xlim([0, 1.02])
ylim([0.7 1.2])
% yTickLabel(90);
print(figname1,'-depsc2', '-r600') %eps -depsc2
%% 
x_nodes=(0:1/5000:L)';
y_points=0.72*ones(length(x_points));
Y_points=0.73*ones(length(x_location));
num=size(x_results,1);
Ef_matrix=zeros(num,length(x_nodes));
locate=zeros(num,8);
for simulate=1:num
    x=x_results(simulate,:);
    locate(simulate,:)=x_results(simulate,1:end);
    Ef=errf(x_points,scale,x',L,condition,appro);
    Ef_matrix(simulate,:)=Ef(x_nodes);
end
Ef_std=std(Ef_matrix);
Ef_mean=mean(Ef_matrix);
locate_mean=mean(locate);
locate_std=std(locate);
% x_location=locate_mean;
x_left=locate_mean-locate_std;
x_right=locate_mean+locate_std;

a=(Ef_mean-Ef_std)';
b=(Ef_mean+Ef_std)';


patch([x_nodes;flipud(x_nodes)],[a;flipud(b)],'c')
hold on
plot(x_nodes,Ef_mean,'r-',x_nodes(1:end),E_true(x_nodes(1:end)),'b-','LineWidth',1)
hold on
plot(sensor_x,sensor_y,'k.','MarkerSize',30)
% hold on
% plot([x_points(length(x_points)),x_points(length(x_points))],[0.7,y_points(length(x_points))],...
%     'ms:','LineWidth',2,'MarkerSize',4)
% hold on
% plot([x_location(length(x_location)),x_location(length(x_location))],[0.7,Y_points(length(x_location))],'g-','LineWidth',2)
% hold on

% for i=1:length(x_points)-1
% plot([x_points(i),x_points(i)],[0.7,y_points(i)],'ms:','LineWidth',2,'MarkerSize',4)
% hold on
% end
% 
% hold on
% for i=1:length(x_location)-1
% plot([x_location(i),x_location(i)],[0.7,Y_points(i)],'g-','LineWidth',2)
% hold on
% end
% hold off
xlim([0 ,1])
k1=xlabel('$x$');
k2=ylabel('$\kappa$','Rotation',0);
k3=legend('$\kappa_{std}$','$\kappa_{mean} $','$\kappa_{true}$',...
    'sensor','FEM mesh','$\kappa_{pred}$ mesh','location','southwest');
set([k1 k2 k3],'Interpreter','latex') 

xlim([0, 1.02])
ylim([0.5 1.2])

figure(3)
semilogy(error_functional(1:num_iter),'-*')
xlabel('iteration')
k1=ylabel(' $\mathcal{F}$','Rotation',0)
set([k1 ],'Interpreter','latex') 

print(figname2,'-depsc2', '-r600') %eps -depsc2
% plot(errorL2L2(:),'-*')
% print(figname3,'-depsc2', '-r600') %eps -depsc2
error_functional=error_functional(1:num_iter);
eval(['ad1=''C:/Users/liyiq/Documents/MATLAB/CODE/beam/euler/Euler Bernoulli/delamination_server/data''']);
figname4 = [ad1,'/',appro,condition,'_x_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location];

mat_name=strcat(figname4,'.mat');

% eval(['save', ' ', appro,condition,'_x_',num2str(numberElements),'_',num2str(TT),...
%     '_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',...
%     num2str(num),'_',location,'.mat',' x','error_functional(1:num_iter)']);

% eval([appro,condition,'_x_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'=x']);
% eval([appro,condition,'_errfunc_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'=error_functional(1:num_iter)']);
% eval(['save', ' ', appro,condition,'_x_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'.mat']);
% eval([' save', ' ',appro, condition,'errfunc_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'.mat']);
% eval([appro,condition,'_time_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'=x']);
% eval(['save', ' ',appro,condition,'_time_',num2str(numberElements),'_',num2str(TT),'_',num2str(M0),'_',num2str(space_point),'_',num2str(time_gap),'_omega',num2str(num),'_',location,'.mat']);
save(mat_name,'x','error_functional','u_healthy','u_damaged')
end



