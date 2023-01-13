function Ef=errf(x_points,scale,E,l,condition,appro)

if strcmp(condition,'fix')
    N=length(x_points)+1;
%     N=length(E);
   x_points=[0;x_points;l];
elseif strcmp(condition,'free')
    N=(length(E)+1)/2;
    x_points=[0;E(N+1:end);l];
    E=E(1:N);
end
 Ef=@(x)0;
%  E_0=1.99948e11; 
for i=1:N
        x0= x_points(i);
        x1= x_points(i+1);
        if i==1
            x0=x0-0.1;
        end
        if i==N
            x1=x1+0.1;
        end
        if strcmp(appro,'erf')
        f1=@(x)E(i).*(erf((x-x0)/scale)-erf((x-x1)/scale))/2; %error function
        elseif  strcmp(appro,'pc') &&  (i<N)
        f1=@(x)E(i).*(x >= x_points(i) ) .*(x < x_points(i+1) );
% %         %piecewise constant
        elseif strcmp(appro,'pc') &&  (i==N)
             f1=@(x)E(i).*(x >= x_points(i) ) .*(x <= x_points(i+1) );
        end
        Ef=@(x)Ef(x)+f1(x);
end

end