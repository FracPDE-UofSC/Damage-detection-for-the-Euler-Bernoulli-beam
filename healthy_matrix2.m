function [K1]=healthy_matrix2(a,b,N,Ef,I)
k=zeros(4,4);
K1=zeros(2*N+2,2*N+2);
h=(b-a)/N;
x_points= a:h:b;
 for cnt=1:N
     for i=1:4
         for j=1:4
     xa=x_points(cnt);
     xb=x_points(cnt+1);
     px4_H1=@(x)(4*(2*x - 2*xb))/h^3 + (2*(h + 2*x - 2*xa))/h^3;
     px4_H2=@(x)(2*(2*x - 2*xb))/h^2 + (2*(x - xa))/h^2;
     px4_H3=@(x)(2*(h - 2*x + 2*xb))/h^3 - (4*(2*x - 2*xa))/h^3;
     px4_H4=@(x)(2*(2*x - 2*xa))/h^2 + (2*(x - xb))/h^2;
     eval(['f2=@(x)I*Ef(x).*px4_H',num2str(i),'(x)','.*px4_H',num2str(j),'(x)']);
     k(i,j,cnt)=quadgk(f2,xa,xb);
         end
     end
 end

% for i=1:N  
%   E_modulus= Ef((x_points(i)+x_points(i+1))/2);
%   k(:,:,i)=E_modulus*I/h^3*[12 6*h -12 6*h;6*h 4*h^2 -6*h 2*h^2;-12 -6*h 12 -6*h;6*h 2*h^2 -6*h 4*h^2];
% end


for cnt=1:N
    K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) = K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) + ...
        k(:,:,cnt);
end


end