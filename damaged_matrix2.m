function [K1]=damaged_matrix2(a,b,N)
k=zeros(4,4);
K1=zeros(2*N+2,2*N+2);
h=(b-a)/N;
for i=1:N
  k(:,:,i)=1/h^3*[12 6*h -12 6*h;6*h 4*h^2 -6*h 2*h^2;-12 -6*h 12 -6*h;6*h 2*h^2 -6*h 4*h^2];
end
for cnt=1:47 
    K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) = K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) + ...
        k(:,:,cnt);
end

for cnt=49:N
    K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) = K1(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) + ...
        k(:,:,cnt);
end

    K1(2*48-1:2*48+2,2*48-1:2*48+2) = K1(2*48-1:2*48+2,2*48-1:2*48+2) + ...
       0.9* k(:,:,48);



end