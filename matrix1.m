function K=matrix1(a,b,N)
Ke=zeros(4,4);
K=zeros(2*N+2,2*N+2);
h=(b-a)/N;
for i=1:N
  Ke(:,:,i)=h/420*[156 22*h 54 -13*h;22*h 4*h^2 13*h -3*h^2;54 13*h 156 -22*h; -13*h -3*h^2 -22*h 4*h^2];
end

for cnt=1:N
    K(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) = K(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) + ...
        Ke(:,:,cnt);
end


end