function F=vecload(q,a,b,N,type)
L=b-a;
h=L/N;
F=zeros(2*N+2,1);
x_nodes=0:h:L;
if type=='uniform'
    for cnt=1:N
        xa=x_nodes(cnt);
        xb=x_nodes(cnt+1);
        H1=@(x)(-(x-xb).^2.*(-h+2*(xa-x)))/h^3;
        H2=@(x)(x-xa).*(x-xb).^2/h^2;
        H3=@(x)(x-xa).^2.*(h+2*(xb-x))/h^3;
        H4=@(x)(x-xa).^2.*(x-xb)/h^2;
        f1=@(x)H1(x).*q(x);
        f2=@(x)H2(x).*q(x);
        f3=@(x)H3(x).*q(x);
        f4=@(x)H4(x).*q(x);
        fe(1,1)=quadgk(f1,xa,xb);
        fe(2,1)=quadgk(f2,xa,xb);
        fe(3,1)=quadgk(f3,xa,xb);
        fe(4,1)=quadgk(f4,xa,xb);
        F(2*cnt-1:2*cnt+2)=F(2*cnt-1:2*cnt+2)+fe;
    end
end
end


