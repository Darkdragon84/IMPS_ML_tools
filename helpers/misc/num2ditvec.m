function vec=num2ditvec(x,d,N)

assert(x<=d^N,'value of x is too damn high!!');
if N<1
    vec=[];
else
    vec=zeros(1,N);
    x=x-1;
    div=d.^(N-1:-1:1);
    for kk=1:N-1
        vec(kk)=floor(x/div(kk))+1;
        x=mod(x,div(kk));
    end;
    vec(N)=x+1;
end;