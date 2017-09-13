function x=ditvec2num(vec,d)

if isempty(vec),x=1;
else
    N=length(vec);
    assert(all(vec<=d),'vec contains elements of too high value');
    tmp=d.^(N-1:-1:0);
    x=dot(vec-1,tmp)+1;
end