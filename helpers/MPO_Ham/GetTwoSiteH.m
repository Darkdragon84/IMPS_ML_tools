function H=GetTwoSiteH(param,d)

assert(length(param)==5 || length(param)==7,'# of parameters must be 5 or 7');

[sx,~,sz,sp,sm] = su2gen(d);
% switch d
%     case 2 % S=1/2
%         sx=0.5*[0,1;1,0];
%         sp=[0,1;0,0];
%         sm=[0,0;1,0];
%         sz=0.5*[1,0;0,-1];
%   
%     case 3 % S=1
%         sx=[0,1,0;1,0,1;0,1,0]/sqrt(2);
%         sp=sqrt(2)*[0,1,0;0,0,1;0,0,0];
%         sm=sqrt(2)*[0,0,0;1,0,0;0,1,0];
%         sz=[1,0,0;0,0,0;0,0,-1];
%         
%     otherwise
%         error('dimension not implemented');
% end;

unity=eye(d);


    
if length(param)==5
    hx1=param(4);
    hx2=param(4);
    hz1=param(5);
    hz2=param(5);
else
    hx1=param(4);
    hx2=param(5);
    hz1=param(6);
    hz2=param(7);
end;
    
    
    Jx=param(1);
    Jy=param(2);
    Jz=param(3);
    
    Jp=Jx+Jy;
    Jm=Jx-Jy;
    
    H=-0.25*(Jp*(kron(sp,sm)+kron(sm,sp)) + Jm*(kron(sp,sp) + kron(sm,sm))) - Jz*kron(sz,sz)...
        -0.5*(hx1*kron(sx,unity) + hx2*kron(unity,sx) + hz1*kron(sz,unity) + hz2*kron(unity,sz));
% end;