function A = randMPS(d,ml,mr,cmplx,dir)
% dir specifies if A is a row or column cell
if nargin<4||isempty(cmplx),cmplx=true;end
if nargin<5||isempty(dir),dir='l';end

if dir=='l',A=cell(d,1);
elseif dir=='r',A=cell(1,d);
else error('wrong direction specified');
end

nrm=sqrt(d*sqrt(ml*mr));
if cmplx
    for kk=1:d
        A{kk} = (randn(ml,mr) + 1i*randn(ml,mr))/nrm;
    end
else
    for kk=1:d
        A{kk} = randn(ml,mr)/nrm;
    end
end
