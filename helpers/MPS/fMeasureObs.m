function [expm] = fMeasureObs(A,L,R,obs,verbose)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5 || isempty(verbose),verbose=false;end

Nobs = length(obs);

%% multisite case
if iscell(A{1}) 
    
    N = length(A);
    PBC = @(n)(mod(n+N-1,N)+1);
    if isempty(L),L=cell(1,N);end
    if isempty(R),R=cell(1,N);end
    
    expm = zeros(N,Nobs);
    
    for kk=1:Nobs
        
%         NS = obs(kk).sites;
%         if NS>1
%             for nn=1:N
%                 Atmp = A{nn};
%                 for ll=1:NS-1, Atmp = concatMPS(Atmp,A{PBC(nn+ll)});end
%                 
%                 if isempty(R{PBC(nn+NS-1)}),expm(nn,kk) = trace(ApplyOpTM(Atmp,Atmp,L{PBC(nn-1)},obs(kk).op,'l'));
%                 else expm(nn,kk) = trace(ApplyOpTM(Atmp,Atmp,L{PBC(nn-1)},obs(kk).op,'l')*R{PBC(nn+NS-1)});
%                 end
%             end
%         else
%             for nn=1:N
%                 if isempty(R{nn}),expm(nn,kk) = trace(ApplyOpTM(A{nn},A{nn},L{PBC(nn-1)},obs(kk).op,'l'));
%                 else expm(nn,kk) = trace(ApplyOpTM(A{nn},A{nn},L{PBC(nn-1)},obs(kk).op,'l')*R{nn});
%                 end
%             end
%         end
        
        if obs(kk).sites==1
            for nn=1:N
                if isempty(R{nn}),expm(nn,kk) = trace(ApplyOpTM(A{nn},A{nn},L{PBC(nn-1)},obs(kk).op,'l'));
                else expm(nn,kk) = trace(ApplyOpTM(A{nn},A{nn},L{PBC(nn-1)},obs(kk).op,'l')*R{nn});
                end
            end
        elseif obs(kk).sites==2
            for nn=1:N
                Atmp = concatMPS(A{nn},A{PBC(nn+1)});
                if isempty(R{PBC(nn+1)}),expm(nn,kk) = trace(ApplyOpTM(Atmp,Atmp,L{PBC(nn-1)},obs(kk).op,'l'));
                else expm(nn,kk) = trace(ApplyOpTM(Atmp,Atmp,L{PBC(nn-1)},obs(kk).op,'l')*R{PBC(nn+1)});
                end
            end
        else error('only one or two site operators supported');
        end
    end
%% single site case
else 
    N = 1;
    expm = zeros(1,Nobs);
    
    for kk=1:Nobs
        if obs(kk).sites==1
            if isempty(R),expm(kk) = trace(ApplyOpTM(A,A,L,obs(kk).op,'l'));
            else expm(kk) = trace(ApplyOpTM(A,A,L,obs(kk).op,'l')*R);
            end
        elseif obs(kk).sites==2
            Atmp = concatMPS(A,A);
            if isempty(R),expm(kk) = trace(ApplyOpTM(Atmp,Atmp,L,obs(kk).op,'l'));
            else expm(kk) = trace(ApplyOpTM(Atmp,Atmp,L,obs(kk).op,'l')*R);
            end
        else error('only one or two site operators supported');
        end
    end
end

if verbose 
    disp({obs.name});
    disp(expm)
    if N>1
        disp('means');
        disp(mean(expm));
    end
end

end

