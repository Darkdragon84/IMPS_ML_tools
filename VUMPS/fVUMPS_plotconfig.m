function [fh,ah,lh] = fVUMPS_plotconfig(params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = params.N;
if params.plotlam
    if isfield(params,'ahlam')
        ah.ahlam = params.ahlam;
        fh.fhlam = get(ah.ahlam,'parent');
    else
        fh.fhlam = figure;
        ah.ahlam = axes('yscale','log','parent',fh.fhlam);
        xlabel('# of Schmidt Value \lambda');
        ylabel('log(\lambda)');
        title('Schmidt Values');
    end
    
    lh.lhlam = zeros(1,N);
    lstr = cell(1,N);
    for nn=1:N
        lh.lhlam(nn) = line(0,0,'linestyle','none','marker','.','color',[nn/N,1-nn/N,0],'parent',ah.ahlam);
        lstr{nn} = ['\lambda_{',int2str(nn),'}'];
    end
    legend(lh.lhlam,lstr,'location','best');
end

if params.plotnorm
    if isfield(params,'ahnrm')
        ah.ahnrm = params.ahnrm;
        fh.fhnrm = get(ah.ahnrm,'parent');
    else
        fh.fhnrm = figure;
        ah.ahnrm = axes('yscale','log','parent',fh.fhnrm);
        if params.plotvst,xlabel('t [s]');
        else xlabel('iterations');
        end
        ylabel('|B|','rotation',0);
        title('Gradient Norm');
    end
    lh.lhnrm = line(1,1,'marker','.','color','k','parent',ah.ahnrm);
    drawnow;
end

if params.plotex
    if isfield(params,'ahex')
        ah.ahex = params.ahex;
        fh.fhex = get(ah.ahex,'parent');
    else
        fh.fhex = figure;
        ah.ahex = axes('yscale','log','parent',fh.fhex);
        if params.plotvst,xlabel('t [s]');
        else xlabel('iterations');
        end
        ylabel('\Deltae','rotation',0);
        title('Energy Error');
    end
    lh.lhex = line(1,1,'marker','.','color','b','parent',ah.ahex);
end

if params.plotdlam
    if isfield(params,'ahdlam')
        ah.ahdlam = params.ahdlam;
        fh.fhdlam = get(ah.ahdlam,'parent');
    else
        fh.fhdlam = figure;
        ah.ahdlam = axes('yscale','log','parent',fh.fhdlam);
        if params.plotvst,xlabel('t [s]');
        else xlabel('k');
        end
        ylabel('|\Delta\lambda|','rotation',0);
        title('Schmidt Value Changes');
    end
    lh.lhdlam = line(1,1,'marker','.','color','k','parent',ah.ahdlam);
    drawnow;
end

if params.plotxi
    if isfield(params,'ahxi')
        ah.ahxi = params.ahxi;
        fh.fhxi = get(ah.ahxi,'parent');
    else
        fh.fhxi = figure;
        ah.ahxi = axes('parent',fh.fhxi);
        if params.plotvst,xlabel('t [s]');
        else xlabel('k');
        end
        ylabel('\xi','rotation',0);
        title('Correlation Lengths');
    end
    lh.lhxi = zeros(params.nxi,1);
    for kk=1:params.nxi
        lh.lhxi(kk) = line(1,1,'marker','.','color','k','parent',ah.ahxi);
    end
    drawnow;
end

if exist('fh','var')~=1,fh=[];end
if exist('ah','var')~=1,ah=[];end
if exist('lh','var')~=1,lh=[];end


end

