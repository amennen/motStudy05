function [] = plotSig(xt,yt,ps,inter)
ht=1.05;ht2=1.1;ht3=1.15;
big=10;small=7;
sigs=[0.01 0.05 0.1];
dx = 1; %xt(2) - xt(1);
if ~inter %no interactions
    for pc=1:length(ps)
        if ps(pc)<sigs(1)
            plot(xt(pc)+.5*dx,max(yt)*ht,'*k','MarkerSize',big)
        elseif ps(pc)<sigs(2)
            plot(xt(pc)+.5*dx,max(yt)*ht,'dk','MarkerSize',small)
        elseif ps(pc)<sigs(3)
            plot(xt(pc)+.5*dx,max(yt)*ht,'^k','MarkerSize',small)
        end;
    end;
else %interactions too
    j=1;k=1;
    for i=1:length(ps)
        q=mod(i,3);
        if q;pst(j)=ps(i);j=j+1;
        else psi(k)=ps(i);k=k+1;
        end
    end;
    j=1;k=1;
    if length(xt)>length(ps)/3 %indicates there are > points than ROIs and we should plot these.
        for pc=1:length(ps)
            gb=floor(pc/3);
            item=pc-gb;
            regp=mod(pc,3);
            if regp
                if pst(j)<sigs(1)
                    plot(xt(j),max(yt)*ht,'*k','MarkerSize',big)
                elseif pst(j)<sigs(2)
                    plot(xt(j),max(yt)*ht,'dk','MarkerSize',small)
                elseif pst(j)<sigs(3)
                    plot(xt(j),max(yt)*ht,'^k','MarkerSize',small)
                end;
                j=j+1;
            else %interaction p
                if psi(k)<sigs(1)
                    plot(xt([item-1:item]),[1 1]*max(yt)*ht2, '-k',mean(xt([item-1:item])),max(yt)*ht3,'*b','MarkerSize',big)
                elseif psi(k)<sigs(2)
                    plot(xt([item-1:item]),[1 1]*max(yt)*ht2, '-k',mean(xt([item-1:item])),max(yt)*ht3,'db','MarkerSize',small)
                elseif psi(k)<sigs(3)
                    plot(xt([item-1:item]),[1 1]*max(yt)*ht2, '-k',mean(xt([item-1:item])),max(yt)*ht3,'^b','MarkerSize',small)
                end;
                k=k+1;
            end;
        end;
    else %indicates there are only simple comparisons to plot
        for pc=1:length(psi)
            if psi(pc)<sigs(1)
                plot(xt(pc),max(yt)*ht,'*k','MarkerSize',big)
            elseif psi(pc)<sigs(2)
                plot(xt(pc),max(yt)*ht,'dk','MarkerSize',small)
            elseif psi(pc)<sigs(3)
                plot(xt(pc),max(yt)*ht,'^k','MarkerSize',small)
            end;
        end;
    end;
end;
end