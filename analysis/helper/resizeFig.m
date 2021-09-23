function resizeFig(hF,t,axs)

function chSize(~,~)
    try
        t.Position(3)=t.Parent.Position(3);
        t.Position(4)=t.Extent(4);
        t.Position(1:2)=[5 t.Parent.Position(4)-t.Position(4)];
        drawnow;
        
        if nargin ==3 && ~isempty(axs)
        
            for kk=1:length(axs)
                axs(kk).Units='pixels';
                axs(kk).Position(2)=55;
                axs(kk).Position(4)=...
                    (axs(kk).Parent.Position(4)-25-...
                    t.Position(4))-axs(kk).Position(2);
                axs(kk).Units='normalized'; 
            end
        end
    end 
end
chSize;

hF.SizeChangedFcn=@chSize;
end

