function h = plotHingeMarkers(mrk,RF,col)

labels = fieldnames(mrk);
nP = numel(labels);
for k=1:nP
    
    h = plot3( mrk.(labels{k})(RF,1) , mrk.(labels{k})(RF,2),mrk.(labels{k})(RF,3),[col,'o'],'LineWidth',1 );
end
