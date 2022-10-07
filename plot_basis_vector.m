function []= plot_basis_vector(G_skel1, v,e)
th=.001;m=25;
k=3;l=2;no_pl=25;
rset=[];bset=[];gset=[];
e1=diag(e);

h=figure('units','normalized','outerposition',[0 0 1 1])
% if z==25
%     k=5;l=5;
%     no_pl=25;
% end
no_pl=6;
pl=1:6;
for k1=1:no_pl
    j=pl(k1);
    rset=[];bset=[];gset=[];p1=0;p2=0;p3=0;
    for i=1:m
        if abs(v(i,j))<th
            nodelabel{i}= 'o';
            p1=p1+1;bset(p1)=i;
        elseif (v(i,j))<0
            nodelabel{i}= '-';
            p2=p2+1;rset(p2)=i;
        elseif (v(i,j))>0
            nodelabel{i}= '+';
           p3=p3+1; gset(p3)=i;
        end
    end
   
    subplot(l,k,k1)
%     subplot(l,k,j)
% h=plot(G_skel1,'Layout','force','LineWidth',3);
    h=plot(G_skel1,'LineWidth',3);
%     h.NodeFontSize = 20;
%     h.NodeLabel=nodelabel;
    highlight(h,rset,'NodeColor','r','MarkerSize',7)
    highlight(h,gset,'NodeColor','g','MarkerSize',7)
    highlight(h,bset,'NodeColor','b','MarkerSize',7)
    h.NodeLabel={};
    h.XData= [0 0 0 0 -1 -1.5 -1 -1.2 1 1.5 1 1.2 -.8 -0.8 -0.8 -1 .8 .8 .8 1 0 -1.3 -1.1 1.3 1.1];
    h.YData=[0 0.3 0.8 1 0.7 0.3 0.1 0 0.7 0.3 0.1 0 -.2 -.5 -.9 -1.1 -.2 -.5 -.9 -1.1 0.5 -0.2 -.1 -0.2 -.1]; 
    h.ZData= zeros(1,25);
    
%     h.XData= 1:1:25;
%     h.YData= zeros(1,25);
%     h.ZData= zeros(1,25);
    
%     set(gca,'LooseInset',get(gca,'TightInset'));
    t=strcat('{\lambda}=',num2str(e1(j)));
    title(t,'FontSize',30)
end
% nam=strcat('basis_linegraph25.png');
% nam=strcat('basis_skel_kinectv2.png');
% saveas(h,nam)
    

end