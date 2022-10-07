function [cw, C] = codebook_generation(X,k)
[idx,C]=kmeans(X,k);

% figure;
% plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
% hold on
% plot(X(idx==3,1),X(idx==3,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off
C=C';


vqenc = dsp.VectorQuantizerEncoder(...
    'Codebook', C, ...
    'CodewordOutputPort', true, ...
    'QuantizationErrorOutputPort', true, ...
    'OutputIndexDataType', 'uint8');


[ind, cw, err] = vqenc(X');
% plot(cw(1,:), cw(2,:), 'rO', x(1,:), x(2,:), 'g.');
% legend('Quantized', 'Inputs', 'location', 'best');

% mean(mean(((cw-X').^2).^(0.5)))
end