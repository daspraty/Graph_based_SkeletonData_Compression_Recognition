function [ output_mean1 ] = box_mean_pulling( data )
L=size(data,2);
L1=floor(size(data,2)/2);
L2=floor(size(data,2)/3);
L3=floor(size(data,2)/4);
L4=floor(size(data,2)/5);
%mean feature value blockwise

    output_mean=[mean(data,2),...
        mean(data(:,1:L1),2), mean(data(:,L1+1:end),2),...
        mean(data(:,1:L2),2),mean(data(:,L2+1:2*L2),2), mean(data(:,2*L2+1:end),2), ...
        mean(data(:,1:L3),2),mean(data(:,L3+1:2*L3),2), mean(data(:,2*L3+1:3*L3),2), mean(data(:,3*L3+1:end),2),...
        mean(data(:,1:L4),2),mean(data(:,L4+1:2*L4),2), mean(data(:,2*L4+1:3*L4),2), mean(data(:,3*L4+1:4*L4),2), mean(data(:,4*L4+1:end),2)];
    
    output_mean1=reshape(output_mean,[],1);
end

