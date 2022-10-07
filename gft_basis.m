function [ v,e ] = gft_basis( G_skel,no_joints)


A=adjacency(G_skel);
D=diag(degree(G_skel));
d1=D^(-0.5);
d2=D^(-1);
L=eye(no_joints,no_joints)-d1*A*d1; %laplacian of graph(symmetric)
An=eye(no_joints,no_joints)*d2+d1*A*d1;

% [va,ea]=eig(An);  %eigen vectors of adjacency


[v,e]=eig(L);
end

