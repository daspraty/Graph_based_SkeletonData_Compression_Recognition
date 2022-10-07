clc;
clear all;
close all;
source_path='./skeleton_compression/sample_file/';%S017C003P020R002A060.skeleton';
output_path='./output';
file_dir=dir(source_path);
nbfiles=length(file_dir)-2;
window_size=20;

 joint_names={'SpineBase','SpineMid','Neck','Head','ShoulderLeft','ElbowLeft','WristLeft','HandLeft','ShoulderRight','ElbowRight','WristRight',...
   'HandRight','HipLeft','KneeLeft','AnkleLeft','FootLeft','HipRight','KneeRight','AnkleRight','FootRight','SpineShoulder','HandTipLeft',...
   'ThumbLeft','HandTipRight','ThumbRight'};
Fs=30;
G_skel_kinect=graph([1 1 1 2 21 21 21 3 9 10 11 12 12 5 6 7 8 8 17 18 19 13 14 15],...
    [2 13 17 21 3 5 9 4 10 11 12 24 25 6 7 8 22 23 18 19 20 14 15 16]);

no_joints=25;
[ v,e ] = gft_basis(G_skel_kinect, no_joints  );
fs=1; % sampling frequency


for s=1:nbfiles
    file_name=file_dir(s+2).name
    [pathstr,fname,ext] = fileparts(file_name);
    
    file_path_complete = strcat(source_path,file_name);
    

    
    bodyinfo = read_skeleton_file(file_path_complete);
    no_frames=length(bodyinfo);
    xr_data=[];yr_data=[];zr_data=[];
    x_data=[];y_data=[];z_data=[];
    flag=0;
    p=0;
    for k=1: no_frames
        if isempty(bodyinfo(k).bodies)
            x_data(k,1:no_joints)=0;
            y_data(k,1:no_joints)=0;
            z_data(k,1:no_joints)=0;
        else
            flag=1;p=p+1;
            for j=1:no_joints
                x_data(p,j)=bodyinfo(k).bodies(1).joints(j).x;
                y_data(p,j)=bodyinfo(k).bodies(1).joints(j).y;
                z_data(p,j)=bodyinfo(k).bodies(1).joints(j).z;
            end
        end
    end
    cf=2;
    if ~isempty(x_data) && size(x_data,1)>5
        if any(x_data)
        [recon_sig]=comp(x_data,y_data,z_data,cf);
        [recon_sig_gft,rs1,rs2,rs3]=comp_gft(x_data,y_data,z_data,1.5,v,16);
        
        joint_number=10; %plot the data of joint_number
        figure
        plot(x_data(:,joint_number));hold on
        plot(recon_sig(1,:,joint_number),'r');hold on;
        plot(recon_sig_gft(1,:,joint_number),'k');
        legend('original','dct','gft')
        
         
         tit2=strcat(output_path, '/com_skel/',fname,'.mat');
%          save(tit2,'recon_sig')
        end
    end

end

function [recon_sig]=comp(y_x,y_y,y_z,cf)

%compute coefficients only using DCT
    data_dct=[];
    no_joints=size(y_x,2);
    for i=1:no_joints
        data_dct(1,:,i)=dct(y_x(:,i))';
        data_dct(2,:,i)=dct(y_y(:,i))';
        data_dct(3,:,i)=dct(y_z(:,i))';
    end
    data_dct=round(data_dct,4);
    [s1,s2,s3]=size(data_dct); 
    nf=ceil(s2/cf);
    data_lg_=data_dct;
    data_lg_(:,s2/2+1:end,:)=0;
    data_all=[];
    data_all=reshape(data_lg_,[s1,s2*s3]);

    recon_sig=[];
    for i=1:no_joints
        cw=[];temp=[];codebook=[];
        cw(:,:)=data_lg_(:,:,i);
        codebook=unique(data_lg_(:,:,i)','rows')';
        [dsig,bit_rate(i)]=hoffman_skel(cw,codebook);  
        xu_x=[];xu_y=[];xu_z=[];
        xu_x=idct(dsig(1,:)');
        xu_y=idct(dsig(2,:)');
        xu_z=idct(dsig(3,:)');

        recon_sig(1,:,i)=xu_x;
        recon_sig(2,:,i)=xu_y;
        recon_sig(3,:,i)=xu_z;
    end
end


function [recon_sig,rs1,rs2,rs3]=comp_gft(x_data,y_data,z_data,cf,v,no_joints)\
%compute ST graph coefficients 
tnj=25;
    vinv=inv(v);
    vx=[];vy=[];vz=[];
    vx=diff(x_data);
    vy=diff(y_data);
    vz=diff(z_data);


    gft_x=v'*vx';
    gft_y=v'*vy';
    gft_z=v'*vz';
    data_dct=[];
    
    gft_x1=gft_x';
    gft_y1=gft_y';
    gft_z1=gft_z';
    
    gft_x(no_joints+1:end,:)=[];
    gft_y(no_joints+1:end,:)=[];
    gft_z(no_joints+1:end,:)=[];
    
    %compute spatial graph coefficient
    
    gft_x=gft_x';
    gft_y=gft_y';
    gft_z=gft_z';
    no_jt=size(gft_x,2);
   
    %compute temporal graph coefficient
    for i=1:no_jt
        data_dct(1,:,i)=dct(gft_x(:,i))';
        data_dct(2,:,i)=dct(gft_y(:,i))';
        data_dct(3,:,i)=dct(gft_z(:,i))';
    end
%     data_dct=round(data_dct,4);

%% Quantization and Encoding
    [s1,s2,s3]=size(data_dct); 
    nf=ceil(s2/cf);
    data_lg_=data_dct;
    data_lg_(:,nf+1:end,:)=0;
    data_all=[];
    data_all=reshape(data_lg_,[s1,s2*s3]);

    recon_sig=[];x=[];y=[];z=[];
    for i=1:no_joints
        cw=[];temp=[];codebook=[];
        cw(:,:)=data_lg_(:,:,i);
        codebook=unique(data_lg_(:,:,i)','rows')';
        [dsig,bit_rate(i)]=hoffman_skel(cw,codebook);  
        xu_x=[];xu_y=[];xu_z=[];
        
        %temporal recostruction 
        xu_x=idct(dsig(1,:)');
        xu_y=idct(dsig(2,:)');
        xu_z=idct(dsig(3,:)');
        x=[x,xu_x];
        y=[y,xu_y];
        z=[z,xu_z];
        
    end
    x=[x,zeros(s2,tnj-no_jt)];
    y=[y,zeros(s2,tnj-no_jt)];
    z=[z,zeros(s2,tnj-no_jt)];
    
    x(:,1)=gft_x(:,1);
    y(:,1)=gft_y(:,1);
    z(:,1)=gft_z(:,1);
    
    %spatial recostruction 
    rs1=(vinv'*x')';
    rs2=(vinv'*y')';
    rs3=(vinv'*z')';


    rs11=[x_data(1,:);rs1];
    rs22=[y_data(1,:);rs2];
    rs33=[z_data(1,:);rs3];
    rs11=cumsum(rs11);
    rs22=cumsum(rs22);
    rs33=cumsum(rs33);
    
    recon_sig(1,:,:)=rs11;
    recon_sig(2,:,:)=rs22;
    recon_sig(3,:,:)=rs33;
end

