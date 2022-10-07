clc;
clear all;
close all;
% source_path="/Users/pratyushadas/Desktop/usc_work/Dataset/open_pose_dataset/";
source_path="/Users/pratyushadas/Desktop/sports_output/";

folder_dir=dir(source_path);
nbfolder=length(folder_dir)-2;
window_size=30;

 joint_names={'SpineBase','SpineMid','Neck','Head','ShoulderLeft','ElbowLeft','WristLeft','HandLeft','ShoulderRight','ElbowRight','WristRight',...
   'HandRight','HipLeft','KneeLeft','AnkleLeft','FootLeft','HipRight','KneeRight','AnkleRight','FootRight','SpineShoulder','HandTipLeft',...
   'ThumbLeft','HandTipRight','ThumbRight'};
Fs=30;
G_skel_openpose=graph([1 1 1 2 2 2 2 3 4 6 7 9 10 12 13 15 16],...
    [2 15 16 3 6 9 12 4 5 7 8 10 11 13 14 17 18]);
no_hand_joint=21;
no_joints=18;
[ v,e ] = gft_basis(G_skel_openpose, no_joints  );
fs=1; % sampling frequency
%define the nonlinear quantization matrix
inc=500;
qt=ones(window_size,no_joints);
for i=1:window_size   
    qt(i,:)=inc.*qt(i,:);
    inc=inc+50;
end
for n=1:nbfolder

    folder_name=folder_dir(n+2).name
    folder_path_complete = strcat(source_path,folder_name);

    file_dir=dir(folder_path_complete);
    nbfiles=length(file_dir)-2;
    xdata=[];ydata=[];cdata=[];xdataHL=[];ydataHL=[];xdataHR=[];ydataHR=[];xdataHLR=[];ydataHLR=[];
    
    for i=1:nbfiles
        file_name=file_dir(i+2).name;
        file_path_complete = strcat(folder_path_complete,'/',file_name);

        fid = fopen(file_path_complete);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
        data=val.people(1).pose_keypoints_2d;
%         dataHL=val.people(1).hand_left_keypoints;
%         dataHR=val.people(1).hand_right_keypoints;

        p=1;
        for j=1:3:(54)
            x(p)=data(j);
            y(p)=data(j+1);
           p=p+1;
        end
%         p=1;
%         for j=1:3:(no_hand_joint*3)
%             xHL(p)=dataHL(j);
%             yHL(p)=dataHL(j+1);
%             xHR(p)=dataHR(j);
%             yHR(p)=dataHR(j+1);
% %             c(p)=data(j+2);
%             p=p+1;
%         end

        xdata=[xdata x'];
        ydata=[ydata y'];
%         xdataHL=[xdataHL xHL'];
%         ydataHL=[ydataHL yHL'];
%         xdataHR=[xdataHR xHR'];
%         ydataHR=[ydataHR yHR'];
%         cdata=[cdata c'];

    end

%     xdataHLR=[xdataHL;xdataHR];
%     ydataHLR=[ydataHL;ydataHR];

    vj_xdata=[];vj_ydata=[];xj_data_HL=[];yj_data_HL=[];vj_xdata_HL=[];vj_ydata_HL=[];xj_data_HR=[];yj_data_HR=[];vj_xdata_HR=[];vj_ydata_HR=[];
     for i=1:no_joints
         vj_xdata(i,:)=diff(xdata(i,:));
         vj_ydata(i,:)=diff(ydata(i,:));
     end
%      for i=1:no_hand_joint
%          xj_data_HL(i,:)=sample_filter( xdataHL(i,:),0);
%          yj_data_HL(i,:)=sample_filter( ydataHL(i,:),0);
%          vj_xdata_HL(i,:)=diff(xj_data_HL(i,:));
%          vj_ydata_HL(i,:)=diff(yj_data_HL(i,:));
% 
%          xj_data_HR(i,:)=sample_filter( xdataHR(i,:),0);
%          yj_data_HR(i,:)=sample_filter( ydataHR(i,:),0);
%          vj_xdata_HR(i,:)=diff(xj_data_HR(i,:));
%          vj_ydata_HR(i,:)=diff(yj_data_HR(i,:));
%      end


  
    
    %processing only the first window, use a loop to process the entire
    %data
%     for i=1:wn
    
    cf=1;
    if ~isempty(xdata) && size(xdata',1)>5
        if any(xdata)

%           [mx(ppp,:),mn(ppp,:)]=comp_gft_mxn(x_data,y_data,z_data,v);
%         [recon_sig,bit_rate]=comp(x_data,y_data,z_data,4);

        [recon_sig_gft,rs1,rs2,rs3,vx,vy,bit_rate_gft]=comp_gft(xdata(:,1:window_size)',ydata(:,1:window_size)',v,no_joints,cf,qt);
%         [s1,s2,s3]=size(recon_sig);
%         end

%%compute the mean square error
%         mse_dct(1)=mse(x_data,reshape(recon_sig(1,:,:),s2,s3));
%         mse_dct(2)=mse(y_data,reshape(recon_sig(2,:,:),s2,s3));
%         mse_dct(3)=mse(z_data,reshape(recon_sig(3,:,:),s2,s3));
%         
%         mse_gft(1)=mse(x_data,reshape(recon_sig_gft(1,:,:),s2,s3));
%         mse_gft(2)=mse(y_data,reshape(recon_sig_gft(2,:,:),s2,s3));
%         mse_gft(3)=mse(z_data,reshape(recon_sig_gft(3,:,:),s2,s3));
%         figure;bar([mse_dct;mse_gft]');legend('dct','gft')
        
        
%         figure;
%         plot(ceil(bit_rate));hold on;plot(ceil(bit_rate_gft),'r')
%         ylabel('bit rate')
%         legend('dct','gft-dct')
        
%         figure
%         subplot(3,1,1)
%         plot(x_data(:,10));hold on
%         plot(recon_sig(1,:,10),'r');hold on;
%         plot(recon_sig_gft(1,:,10),'k');
%         legend('original','dct','gft')
%         title('X')
%         subplot(3,1,2)
%         plot(y_data(:,10));hold on
%         plot(recon_sig(2,:,10),'r');hold on;
%         plot(recon_sig_gft(2,:,10),'k');
%         legend('original','dct','gft')
%         title('Y')
%         subplot(3,1,3)
%         plot(z_data(:,10));hold on
%         plot(recon_sig(3,:,10),'r');hold on;
%         plot(recon_sig_gft(3,:,10),'k');
%         legend('original','dct','gft')
%         title('Z')
        
         
%          tit2=strcat('./com_skel/',fname,'.mat');
%          save(tit2,'recon_sig')
        end
    end

end

function [xu,bit_rate]=comp(y_x,y_y,y_z,cf)
    data_dct=[];
    no_joints=size(y_x,2);
    vx=[];vy=[];vz=[];
    vx=(y_x);
    vy=(y_y);
    vz=(y_z);
    for i=1:no_joints
        data_dct(1,:,i)=dct(vx(:,i))';
        data_dct(2,:,i)=dct(vy(:,i))';
        data_dct(3,:,i)=dct(vz(:,i))';
    end
    data_lg_=round(data_dct,cf);
%     data_lg_=data_dct;
    [s1,s2,s3]=size(data_dct); 
    

    recon_sig=[];x=[];y=[];z=[];
    xu=[];
    for i=1:s1
        cw=[];temp=[];codebook=[];
        cw=reshape(data_lg_(i,:,:),[s2*s3,1]);
%         codebook=unique(cw,'rows')';
        %%vector quantization
        [idx,C]=kmeans(cw,16);
        vqenc = dsp.VectorQuantizerEncoder(...
        'Codebook', C', ...
        'CodewordOutputPort', true, ...
        'QuantizationErrorOutputPort', true);
        codebook=C';

        cwq=[];
        [ind, cwq, err] = vqenc(cw');
        [dsig,bit_rate(i)]=hoffman_skel1d(cwq,codebook);  
        dsig_=reshape(dsig,[s2,s3]);
        for j=1:no_joints
            xu(i,:,j)=idct(dsig_(:,j));          
        end
    end
     
    
end


function [recon_sig,rs1,rs2,rs3,vx,vy,bit_rate]=comp_gft(x_data,y_data,v,no_joints,cf,qt)
    tnj=18;
    vinv=inv(v);
    vx=[];vy=[];
    vx=(x_data);
    vy=(y_data);
    no_frames=size(vx,1);
    


    gft_x=v'*vx';
    gft_y=v'*vy';
    
    data_dct=[];
     
    
    gft_x=gft_x';
    gft_y=gft_y';

    no_jt=size(gft_x,2);
    for i=1:no_jt
        data_dct(1,:,i)=dct(gft_x(:,i))';
        data_dct(2,:,i)=dct(gft_y(:,i))';
        
    end
    temp(:,:)=data_dct(1,:,:);
    data_dct(1,:,:)=temp./qt;
    temp(:,:)=data_dct(2,:,:);
    data_dct(2,:,:)=temp./qt;
%     
 
    
    [s1,s2,s3]=size(data_dct); 

%     min_data_dct_x = min(data_dct(1,:,:));
%     min_data_dct_x=min(min_data_dct_x(:));
%     data_dct(1,:,:)=data_dct(1,:,:)-min_data_dct_x;
%     
%     min_data_dct_y = min(data_dct(2,:,:));
%     min_data_dct_y=min(min_data_dct_x(:));
%     data_dct(2,:,:)=data_dct(2,:,:)-min_data_dct_y;
    
    
    data_lg_=round(data_dct,3);
%     data_lg_=data_dct;
    data_all=[];
%     dsig_=data_lg_;
    
    
    recon_sig=[];x=[];y=[];rs=[];
    for i=1:s1
        cw=[];temp=[];codebook=[];
        cw=reshape(data_lg_(i,:,:),[s2*s3,1]);
%         codebook=unique(cw,'rows')';
        
        %%vector quantization
        [idx,C]=kmeans(cw,64);
        vqenc = dsp.VectorQuantizerEncoder(...
        'Codebook', C', ...
        'CodewordOutputPort', true, ...
        'QuantizationErrorOutputPort', true);
        codebook=C';

        cwq=[];
        [ind, cwq, err] = vqenc(cw');
        
        [dsig,bit_rate(i)]=hoffman_skel1d(cwq,codebook,no_frames,no_joints);  
        dsig_=reshape(dsig,[s2,s3]);
        temp(:,:)=dsig_(:,:);
        temp=temp.*qt;
        for j=1:no_joints
            rs(i,:,j)=idct(temp(:,j));          
        end
        
        
     end
    x(:,:)=rs(1,:,:);
    y(:,:)=rs(2,:,:);
    
    
    rs1=(v*x')';
    rs2=(v*y')';
    e1=norm(rs1-x_data);
    e2=norm(rs2-y_data);
    er=mean([mean(e1(:)),mean(e2(:))])
    bxy=mean(bit_rate)
    
    rs11=rs1;
    rs22=rs2;
    
%     rs11=[x_data(1,:);rs1];
%     rs22=[y_data(1,:);rs2];
%     rs33=[z_data(1,:);rs3];
%     rs11=cumsum(rs11);
%     rs22=cumsum(rs22);
%     rs33=cumsum(rs33);
%     
    recon_sig(1,:,:)=rs11;
    recon_sig(2,:,:)=rs22;

end


