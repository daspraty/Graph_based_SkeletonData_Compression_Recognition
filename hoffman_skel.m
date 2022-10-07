function [decod_sig, bit_rate] = hoffman_skel(sig,codebook)
% codebook is a kXn sized atrix
% sig is a kXt sized matrix

[k,n]=size(codebook);
[k,t]=size(sig);
cnt=zeros(1,n);
sigma=0;
for i=1:n
    for j=1:t
        cnt(i)=cnt(i)+isequal(sig(:,j),codebook(:,i));
    end
    prb(i)=cnt(i)/t;
    sigma=sigma + prb(i);
    cumpro(i)=sigma;
end

for i=1:n
    sym{i}=strcat(num2str(codebook(1,i)),',',num2str(codebook(2,i)),',',num2str(codebook(3,i)));
end

for i=1:t
    signl{i}=strcat(num2str(sig(1,i)),',',num2str(sig(2,i)),',',num2str(sig(3,i)));
end


dict = huffmandict(sym,prb);
no_bit=0;
for i=1:length(dict)
    no_bit=no_bit+length(dict{2,2});
end
bit_rate=no_bit/length(dict);

hcode = huffmanenco(signl,dict);
dhsig = huffmandeco(hcode,dict);

% isequal(signl,dhsig)

for i=1:t
    temp=strsplit(dhsig{i},',') ;
    decod_sig(1,i)=str2num(temp{1});
    decod_sig(2,i)=str2num(temp{2});
    decod_sig(3,i)=str2num(temp{3});
end

end