clc;clear; format Long G

fid=fopen("Mereni.txt","r");
data=fscanf(fid,"%f %f %f %f",[4,inf])';
fclose(fid);

%% poznani polohy
r=size(data,1);
data=[data,zeros(r,1)];
m24=[]; m26=[]; m45=[]; m46=[];
for n=1:r
    if data(n,3)>200
        data(n,end)=2;
    else
        data(n,end)=1;
    end
    if data(n,1)>12400 && data(n,1)<12500
        m24=[m24;data(n,:)];
    elseif data(n,1)>12600 && data(n,1)<12700
        m26=[m26;data(n,:)];
    elseif data(n,1)>14500 && data(n,1)<14600
        m45=[m45;data(n,:)];
    else
        m46=[m46;data(n,:)];
    end
end

fid=fopen("M24.txt","w");
data=fprintf(fid,"%d    %9.4f    %9.4f    %9.3f   %d\n",m24');
fclose(fid);

fid=fopen("M26.txt","w");
data=fprintf(fid,"%d    %9.4f    %9.4f    %9.3f   %d\n",m26');
fclose(fid);

fid=fopen("M45.txt","w");
data=fprintf(fid,"%d    %9.4f    %9.4f    %9.3f   %d\n",m45');
fclose(fid);

fid=fopen("M46.txt","w");
data=fprintf(fid,"%d    %9.4f    %9.4f    %9.3f   %d\n",m46');
fclose(fid);




