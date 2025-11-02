%%
%T1

sub=dir('D:\NEW\ADNI');
sub(1:2)=[];
T1={};
num=0 
    folders={sub_sub.name};
    ifcontains=cellfun(@(folder) contains(folder,'MPRAGE_SENSE2'),folders);
    T1{i,1}=folders(ifcontains);
    temp=folders(ifcontains);
    for j=1:1
        files=dir([sub(i).folder,'\',sub(i).name,'\',temp{j}]);
        files(1:2)=[];
        for k=1:length(files)
            num=num+1;
            temp1=[sub(i).name,'_',files(k).name];
            T1all{num,1}=temp1;
            [a1,b1]=ismember(temp1,fmrisub);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if a1
                mkdir(['J:\zhengzihao_ELMCI_20241120\NC\resting_nii\T1\',fmrisubnew{b1}]);
                copyfile([sub(i).folder,'\',sub(i).name,'\',temp{j},'\',files(k).name,'\*'],...
                    ['J:\zhengzihao_ELMCI_20241120\NC\resting_nii\T1\',fmrisubnew{b1}]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        end
        
    end
    
end



sub=dir('J:\zhengzihao_ELMCIT1all\ADNI');
sub(1:2)=[];
T1={};
num=0
for i=1:length(sub)
    sub_sub=dir([sub(i).folder,'\',sub(i).name]);
    sub_sub(1:2)=[];
    folders={sub_sub.name};
    ifcontains=cellfun(@(folder) contains(folder,'MPRAGE'),folders);
    T1{i,1}=folders(ifcontains);
    temp=folders(ifcontains);
    if ~isempty(temp)
    for j=1:1
        files=dir([sub(i).folder,'\',sub(i).name,'\',temp{j}]);
        files(1:2)=[];
        for k=1:length(files)
            num=num+1;
            temp1=[sub(i).name,'_',files(k).name];
            T1all{num,1}=temp1;
  

                mkdir(['J:\ELMCI_T1\',temp1]);
                copyfile([sub(i).folder,'\',sub(i).name,'\',temp{j},'\',files(k).name,'\*'],...
                    ['J:\zhengzihao_ELMCI_20241120\NC\resting_nii\T1\',temp1]);

        end
        
    end
    end
    
end

sub=dir('D:\NEW\ADNI');
sub(1:2)=[];
T1={};
num=0;
for i=1:length(sub)
    sub_sub=dir([sub(i).folder,'\',sub(i).name]);
    sub_sub(1:2)=[];
    folders={sub_sub.name};
    ifcontains=cellfun(@(folder) contains(folder,'MPRAGE_SENSE2'),folders);
    T1{i,1}=folders(ifcontains);
    temp=folders(ifcontains);
    if ~isempty(temp)
    for j=1:length(temp)
        files=dir([sub(i).folder,'\',sub(i).name,'\',temp{j}]);
        files(1:2)=[];
        for k=1:length(files)
            num=num+1;
            temp1=[sub(i).name,'_',files(k).name];
            T1all{num,1}=temp1;
  

                mkdir(['J:\zhengzihao_ELMCI_20241120\NC\resting_nii\T1\',temp1]);
                copyfile([sub(i).folder,'\',sub(i).name,'\',temp{j},'\',files(k).name,'\*'],...
                    ['J:\zhengzihao_ELMCI_20241120\NC\resting_nii\T1\',temp1]);

        end
        
    end
    end
    
end







pathdir1=dir('J:\MCI_fmri');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    fmrisub{i,1}=[pathdir1(i).name(1:10),'_',s{5},'_',s{6},'_',s{7},'_',s{8}];
    
end
fmrisubnew={pathdir1.name}';
[a,b]=ismember(fmrisub,T1all);

%%

load('J:\adni.mat');
NC_excel(:,1)=[];

MCI_NC_excel=cat(1,MCI_excel,NC_excel);

for i=1:length(MCI_NC_excel)
    s1=MCI_NC_excel{i,1};
    s2=MCI_NC_excel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
    if i>259 & str2num(s2{2})<10;
        s2{2}=['0',s2{2}];
    end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    MCI_NC_COV{i,1}=[s1,'_',s2new];
    if  strcmp(MCI_NC_excel{i,2},'MCI')
        MCI_NC_COV{i,2}=1;
    else
        MCI_NC_COV{i,2}=-1;
    end
    
    if  MCI_NC_excel{i,3}=='M'
        MCI_NC_COV{i,3}=1;
    else
        MCI_NC_COV{i,3}=2;
    end
    MCI_NC_COV{i,4}=MCI_NC_excel{i,4};
    
    
    
end


MCI_head_path='J:\MCI_fmri_nii\preprocess\headmotion';
MCI_head_files=dir(MCI_head_path);
MCI_head_files(1:2)=[];

NC_head_path='J:\NC_fmri_nii\preprocess\headmotion';
NC_head_files=dir(NC_head_path);
NC_head_files(1:2)=[];
num=0;
MCI_NC_COV_new={};
for i=1:length(MCI_head_files)
    
    
    num=num+1;
    cd([MCI_head_path,'\',MCI_head_files(i).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];
    meanFD(num)=mean(FD);
    s=strsplit(MCI_head_files(i).name,'_');
    
    sub=[MCI_head_files(i).name(1:10),'_',s{5}];
    [id1,id2]=ismember(sub,MCI_NC_COV(:,1));
    MCI_NC_COV_new{num,1}=sub;
    MCI_NC_COV_new{num,2}=MCI_NC_COV{id2,2};
    MCI_NC_COV_new{num,3}=MCI_NC_COV{id2,3};
    MCI_NC_COV_new{num,4}=MCI_NC_COV{id2,4};
    MCI_NC_COV_new{num,5}=mean(FD);
    
end

for i=1:length(NC_head_files)
    num=num+1;
    cd([NC_head_path,'\',NC_head_files(i).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);
  
    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];
    meanFD(num)=mean(FD);
    
    s=strsplit(NC_head_files(i).name,'_');
    
    sub=[NC_head_files(i).name(1:10),'_',s{5}];
    [id1,id2]=ismember(sub,MCI_NC_COV(:,1));
    MCI_NC_COV_new{num,1}=sub;
    MCI_NC_COV_new{num,2}=MCI_NC_COV{id2,2};
    MCI_NC_COV_new{num,3}=MCI_NC_COV{id2,3};
    MCI_NC_COV_new{num,4}=MCI_NC_COV{id2,4};
    MCI_NC_COV_new{num,5}=mean(FD);
    
end
% MCI_exclude=find(meanFD>0.25);
%%
load('J:\adni.mat');
load('J:\adni_nc_new.mat');
NC_excel1=NC_excel_newCOV1;

% MCI_NC_excel=cat(1,MCI_excel,NC_excel1);
% 
% for i=1:length(MCI_NC_excel)
%     s1=MCI_NC_excel{i,1};
%     s2=MCI_NC_excel{i,9};
%     s2=s2(2:end);
%     s2=strsplit(s2,'-');
%     if str2num(s2{1})<10;
%         s2{1}=['0',s2{1}];
%     end
%     
% %     if i>259 & str2num(s2{2})<10;
% %         s2{2}=['0',s2{2}];
% %     end
%     s2new=[s2{3},'-',s2{1},'-',s2{2}];
%     MCI_NC_COV{i,1}=[s1,'_',s2new];
%     if  strcmp(MCI_NC_excel{i,2},'MCI')
%         MCI_NC_COV{i,2}=1;
%     else
%         MCI_NC_COV{i,2}=-1;
%     end
%     
%     if  MCI_NC_excel{i,3}=='M'
%         MCI_NC_COV{i,3}=1;
%     else
%         MCI_NC_COV{i,3}=2;
%     end
%     MCI_NC_COV{i,4}=MCI_NC_excel{i,4};
%     
%     
%     
% end


% MCI_head_path='J:\MCI_fmri_nii\preprocess\headmotion';
% MCI_head_files=dir(MCI_head_path);
% MCI_head_files(1:2)=[];
% 
% NC_head_path='J:\NC_fmri_nii\preprocess\headmotion';
% NC_head_files=dir(NC_head_path);
% NC_head_files(1:2)=[];
% num=0;
% MCI_NC_COV_new={};
% for i=1:length(MCI_head_files)
%     
%     
%     num=num+1;
%     cd([MCI_head_path,'\',MCI_head_files(i).name]);
%     head_file=dir(['rp*']);
%     rp1=importdata([head_file.name]);
% 
%     rp = rp1;
%     rp(:,4:6) = rp1(:,4:6)*180/pi;
%     rpMax = max(abs(rp));
%     rpDiff = diff(rp1);
%     rpDiff(:,4:6) = rpDiff(:,4:6)*50;
%     FD = [0;sum(abs(rpDiff),2)];
%     meanFD(num)=mean(FD);
%     s=strsplit(MCI_head_files(i).name,'_');
%     
%     sub=[MCI_head_files(i).name(1:10),'_',s{5}];
%     [id1,id2]=ismember(sub,MCI_NC_COV(:,1));
%     MCI_NC_COV_new{num,1}=sub;
%     MCI_NC_COV_new{num,2}=MCI_NC_COV{id2,2};
%     MCI_NC_COV_new{num,3}=MCI_NC_COV{id2,3};
%     MCI_NC_COV_new{num,4}=MCI_NC_COV{id2,4};
%     MCI_NC_COV_new{num,5}=mean(FD);
%     
% end
% 
% for i=1:length(NC_head_files)
%     num=num+1;
%     cd([NC_head_path,'\',NC_head_files(i).name]);
%     head_file=dir(['rp*']);
%     rp1=importdata([head_file.name]);
%   
%     rp = rp1;
%     rp(:,4:6) = rp1(:,4:6)*180/pi;
%     rpMax = max(abs(rp));
%     rpDiff = diff(rp1);
%     rpDiff(:,4:6) = rpDiff(:,4:6)*50;
%     FD = [0;sum(abs(rpDiff),2)];
%     meanFD(num)=mean(FD);
%     
%     s=strsplit(NC_head_files(i).name,'_');
%     
%     sub=[NC_head_files(i).name(1:10),'_',s{5}];
%     [id1,id2]=ismember(sub,MCI_NC_COV(:,1));
%     MCI_NC_COV_new{num,1}=sub;
%     MCI_NC_COV_new{num,2}=MCI_NC_COV{id2,2};
%     MCI_NC_COV_new{num,3}=MCI_NC_COV{id2,3};
%     MCI_NC_COV_new{num,4}=MCI_NC_COV{id2,4};
%     MCI_NC_COV_new{num,5}=mean(FD);
%     
% end



%%

age_MCI=cell2mat(MCI_excel(:,4));
[temp,id1]=sort(age_MCI);

age_NC=cell2mat(NC_excel(:,4));
[temp1,id2]=sort(age_NC);

age_NC1=cell2mat(NC_excel1(:,4));
[temp2,id3]=sort(age_NC1);


NC_index1=randsample(1:221,115);
NC_index2=randsample(1:261,129);


MCI_excel_sort=MCI_excel(id1,:);
NC_excel_sort=NC_excel(id2,:);
NC_excel1_sort=NC_excel1(id3,:);

NC_excelnew=cat(1,NC_excel_sort(:,:),NC_excel1_sort(:,:));



MCI_NC_excel=cat(1,MCI_excel_sort,NC_excelnew);
for i=1:length(MCI_NC_excel)
    s1=MCI_NC_excel{i,1};
    s2=MCI_NC_excel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
    if i>259 & i<489 & str2num(s2{2})<10;
        s2{2}=['0',s2{2}];
    end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    MCI_NC_COV{i,1}=[s1,'_',s2new];
    if  strcmp(MCI_NC_excel{i,2},'MCI')
        MCI_NC_COV{i,2}=1;
    else
        MCI_NC_COV{i,2}=-1;
    end
    
    if  MCI_NC_excel{i,3}=='M'
        MCI_NC_COV{i,3}=1;
    else
        MCI_NC_COV{i,3}=2;
    end
    MCI_NC_COV{i,4}=MCI_NC_excel{i,4};
    
    
    
end


pathdir1=dir('I:\ADNI\MCI\Spreading_Smooth_fc_fdr');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    MCIsub{i,1}=[pathdir1(i).name(1:10),'_',s{5}];
    
end

pathdir1=dir('I:\ADNI\NC\BOLD_Smooth_new');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    NCsub1{i,1}=[pathdir1(i).name(1:10),'_',s{5}];
    
end

pathdir1=dir('I:\ADNI\NC_NEW\BOLD_Smooth_new');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    NCsub2{i,1}=[pathdir1(i).name(1:10),'_',s{5}];
    
end

[a1,b1]=ismember(MCIsub,MCI_NC_COV(:,1));
[a2,b2]=ismember(NCsub1,MCI_NC_COV(:,1));
[a3,b3]=ismember(NCsub2,MCI_NC_COV(:,1));

pathdir1=dir('I:\ADNI\NC\BOLD_Smooth_new');
pathdir1(1:2)=[];
for i=1:length(a2)
    if a2(i)
    copyfile([pathdir1(i).folder,'\',pathdir1(i).name],['I:\ADNI\NC_All_BOLD_Smooth_new'])
    end
    
end

pathdir1=dir('I:\ADNI\NC_NEW\BOLD_Smooth_new');
pathdir1(1:2)=[];
for i=1:length(a3)
    if a3(i)
    copyfile([pathdir1(i).folder,'\',pathdir1(i).name],['I:\ADNI\NC_All_BOLD_Smooth_new'])
    end
    
end



pathdir1=dir('I:\ADNI\MCI\Spreadingnew');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    MCIsub{i,1}=[pathdir1(i).name(1:10),'_',s{5}];
    
end

pathdir1=dir('I:\ADNI\NC_allnew');
pathdir1(1:2)=[];
for i=1:length(pathdir1)
    s=strsplit(pathdir1(i).name,'_');
    
    NCsuball{i,1}=[pathdir1(i).name(1:10),'_',s{5}];
    
end

[a1,b1]=ismember(MCIsub,MCI_NC_COV(:,1));
[a2,b2]=ismember(NCsuball,MCI_NC_COV(:,1));

MCI_NC_COV_new=MCI_NC_COV([b1;b2],:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a1,b1]=ismember(MCIsub,MCI_NC_COV_new(:,1));
[a2,b2]=ismember(NCsub1,MCI_NC_COV_new(:,1));
[a3,b3]=ismember(NCsub2,MCI_NC_COV_new(:,1));

pathdir1=dir('J:\MCI_T1_nii\T1');
pathdir1(1:2)=[];
for i=1:length(a1)
    if a1(i)
%         mkdir(['J:\NC_T1_nii\T1\',pathdir1(i).name])
      files=dir([pathdir1(i).folder,'\',pathdir1(i).name,'\smwrc2*']);
    copyfile([pathdir1(i).folder,'\',pathdir1(i).name,'\',files.name],['J:\MCI_T1_WM_Smooth\',pathdir1(i).name,'.nii'])
    end
    
end






pathdir1=dir('J:\NC_T1_nii\T1');
pathdir1(1:2)=[];
for i=1:length(a2)
    if a2(i)
%         mkdir(['J:\NC_T1_nii\T1\',pathdir1(i).name])
      files=dir([pathdir1(i).folder,'\',pathdir1(i).name,'\smwrc2*']);
    copyfile([pathdir1(i).folder,'\',pathdir1(i).name,'\',files.name],['J:\NC_T1_WM_Smooth\',pathdir1(i).name,'.nii'])
    end
    
end

pathdir1=dir('J:\ADNI_NC_gt75\NC_T1_nii\T1');
pathdir1(1:2)=[];
for i=1:length(a3)
    if a3(i)
%         
%   mkdir(['J:\ADNI_NC_gt75\NC_T1_nii\T1\',pathdir1(i).name])
      files=dir([pathdir1(i).folder,'\',pathdir1(i).name,'\smwrc2*']);
    copyfile([pathdir1(i).folder,'\',pathdir1(i).name,'\',files.name],['J:\NC_T1_WM_Smooth\',pathdir1(i).name,'.nii'])
    end
    
end


MCI_NC_COV_new=MCI_NC_COV([b1;b2],:);









%%
T1=dir('J:\MCI_T1_nii\T1');
fmri=dir('J:\MCI_fmri_nii\preprocess\Normalise');
[a,b]=ismember({fmri.name},{T1.name});


TIV_mci={};
TIV_nc={};
num=0;
for i=1:479
    num=num+1;
    temp=TIV_all{i,1};
    temps=strsplit(temp,'\');
    name(num,1)=temps(4);
    
    
end
for i=480:740
    num=num+1;
    temp=TIV_all{i,1};
    temps=strsplit(temp,'\');
    name(num,1)=temps(5);
    
    
end

files1=dir('I:\ADNI\MCI\Spreadingnew\*.mat');
files2=dir('I:\ADNI\NC_allnew\*.mat');

num=0;
for i=1:length(files1)

    num=num+1;
    
    name_fmri(num,1)=cellstr(files1(i).name(1:end-21));

end
for j=1:length(files2)

   num=num+1;
  name_fmri(num,1)=cellstr(files2(j).name(1:end-21));

end


[a,b]=ismember(name_fmri,name);
b(b==0)=[];

TIV=TIV_all(b,:);



for i=1:length(name_fmri)
    i
    path1=['J:\MCI_T1\',name_fmri{i}];
    path2=['J:\NC_T1\',name_fmri{i}];
    path3=['J:\ADNI_NC_gt75\NC_T1\',name_fmri{i}];
    if exist(path1)
        cd(path1)
        
        
    elseif exist(path2)
        cd(path2)
    else
         cd(path3)
    end
    files=dir('*\*');
    info=dicominfo([files(3).folder,'\',files(3).name]);    
    
%     TIV{i,5}=info.Manufacturer;
    if strcmp(info.Manufacturer,'SIEMENS')
        TIV{i,5}=1;
    else
        TIV{i,5}=2;
    end
    
    
end
[Data_MCI, VoxelSize, FileList, Header] = y_ReadAll('J:\MCI_T1_smooth');
Data=cat(4,Data_MCI,Data_NC);
graymaskpath='D:\DPABI_V6.0_210501\Templates\GreyMask_02_91x109x91.img';
graymask=spm_read_vols(spm_vol(graymaskpath));
data=[];
for i=1:size(Data,4)
temp=Data(:,:,:,i);
data(:,i)=temp(graymask~=0);
end

clear b r SSE SSR T TF_ForContrast  Cohen_f2 TF_ForContrast_brain1 TF_ForContrast_brain2
for i=1:size(data,1)
DependentVariable=data(i,:)';
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,[Regressors_all],Contrast,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
TF_ForContrast_brain1(i)=TF_ForContrast;
end

Regressors_allnewTemp=MCI_NC_COV_new(:,2:7);
Regressors_allnewTempnew=Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,4,6]);
for i=1:size(data,1)
DependentVariable=data(i,Regressors_allnewTemp(:,5)==1)';
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,Regressors_allnewTempnew,[1,0,0,0,0],'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
TF_ForContrast_brain1(i)=TF_ForContrast;
end


vimage=spm_vol(graymaskpath);
vimage.descrip=sprintf('DPABI{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',df,0,0,0,0);
vimage.dt=[16,0];
result=zeros(91,109,91);
result(graymask~=0)=TF_ForContrast_brain1;
vimage.fname='J:\T2_GM.nii';


%%



files1=dir('I:\ADNI\MCI\Spreading_Smooth_fc_fdr\*.mat');
files2=dir('I:\ADNI\NC_all_Spreading_Smooth_fc_fdr\*.mat');
S_all=[];
num=0;
for i=1:length(files1)
    num
    num=num+1;
    load([files1(i).folder,'\',files1(i).name]);
    for j=1:size(S,3);
       S_all(num,:,:,j)=S(:,:,j);
    end
end
for j=1:length(files2)
    num
   num=num+1;
   load([files2(j).folder,'\',files2(j).name]);
   for k=1:size(S,3);
    S_all(num,:,:,k)=S(:,:,k);
   end
end
S_mci_nc=S_all(:,:,:,7);

for i=1:size(S_mci_nc,1)
    
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==0)=max(temp(:));
    S_mci_nc(i,:,:)=temp;
end


mask=(ones(416,416))-eye(416);
for i=1:size(S_mci_nc,1)
temp=squeeze(S_mci_nc(i,:,:));
Y(:,i)=temp(mask~=0);
end
Y=combat(Y,Regressors_allnewTemp(:,5),[],1);
for i=1:size(S_mci_nc,1)
temp=zeros(416,416);
temp(mask~=0)=Y(:,i);
S_mci_ncTemp(i,:,:)=temp;
end

for i=1:size(S_all,1)
    for j=1:size(S_all,4)
        temp=S_all(i,:,:,j);
        temp=squeeze(temp);
        temp(temp==0)=max(temp(:));
        temp=temp-diag(diag(temp));
       S_all(i,:,:,j)=temp;
            
    end
end




[a,b]=ind2sub([490,416],find(times_in_Temp==0));
unique(a)

[a1,b1]=ind2sub([490,416],find(times_in_Temp==0));

 


S_mci_ncTemp=S_mci_nc;
S_mci_ncTemp(unique(a),:,:)=[];
avgGroupZ = squeeze(nanmean(S_mci_ncTemp(Regressors_allnewTemp(:,1)==1,:,:)));
avgGroupZ1 = squeeze(nanmean(S_mci_ncTemp(Regressors_allnewTemp(:,1)==-1,:,:)));
% for i=1:416
%     for j=1:size(S_all,4)
%     DependentVariable= times_out(:,i,j);
%     [b,r,SSE,SSR_1, T, TF_ForContrast_time_out(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regressors_allnew(:,:),Contrastnew,'T'); 
%     DependentVariable= times_in(:,i,j);
%     [b,r,SSE,SSR_1, T, TF_ForContrast_time_in(i,j), Cohen_f2] = y_regress_ss(DependentVariable,Regressors_allnew(:,:),Contrastnew,'T');
%     end
%     
% end
Regressors_allnewTemp(unique(a),:)=[];


Contrastnew=zeros(1,size(Regressors_allnewTemp,2));
Contrastnew(1)=1;
for j=1:416
    for k=1:416

    if j~=k
        DependentVariable=S_mci_ncTemp(:,j,k);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp,Contrastnew,'T'); 
         
    else
 
         TF_ForContrast(j,k)=0;

    end
end

end
T=TF_ForContrast;
T(isnan(T))=0;

mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regressors_allnewTemp,1)-size(Regressors_allnewTemp,2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
%     P   = 0.00001;
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,416);
Thresholdedmask(mask~=0)=Thresholded;

T_fdr=T.*Thresholdedmask;   




files1=dir('I:\ADNI\MCI\Spreading_Smooth_new\*.mat');
files2=dir('I:\ADNI\NC_all_Spreading_Smooth_new\*.mat');
T_all=[];
num=0;
for i=1:length(files1)
    num=num+1;
    load([files1(i).folder,'\',files1(i).name]);
    for j=1:size(T,3);
       T_all(num,:,:,j)=T(:,:,j);
    end
end
for j=1:length(files2)
    num=num+1;
   load([files2(j).folder,'\',files2(j).name]);
   for j=1:size(T,3);
    T_all(num,:,:,j)=T(:,:,j);
   end
end

files1=dir('I:\ADNI\MCI\Spreading_Smooth_new\*.mat');
files2=dir('I:\ADNI\NC_all_Spreading_Smooth_new\*.mat');
S_all=[];
num=0;
for i=1:length(files1)
    num=num+1;
    load([files1(i).folder,'\',files1(i).name]);
    for j=1:size(S,3);
       S_all(num,:,:,j)=S(:,:,j);
    end
end
for j=1:length(files2)
    num=num+1;
   load([files2(j).folder,'\',files2(j).name]);
   for j=1:size(S,3);
    S_all(num,:,:,j)=S(:,:,j);
   end
end


T_mci_ncTemp=T_all(:,:,:,4);
T_mci_ncTemp(unique(a),:,:)=[];

T_mciTemp=T_mci_ncTemp;  
times=unique(T_mciTemp(:));
for i=1:416
    i
    for j=1:416
        temp=T_mciTemp(:,i,j);
        for k=1:length(times)
            count = sum(temp == times(k));
            hypo=1/length(times);
            if count==0
                p_value(k)=1;
            else
                
            p_value(k) = 1 - binocdf(count - 1, length(temp), hypo);
            end
               
        end
        [min_v(i,j),min_ind(i,j)]=min(p_value);
        
        
    end
end






for i=1:size(S_all,4)
    for j=1:size(S_all,1)
        temp=squeeze(S_all(j,:,:,i));
        mean_time(j,i)=sum(sum(temp')./416);
        
        
    end
end

for j=1:size(S_all,4)
 DependentVariable=mean_time(:,j);
 [b,r,SSE,SSR_1, T, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,Regressors_allnewTemp,Contrastnew,'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end

for j=1:size(S_all,4)
 DependentVariable=mean_time(Regressors_allnewTemp(:,5)==1,j);
 mean_timenew(:,j)= zscore(DependentVariable - (Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[2,3,4,6,7]))*(((Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[2,3,4,6,7]))\DependentVariable)));
 end

for j=1:size(S_all,4)
 DependentVariable=mean_time(Regressors_allnewTemp(:,5)==1,j);
 [b,r,SSE,SSR_1, T, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,...
     Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 ,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end

for j=1:size(S_all,4)
 DependentVariable=mean_time(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.2,j);
 [b,r,SSE,SSR_1, T, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,...
     Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.2 ,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end


for j=1:size(S_all,4)
    effect_size(j)=computeCohen_d(mci_mean(:,j),nc_mean(:,j),'independent');
end



Regress_Simens=Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]);
for i=1:9
    [p(i), h, stats] = ranksum(mean_timenew(Regress_Simens(:,1)==1 & Regress_Simens(:,5)<0.5,i),...
        mean_timenew(Regress_Simens(:,1)==-1 & Regress_Simens(:,5)<0.5,i));
    t_meantime(i)=stats.zval;
end
for i=1:size(S_all,4)
    effect_size(i)=computeCohen_d(mean_timenew(Regress_Simens(:,1)==1 & Regress_Simens(:,5)<0.2,i),...
        mean_timenew(Regress_Simens(:,1)==-1 & Regress_Simens(:,5)<0.2,i),'independent');
end

for i=1:9
   temp=(mean_time(Regressors_allnewTemp(:,1)==1 & Regressors_allnewTemp(:,6)<0.5 & Regressors_allnewTemp(:,5)==1,i)-...
       mean(mean_time(Regressors_allnewTemp(:,1)==-1 & Regressors_allnewTemp(:,6)<0.5 & Regressors_allnewTemp(:,5)==1,i)))...
       ./std(mean_time(Regressors_allnewTemp(:,1)==-1 & Regressors_allnewTemp(:,6)<0.5 & Regressors_allnewTemp(:,5)==1,i));
   
   mean(temp)
   hold on;
   scatter(ones(size(temp))*i,temp);
   scatter(i,mean(temp));
   
end

Regress_Simens=Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]);
cInd=find(Regress_Simens(:,1)==-1);
pInd=find(Regress_Simens(:,1)==1);
for i=1:9
    i
   data=reshape(squeeze(S_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:,i)),...
       size(squeeze(S_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:,i)),1),[]);
   
%    for j=1:size(data,1)
%        data(j,:)=(data(j,:)-min(data(j,:)))./(max(data(j,:))-min(data(j,:)));
%    end
   
%    z = myzscore(data,0,2);   
%    m_z = mean(z,2,'omitnan');
   
%    zc = zScoreToSubset(data,cInd);
%    m_zc = mean(zc,2,'omitnan');
%    
    z = myzscore(data,0,2);
   zz = zScoreToSubset(z,cInd);
   m_zz = mean(zz,2,'omitnan');
     zvalue(i)=mean(m_zz(1:149),1);
     
     
     zvalue_NC(i)=mean(m_zz(150:end),1);
%    zc = zScoreToSubset(data,cInd); % z-scoring to controls the individually z-scored data
%     m_zc = mean(zc,2,'omitnan');
%     m_zc= m_zc - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zc);
%     m_zz= m_zz - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zz);
%      [h,p,c,d]=ttest2(m_zc(pInd),m_zc(cInd),'Alpha',0.05);
%     [b,r,SSE,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(m_zc,...
%     Regress_Simens,[1,0,0,0,0,0],'T');
    
%     [p(i), h, stats] = ranksum(m_zz(pInd),m_zz(cInd));
%     t_meantime(i)=stats.zval;
    [b,r,SSE,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(mean(data,2),...
    Regress_Simens,[1,0,0,0,0,0],'T');
   
end
T=TF_ForContrast';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regress_Simens,1)...
    -size(Regress_Simens,2);
PMap=2*(1-tcdf(abs(Ttemp), df));




cInd=find(Regress_Simens(:,1)==-1);
pInd=find(Regress_Simens(:,1)==1);
for i=1:9
    i
   data=reshape(squeeze(T_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:,i)),...
       size(squeeze(T_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:,i)),1),[]);
   
%    for j=1:size(data,1)
%        data(j,:)=(data(j,:)-min(data(j,:)))./(max(data(j,:))-min(data(j,:)));
%    end
% data(data==22)=NaN;
   
%    z = myzscore(data,0,2);   
%    m_z = mean(data,2,'omitnan');
   
%    zc = zScoreToSubset(data,cInd);
%    m_zc = mean(zc,2,'omitnan');
%    
    z = myzscore(data,0,2);
   zz = zScoreToSubset(z,cInd);
   m_zz = mean(zz,2,'omitnan');
     zvalue(i)=mean(m_zz(1:149),1);
     zvalueNC(i)=mean(m_zz(150:end),1);
          mci_mean(:,i)=m_zz(1:149);
      nc_mean(:,i)=m_zz(150:end);
%    zc = zScoreToSubset(data,cInd); % z-scoring to controls the individually z-scored data
%     m_zc = mean(zc,2,'omitnan');
%     m_zc= m_zc - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zc);
%     m_zz= m_zz - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zz);
%      [h,p,c,d]=ttest2(m_zc(pInd),m_zc(cInd),'Alpha',0.05);
    [b,r,SSE,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(m_zz,...
    Regress_Simens,[1,0,0,0,0,0],'T');
    
%     [p(i), h, stats] = ranksum(m_zz(pInd),m_zz(cInd));
%     t_meantime(i)=stats.zval;
   
end




for i=1:9
    i
   data=reshape(squeeze(T_all( Regressors_allnewTemp(:,6)<0.5,:,:,i)),...
       size(squeeze(T_all(Regressors_allnewTemp(:,6)<0.5,:,:,i)),1),[]);
   
%    for j=1:size(data,1)
%        data(j,:)=(data(j,:)-min(data(j,:)))./(max(data(j,:))-min(data(j,:)));
%    end
% data(data==22)=NaN;
   
%    z = myzscore(data,0,2);   
%    m_z = mean(data,2,'omitnan');
   
%    zc = zScoreToSubset(data,cInd);
%    m_zc = mean(zc,2,'omitnan');
%    
    z = myzscore(data,0,2);
   zz = zScoreToSubset(z,cInd);
   m_zz = mean(zz,2,'omitnan');
     zvalue(i)=mean(m_zz(1:149),1);

%    zc = zScoreToSubset(data,cInd); % z-scoring to controls the individually z-scored data
%     m_zc = mean(zc,2,'omitnan');
%     m_zc= m_zc - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zc);
%     m_zz= m_zz - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zz);
%      [h,p,c,d]=ttest2(m_zc(pInd),m_zc(cInd),'Alpha',0.05);
% %     [b,r,SSE,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(m_zz,...
% %     Regress_Simens,[1,0,0,0,0,0],'T');
    
%     [p(i), h, stats] = ranksum(m_zz(pInd),m_zz(cInd));
%     t_meantime(i)=stats.zval;
   
end


cInd=find(Regressors_allnewTemp(:,1)==-1 & Regressors_allnewTemp(:,6)<0.5);
pInd=find(Regressors_allnewTemp(:,1)==1& Regressors_allnewTemp(:,6)<0.5);
for i=1:9
    i
   data=reshape(squeeze(T_all(:,:,:,i)),...
       size(squeeze(T_all(:,:,:,i)),1),[]);
   
%    for j=1:size(data,1)
%        data(j,:)=(data(j,:)-min(data(j,:)))./(max(data(j,:))-min(data(j,:)));
%    end
% data(data==22)=NaN;
   
%    z = myzscore(data,0,2);   
%    m_z = mean(data,2,'omitnan');
   
%    zc = zScoreToSubset(data,cInd);
%    m_zc = mean(zc,2,'omitnan');
%    
    z = myzscore(data,0,2);
   zz = zScoreToSubset(z,cInd);
   m_zz = mean(zz,2,'omitnan');
     zvalue(i)=mean(m_zz(Regressors_allnewTemp(:,1)==1 & Regressors_allnewTemp(:,6)<0.5),1);
%    zc = zScoreToSubset(data,cInd); % z-scoring to controls the individually z-scored data
%     m_zc = mean(zc,2,'omitnan');
%     m_zc= m_zc - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zc);
%     m_zz= m_zz - Regress_Simens(:,2:end)*((Regress_Simens(:,2:end))\m_zz);
%      [h,p,c,d]=ttest2(m_zc(pInd),m_zc(cInd),'Alpha',0.05);
%     [b,r,SSE,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(m_zz,...
%     Regress_Simens,[1,0,0,0,0,0],'T');
    
%     [p(i), h, stats] = ranksum(m_zz(pInd),m_zz(cInd));
%     t_meantime(i)=stats.zval;
   
end


for j=1:9
    mci_mean(:,j)=mean_time((Regressors_allnewTemp(:,5)==1) & (Regressors_allnewTemp(:,1)==1) ,j);
end
for j=1:9
    nc_mean(:,j)=mean_time((Regressors_allnewTemp(:,5)==1) & (Regressors_allnewTemp(:,1)==-1) ,j);
end



for j=1:9
    mci_mean(:,j)=mean_timenew(1:163,j);
end
for j=1:9
    nc_mean(:,j)=mean_timenew(164:end,j);
end



for i=1:size(mci_mean,2)
             
    sss=mci_mean(:,i);
%     mci_Qmed(i)=median(sss);
    mci_Qmed(i)=mean(sss);
    mci_Q25e(i)=mci_Qmed(i)-quantile(sss,.25);
    mci_Q75e(i)=quantile(sss,.75)-mci_Qmed(i);
end

for i=1:size(nc_mean,2)
             
    sss=nc_mean(:,i);
%     nc_Qmed(i)=median(sss);
     nc_Qmed(i)=mean(sss);
    nc_Q25e(i)=nc_Qmed(i)-quantile(sss,.25);
    nc_Q75e(i)=quantile(sss,.75)-nc_Qmed(i);
end
hold on;
hl=errorbar(0.9:size(mci_mean,2)-1+0.9,mci_Qmed,mci_Q25e,mci_Q75e,'LineStyle', 'None');
hr=errorbar(1.1:size(nc_mean,2)-1+1.1,nc_Qmed,nc_Q25e,nc_Q75e,'LineStyle', 'None');
hl.Marker='o';
hl.MarkerEdgeColor= [.2 .2 .2];
hl.MarkerFaceColor= [208,79,88]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hl.MarkerSize=10;  
hl.LineWidth=1;
hl.Color='k';

hr.Marker='o';
hr.MarkerEdgeColor= [.2 .2 .2];
hr.MarkerFaceColor= [23,112,159]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hr.MarkerSize=10;  
hr.LineWidth=1;
hr.Color='k';

set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold')
set(gca,'LineWidth',2);
set(gca,'box','off');
set(gcf,'color','w')
hold on
plot(1.1:size(nc_mean,2)-1+1.1,nc_Qmed,'-','LineWidth',2,'Color',[23,112,159]./255);
plot(0.9:size(mci_mean,2)-1+0.9,mci_Qmed,'-','LineWidth',2,'Color',[208,79,88]./255);

gretna_plot_violin({mci_mean,nc_mean},{'MCI','NC'},{'Reg1','Reg2','Reg3','Reg4','Reg5','Reg6','Reg7','Reg8','Reg9'},'dot');

for i=1:9
    startIndex=0.9;
    endIndex=1.1;
    scatter((startIndex+(i-1))*ones(size(mci_mean,1),1),mci_mean(:,i),'filled','MarkerFaceColor',[213,135,73]./255)

    hold on;
    scatter((endIndex+(i-1))*ones(size(nc_mean,1),1),nc_mean(:,i),'filled','MarkerFaceColor',[30,152,199]./255)
    
    
end
 plot([1:9],mean(mci_mean),'Marker','^','color',[213,135,73]./255,'MarkerSize',10,'Linewidth',2,'MarkerFaceColor',[1,1,1]);
 plot([1:9],mean(nc_mean),'Marker','o','color',[30,152,199]./255,'MarkerSize',10,'Linewidth',2,'MarkerFaceColor',[1,1,1]);

scatter(0.98*ones(size(mci_mean,1),1),mci_mean(:,1))
hold on;
scatter(1.02*ones(size(nc_mean,1),1),nc_mean(:,1))

scatter(1.98*ones(size(mci_mean,1),1),mci_mean(:,2))

scatter(2.02*ones(size(nc_mean,1),1),nc_mean(:,2))







tempdir=dir('I:\ADNI\MCI\Spreading_Smooth');
tempdir(1:2)=[];
times_MCI=[];
times_MCI(cellfun(@(folder) contains(folder,'T1'),{tempdir.name}))=1;
times_MCI(cellfun(@(folder) contains(folder,'T2'),{tempdir.name}))=2;
times_MCI(cellfun(@(folder) contains(folder,'T3'),{tempdir.name}))=3;
times_MCI(cellfun(@(folder) contains(folder,'T4'),{tempdir.name}))=4;
times_MCI(cellfun(@(folder) contains(folder,'T5'),{tempdir.name}))=5;
times_MCI(cellfun(@(folder) contains(folder,'T6'),{tempdir.name}))=6;

tempdir=dir('I:\ADNI\NC_all_Spreading_Smooth');
tempdir(1:2)=[];
times_NC=[];
times_NC(cellfun(@(folder) contains(folder,'T1'),{tempdir.name}))=1;
times_NC(cellfun(@(folder) contains(folder,'T2'),{tempdir.name}))=2;
times_NC(cellfun(@(folder) contains(folder,'T3'),{tempdir.name}))=3;
times_NC(cellfun(@(folder) contains(folder,'T4'),{tempdir.name}))=4;

times=[times_MCI';times_NC'];
times(208)=[];

times_MCI(208)=[];


tempdir=dir('I:\ADNI\MCI\Spreading_Smooth');
tempdir(1:2)=[];
name_MCI=cellfun(@(folder) folder(1:10),{tempdir.name},'UniformOutput',false);
tempdir=dir('I:\ADNI\NC_all_Spreading_Smooth');
tempdir(1:2)=[];
name_NC=cellfun(@(folder) folder(1:10),{tempdir.name},'UniformOutput',false);
nameall=[name_MCI,name_NC];
nameall(208)=[];
names=nameall(times(1:248)==2);


tempdir=dir('I:\ADNI\MCI\Spreading_Smooth');
tempdir(1:2)=[];
name_MCI=cellfun(@(folder) folder(1:13),{tempdir.name},'UniformOutput',false);
tempdir=dir('I:\ADNI\NC_all_Spreading_Smooth');
tempdir(1:2)=[];
name_NC=cellfun(@(folder) folder(1:13),{tempdir.name},'UniformOutput',false);







Regressors_allnewTemp=[Regressors_allnewTemp,times];

S_mci_nc=[];
for sub=1:size(S_all,1)
    sub
    temp=squeeze(S_all(sub,:,:,7));
    mu = nanmean(temp(:));
    sigma = nanstd(temp(:),0);
    z = bsxfun(@minus,temp(:), mu);
    z = bsxfun(@rdivide, z, sigma);
    S_mci_nc(sub,:,:)=reshape(z,416,416);
end
S_mci_nc=squeeze(S_all(:,:,:,7));
S_mci_nc(208,:,:)=[];
for i=1:size(S_mci_nc,1)
    times_out_Temp(i,:)=sum(squeeze(S_mci_nc(i,:,:)))./416;
    times_in_Temp(i,:)=sum(squeeze(S_mci_nc(i,:,:)),2)./416;
end


load('I:\ADNI\T_all.mat')
S_mci_nc=squeeze(T_all(:,:,:,4));

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end


for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
%     temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end

for i=1:size(S_mci_nc_temp,1)
    times_out_Temp(i,:)=sum(squeeze(S_mci_nc_temp(i,:,:)))./416;
    times_in_Temp(i,:)=sum(squeeze(S_mci_nc_temp(i,:,:)),2)./416;
end



times_in_Temp(times_in_Temp==0)=NaN;
times_out_Temp(times_out_Temp==0)=NaN;

data=times_out_Temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:);
z = myzscore(data,0,2);
zz = zScoreToSubset(z,cInd);

for i=1:416
    DependentVariable=times_in_Temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regress_Simens,Contrastnew(1:6),'T'); 
    
end



for i=1:416
    DependentVariable=zz(:,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regress_Simens,Contrastnew(1:6),'T'); 
    
end
T=TF_ForContrast_time_outTemp';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regress_Simens,1)...
    -size(Regress_Simens,2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T'.*Thresholdedmask;



[~, idx] = sort(TF_ForContrast_time_outTemp');
Temp_out_Grank = sort_back([1:size(TF_ForContrast_time_outTemp,2)]', idx);


[~, idx] = sort(TF_ForContrast_time_outTemp');
Temp_in_Grank = sort_back([1:size(TF_ForContrast_time_outTemp,2)]', idx);




for i=1:416
    DependentVariable=times_in_Temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regressors_allnewTemp...
         (Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    
end

T=TF_ForContrast_time_outTemp';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)...
    -size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T'.*Thresholdedmask;





for i=1:18
    DependentVariable=sum(times_in_Temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5 ,yeoindex17==i)')';
[b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    
     
     PMap(i)=2*(1-tcdf(abs(TF_ForContrast_time_outTemp(i)), df));
    
    
end

for i=1:8
    DependentVariable=sum(times_out_Temp(Regressors_allnewTemp(:,5)==1,yeoindex==i)')';
[b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp, Cohen_f2] =...
         y_regress_ss(DependentVariable,Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,6,7]),Contrastnew(1:5),'T'); 
PMap(i)=2*(1-tcdf(abs(TF_ForContrast_time_outTemp), df));
    
    
end


for i=1:18
    DependentVariable=sum(zz(: ,yeoindex17==i)')';
[b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regress_Simens,Contrastnew(1:6),'T'); 

     PMap(i)=2*(1-tcdf(abs(TF_ForContrast_time_outTemp(i)), df));
    
    
end

for i=1:8
    DependentVariable=sum(times_out_Temp(Regressors_allnewTemp(:,5)==1,yeoindex==i)')';
[b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp, Cohen_f2] =...
         y_regress_ss(DependentVariable,Regress_Simens,Contrastnew(1:6),'T'); 
PMap(i)=2*(1-tcdf(abs(TF_ForContrast_time_outTemp), df));
    
    
end



for i=1:size(S_all,4)
    S_mci_nc=squeeze(S_all(:,:,:,i));
%     for sub=1:size(S_mci_nc,1)
%         temp=squeeze(S_mci_nc(sub,:,:));
%         
%         maxV=max(temp(:));
%         minV=min(temp(:));
%         S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
%     end
        
        
        data=S_mci_nc_temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:);
        data=reshape(data,size(data,1),[]);
        data = myzscore(data,0,2);
        data = zScoreToSubset(data,cInd);
        data=reshape(data,size(data,1),416,416);
        
    for j=1:416
        for k=1:416

        if j~=k
            DependentVariable=data(:,j,k);
            [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
               Regress_Simens,Contrastnew(1:6),'T'); 

        else

             TF_ForContrast(j,k)=0;

        end
    end

    end
    T=TF_ForContrast;
    T(isnan(T))=0;

    mask=T;
    mask(mask~=0)=1;

    Ttemp=T(mask~=0);
    df=size(Regress_Simens,1)-...
        size(Regress_Simens,2);
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    %     P   = 0.00001;
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
        Thresholded(find(PMap<=P))=1;
    end
    Thresholdedmask=zeros(416,416);
    Thresholdedmask(mask~=0)=Thresholded;
    T_fdr=T.*Thresholdedmask;   
    length_edges(i)=length( find(T_fdr~=0));

      
end

clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for step =1:9
S_mci_nc=squeeze(T_all(:,:,:,step));
S_mci_nc=S_mci_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);


for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
%     temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end



% for sub=1:size(S_mci_nc,1)
%     temp=squeeze(S_mci_nc(sub,:,:));
% 

%     maxV=max(temp(:));
%     minV=min(temp(:));
%     S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
% end

new=S_mci_nc_temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:);
new(isnan(new))=0;

for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data = myzscore(data,0,2);
data = zScoreToSubset(data,cInd);
data=reshape(data,size(data,1),18,18);

% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,18,18);
% end

for i=1:18
for j=1:18
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regress_Simens,Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regress_Simens,1)-...
    size(Regress_Simens,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;

T_all_fdr(:,:,step)=T_fdr;
end

new=S_mci_nc_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
% data = myzscore(data,0,2);
% data = zScoreToSubset(data,cInd);
% data=reshape(data,size(data,1),18,18);

% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,18,18);
% end
data_harmonized = combat(data', Regressors_allnewTemp(:,5)', Regressors_allnewTemp(:,[1,2,3,4,6,7]), 1);
data_harmonized = combat(data', Regressors_allnewTemp(:,5)',[], 1);
data=reshape(data_harmonized',size(data_harmonized',1),18,18);
for i=1:18
for j=1:18
        DependentVariable=data(Regressors_allnewTemp(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),1)-...
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;












new=S_mci_nc_temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,...
    yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(:,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regress_Simens,Contrastnew(1:6),'T'); 
    end
    T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size(Regress_Simens,1)-...
        size(Regress_Simens,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end


T2map=[];
for net=1:18
    network=new(:,:,yeoindex17==net);

    network=mean(network,3);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(:,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regress_Simens,Contrastnew(1:6),'T'); 
    end
    T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size(Regress_Simens,1)-...
        size(Regress_Simens,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end


for i=1:18
    
    try
    plot_hemispherewithboundarynew_zhengnew(T_fdr(yeo_17labels_17_to_7,i),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
    print(gcf,['I:\ADNI\yeonet_',num2str(i),'input.png'],'-dpng','-r600');
    close all;
    
    catch
    end
end

new=S_mci_nc_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp(:,5)',Regressors_allnewTemp(:,[2,3,4,6]), 0);
    for i=1:416
        DependentVariable=data_harmonized(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
for i=1:18
    
    try
    plot_hemispherewithboundarynew_zhengnew(T_fdr(yeo_17labels_17_to_7,i),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
%     print(gcf,['I:\ADNI\Combatyeonet_',num2str(i),'.png'],'-dpng','-r600');
%     close all;
    
    catch
    end
end

real_data=T2map(:,yeo_17labels_17_to_7);
real_data=T_fdr(yeo_17labels_17_to_7,:)';
load('H:\light_therapysyy\baselinall\controal\asl_enriched\spin.mat')
[U,S,V]=svd(T2map(:,yeo_17labels_17_to_7));
[U,S,V]=svd(real_data);
for i = 1:5000
    % ? 打乱 T 值矩阵：示例是对每行（网络）打乱列顺序（或更复杂的打乱）
    T_perm = real_data(:,perm_id(:,i));
    
%     T_perm = real_data(:,randperm(416));

    [U_perm, ~, V_perm] = svd(T_perm);
    v1_null(i,:) = V_perm(:,1);  % 存储置换主成分 loading
    u1_null(i,:) = U_perm(:,1);  % 存储置换主成分 loading
end
for i=1:416
 up=prctile(v1_null(:,i),97.5);
 down=prctile(v1_null(:,i),2.5);
real_r=V(i,1);
if real_r<=down
    p(i)=sum(v1_null(:,i)<=real_r)/5000;

elseif real_r>=up
      p(i)=sum(v1_null(:,i)>=real_r)/5000;

else
    p(i)=0.06;
end

end
for i=1:416
 up=prctile(u1_null(:,i),97.5);
 down=prctile(u1_null(:,i),2.5);
real_r=U(i,1);
if real_r<=down
    p(i)=sum(u1_null(:,i)<=real_r)/5000;

elseif real_r>=up
      p(i)=sum( u1_null(:,i)>=real_r)/5000;

else
    p(i)=0.06;
end

end


for i=1:416
    p_values(i) = mean(abs(v1_null(:,i))>= abs(V(i,1)));
end
p_values = mean(abs(v1_null') >= abs(real_r), 2);




for i = 1:416
    if V(i,1)> median(v1_null(:,i))
        p_perm_low(i) = sum(v1_null(:,i) > V(i,1))/5000;
    elseif V(i,1)< median(v1_null(:,i))
         p_perm_low(i) = sum(v1_null(:,i)< V(i,1))/5000;
    else
        
        p_perm_low(i) = 0.06;    
    end
end



FDR_high=pval_adjust(p_perm_low,'fdr');






clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for step =1:9
S_mci_nc=squeeze(S_all(:,:,:,step));
% for sub=1:size(S_mci_nc,1)
%     temp=squeeze(S_mci_nc(sub,:,:));
% 
%     maxV=max(temp(:));
%     minV=min(temp(:));
%     S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
% end
for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end
new=S_mci_nc_temp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:);
new(isnan(new))=0;

for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:8];
    mask=yeoindex;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,8,8);
% end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
% data = myzscore(data,0,2);
% data = zScoreToSubset(data,cInd);
data=reshape(data,size(data,1),8,8);

for i=1:8
for j=1:8
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regress_Simens,Contrastnew(1:6),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regress_Simens,1)-...
    size(Regress_Simens,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(8,8);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;
T_all(:,:,step)=T_fdr;
end




load('J:\MMSE_MCI.mat')
MMSE_MCI=table2array(MMSE);
load('J:\MMSE_NC.mat')
MMSE_NC=table2array(MMSE);
MMSE=[MMSE_MCI;MMSE_NC];
MMSE(208)=[];

MMSE=MMSE(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5);
[r,p]=corr(data(~isnan(MMSE),17,6),MMSE(~isnan(MMSE)))


data_mat_mean_all_temp=data_mat_mean_all;
data_mat_mean_all_temp(isnan(MMSE),:,:)=[];
Regress_Simens_temp=Regress_Simens;
Regress_Simens_temp(isnan(MMSE),:)=[];

mmse=MMSE;
mmse(isnan(mmse))=[];
scatter(mmse(Regress_Simens_temp(:,1)==1),data_mat_mean_all_temp(Regress_Simens_temp(:,1)==1,18,2));
hold on;
scatter(mmse(Regress_Simens_temp(:,1)==-1),data_mat_mean_all_temp(Regress_Simens_temp(:,1)==-1,18,2));





[x,y]=ind2sub([size(T_fdr,1),size(T_fdr,2)],find(T_fdr));
r=zeros(size(T_fdr,1),size(T_fdr,2));
p=zeros(size(T_fdr,1),size(T_fdr,2));
for i=1:length(x)
    [r(x(i),y(i)),p(x(i),y(i))]=corr(data(~isnan(MMSE),x(i),y(i)),MMSE(~isnan(MMSE)));
   
end

cortex=spm_read_vols(spm_vol('E:\WorkSpace\MRI_MASK\r2schaefer400MNI.nii'));
cortex(isnan(cortex))=0;
subcortex=spm_read_vols(spm_vol('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii'));

M_youngcortex=[];M_youngsubcortex=[];
MCIT1dir=dir('J:\MCI_T1_smooth\*.nii');
for i=1:length(MCIT1dir)
    i
    image=spm_read_vols(spm_vol([MCIT1dir(i).folder,'\',MCIT1dir(i).name]));
    for roi=1:400
        M_youngcortex(i,roi)=mean(image(cortex==roi));
    end
    
    for roi=1:16
        M_youngsubcortex(i,roi)=mean(image(subcortex==roi));
    end
    MCI_T1_data=[ M_youngcortex';M_youngsubcortex'];
end

M_youngcortex=[];M_youngsubcortex=[];
clear M_youngcortex M_youngsubcortex
MCIT1dir=dir('J:\NC_T1_smooth\*.nii');
for i=1:length(MCIT1dir)
    image=spm_read_vols(spm_vol([MCIT1dir(i).folder,'\',MCIT1dir(i).name]));
    for roi=1:400
        M_youngcortex(i,roi)=mean(image(cortex==roi));
    end
    
    for roi=1:16
        M_youngsubcortex(i,roi)=mean(image(subcortex==roi));
    end
    NC_T1_data=[ M_youngcortex';M_youngsubcortex'];
end
T1_all=[MCI_T1_data';NC_T1_data'];


for i=1:416
    DependentVariable=T1_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_T1(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regressors_allnewTemp...
         (Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    
end

network = combat(T1_all(Regressors_allnewTemp(:,6)<0.5,:)', Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,5)',[], 0);
    for i=1:416
        i
        DependentVariable=network(i,:)';
  T = table(subname, group, times_new, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Time', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group + Time + Age + Sex + Motion + Tiv +(1|Subject)');
TF_ForContrast(i)=mdl.Coefficients.tStat(2);
p_ForContrast(i)=mdl.Coefficients.pValue(2);
Cohen_d(i)=mdl.Coefficients.Estimate(2)./sqrt(mdl.MSE);

    end




T=TF_ForContrast_T1';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)...
    -size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr_T1=T'.*Thresholdedmask;


T1_twosample=TF_ForContrast_T1;
T1_twosample_fdr=T_fdr;
[T1_twosample_sort, idx] = sort(T1_twosample');
Grank = sort_back([1:size(T1_twosample,2)]', idx);
nbin=100
ratio=[];
for b=1:nbin
    left(b)=round((length(Grank) / nbin) * b);
    right(b)=round((length(Grank) / nbin) * (b - 1));
    idx = ((Grank < round((length(Grank) / nbin) * b))+ (Grank > round((length(Grank) / nbin) * (b - 1)))) == 2;
    temp=yeoindex(idx);
    for i=1:8
        ratio(b,i)=sum(temp==i)*100/length(temp);
    end
end
b=bar(ratio,'stacked','EdgeColor',[0 0 0],'LineWidth',1.5);
b=bar(ratio,'stacked');
for i=1:8
    b(i).FaceColor=col(i,:);
end



for i=1:size(T1_all,2)

Regressors_temp=[Regressors_allnewTemp(:,[2,3,4,6,7]),times_in_Temp(:,i),times_out_Temp(:,i)];
contrastnew=zeros(1,size(Regressors_temp,2));
contrastnew(6)=1;
   DependentVariable=T1_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i);
[b,r,SSE(i),SSR(i), T, TF_ForContrast_time_outTemp(i), Cohen_f2]=y_regress_ss(DependentVariable,Regressors_temp...
         (Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:),contrastnew,'T'); 

end
r2_2=1-SSE./(SSE+SSR);


avgGroupZ=mean(S_mci_nc_temp(Regressors_allnewTemp(:,1)==1 & Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,:,:),1);
distMatrix = pdist(squeeze(avgGroupZ), 'euclidean'); 


load('I:\ADNI\T_all.mat')
S_mci_nc=squeeze(T_all(:,:,:,7));
S_mci_nc(S_mci_nc==1 | S_mci_nc==22)=0;
for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end



MCI_FC=dir('I:\ADNI\MCI\BOLD_Smooth_new\*.mat');
num=0;
for i=1:length(MCI_FC)
    num=num+1;
 load([MCI_FC(i).folder,'\',MCI_FC(i).name]);
    fc=corr(theROITimeCoursesTotal);
% fc=readNPY('I:\ADNI\MCI\fc.npy');
    fc=atanh(fc);
%         fc(fc<=0)=0;
%     fc(isinf(fc))=0;
%     fc(isnan(fc))=0;
%     cutoffValues=prctile(fc(:),90);
    threshde_FC=fc;
%     threshde_FC(threshde_FC<cutoffValues)=0;
    fc_all(num,:,:)=threshde_FC;
%     new=threshde_FC.*squeeze(S_mci_nc_temp(num,:,:));
%     time_in_new=sum(new)';
%     time_out_new=sum(new,2);
%     time_in(:,num)=time_in_new;
%     time_out(:,num)=time_out_new;
     
end

NC_FC=dir('I:\ADNI\NC_All_BOLD_Smooth_new\*.mat');
for i=1:length(NC_FC)
    num=num+1;
 load([NC_FC(i).folder,'\',NC_FC(i).name]);
    fc=corr(theROITimeCoursesTotal);
%  fc=readNPY('I:\ADNI\MCI\fc.npy');

    fc=atanh(fc);
% %         fc(fc<=0)=0;
% %     fc(isinf(fc))=0;
% %     fc(isnan(fc))=0;
% %     cutoffValues=prctile(fc(:),90);
    threshde_FC=fc;
%     threshde_FC(threshde_FC<cutoffValues)=0;
    fc_all(num,:,:)=threshde_FC;
    
%     new=threshde_FC.*squeeze(S_mci_nc_temp(num,:,:));
%     time_in_new=sum(new)';
%     time_out_new=sum(new,2);
%     time_in(:,num)=time_in_new;
%     time_out(:,num)=time_out_new;
     
end

data=time_out(:,Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5)';
z = myzscore(data,0,2);
zz = zScoreToSubset(z,cInd);
for i=1:416
    DependentVariable=zz(:,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,Regress_Simens,Contrastnew(1:6),'T'); 
    
end
T=TF_ForContrast_time_outTemp';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size(Regress_Simens,1)...
    -size(Regress_Simens,2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T'.*Thresholdedmask;











% T_mci_ncTemp=T_mci_ncTemp([MCI_index;NC_index+106],:,:);
% T_mci_ncTemp(unique(a),:,:)=[];
%%
mask=triu(ones(416,416))-eye(416);
Y=[];
for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    Y(:,i)=temp(mask~=0);
    
end

Y=combat(Y,Regressors_allnewTemp(:,5),[],1);

Ynew=[];

for i=1:size(S_mci_nc,1)
   temp=zeros(416,416);
   temp(mask~=0)=Y(:,i);
   temp=temp+temp';
    
    Ynew(i,:,:)=temp;
end








files1=dir('I:\ADNI\MCI\BOLD_Smooth\*.mat');
files2=dir('I:\ADNI\NC_All_BOLD_Smooth\*.mat');
S_all=[];
num=0;
for i=1:length(files1)
    num
    num=num+1;
    load([files1(i).folder,'\',files1(i).name]);
    S_mci_nc(num,:,:)=atanh(corr(theROITimeCoursesTotal));

end
for j=1:length(files2)
    num
   num=num+1;
   load([files2(j).folder,'\',files2(j).name]);
S_mci_nc(num,:,:)=atanh(corr(theROITimeCoursesTotal));
end


Regressors_allnewTemp=cell2mat(MCI_NC_COV_new(:,2:end));
Regressors_allnewTemp=[Regressors_allnewTemp,ones(size(Regressors_allnewTemp,1),1)];


Contrastnew=zeros(1,size(Regressors_allnewTemp,2));
Contrastnew(1)=1;
for j=1:416
    for k=1:416

    if j~=k
        DependentVariable=S_mci_nc(:,j,k);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp,Contrastnew,'T'); 
         
    else
 
         TF_ForContrast(j,k)=0;

    end
end

end

for j=1:416
    for k=1:416

    if j~=k
        DependentVariable=S_mci_nc(Regressors_allnewTemp(:,5)==1,j,k);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,6,7]),Contrastnew(1:5),'T'); 
         
    else
 
         TF_ForContrast(j,k)=0;

    end
end

end

T=TF_ForContrast;
T(isnan(T))=0;

mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);

PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
%     P   = 0.00001;
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T.*Thresholdedmask;   

%%
for i=1:size(S_all,4)
    S_mci_nc=squeeze(S_all(:,:,:,i));
%     for sub=1:size(S_mci_nc,1)
%         temp=squeeze(S_mci_nc(sub,:,:));
%         
%         maxV=max(temp(:));
%         minV=min(temp(:));
%         S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
%     end
    
    
   for sub=1:size(S_mci_nc,1)
        temp=squeeze(S_mci_nc(sub,:,:));
        mu = nanmean(temp(:));
        sigma = nanstd(temp(:),0);
        z = bsxfun(@minus,temp(:), mu);
        z = bsxfun(@rdivide, z, sigma);
        S_mci_nc(sub,:,:)=reshape(z,416,416);
    end
    
    
    
    for j=1:416
        for k=1:416

        if j~=k
            DependentVariable=S_mci_nc(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,j,k);
            [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
                Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 

        else

             TF_ForContrast(j,k)=0;

        end
    end

    end
    T=TF_ForContrast;
    T(isnan(T))=0;

    mask=T;
    mask(mask~=0)=1;

    Ttemp=T(mask~=0);
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1  & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2);
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    %     P   = 0.00001;
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
        Thresholded(find(PMap<=P))=1;
    end
    Thresholdedmask=zeros(416,416);
    Thresholdedmask(mask~=0)=Thresholded;
    T_fdr=T.*Thresholdedmask;   
    length_edges(i)=length( find(T_fdr~=0));

      
end

clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean

for step =1:9
S_mci_nc=squeeze(S_all(:,:,:,step));
S_mci_nc=S_mci_nc(:,yeo_17labels,yeo_17labels);

% for sub=1:size(S_mci_nc,1)
%     temp=squeeze(S_mci_nc(sub,:,:));
% 
%     maxV=max(temp(:));
%     minV=min(temp(:));
%     S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
% end

new=S_mci_nc;
new(isnan(new))=0;

for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
for sub=1:size(data_mat_mean_all,1)
    temp=squeeze(data_mat_mean_all(sub,:,:));
    mu = nanmean(temp(:));
    sigma = nanstd(temp(:),0);
    z = bsxfun(@minus,temp(:), mu);
    z = bsxfun(@rdivide, z, sigma);
    data_mat_mean_all(sub,:,:)=reshape(z,18,18);
end

for i=1:18
for j=1:18
DependentVariable=data_mat_mean_all(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;

T_all(:,:,step)=T_fdr;
end

clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for step =1:9
S_mci_nc=squeeze(S_all(:,:,:,step));
% for sub=1:size(S_mci_nc,1)
%     temp=squeeze(S_mci_nc(sub,:,:));
% 
%     maxV=max(temp(:));
%     minV=min(temp(:));
%     S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
% end

new=S_mci_nc;
new(isnan(new))=0;   

for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:8];
    mask=yeoindex;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
for sub=1:size(data_mat_mean_all,1)
    temp=squeeze(data_mat_mean_all(sub,:,:));
    mu = nanmean(temp(:));
    sigma = nanstd(temp(:),0);
    z = bsxfun(@minus,temp(:), mu);
    z = bsxfun(@rdivide, z, sigma);
    data_mat_mean_all(sub,:,:)=reshape(z,8,8);
end

for i=1:8
for j=1:8
DependentVariable=data_mat_mean_all(Regressors_allnewTemp(:,5)==1  & Regressors_allnewTemp(:,6)<0.5,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1  & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,7]),Contrastnew(1:5),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data', Regressors_allnewTemp(:,5)',[], 1);
data=reshape(data_harmonized',size(data_harmonized',1),8,8);

for i=1:8
for j=1:8
DependentVariable=data(Regressors_allnewTemp(:,6)<0.5,i,j);
[b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T');
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;
end
end

Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,7]),1)-...
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(8,8);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;
T_all(:,:,step)=T_fdr;
end


files1=dir('I:\ADNI\MCI\Shortpath\*.mat');
files2=dir('I:\ADNI\NC_all_Shortpath\*.mat');
num=0;
ShortPath=[];
for i=1:length(files1)
    num
    num=num+1;
    load([files1(i).folder,'\',files1(i).name]);
    ShortPath(num,:,:)=P;
end
for j=1:length(files2)
    num
   num=num+1;
   load([files2(j).folder,'\',files2(j).name]);
    ShortPath(num,:,:)=P;
end
ShortPath(208,:,:)=[];
for sub=1:size(ShortPath,1)
    temp=squeeze(ShortPath(sub,:,:));
    temp(isinf(temp))=0;
    mu = nanmean(temp(:));
    sigma = nanstd(temp(:),0);
    z = bsxfun(@minus,temp(:), mu);
    z = bsxfun(@rdivide, z, sigma);
    ShortPath_zscrore(sub,:,:)=reshape(z,416,416);
end
ShortPath(isinf(ShortPath))=0;


for step=1:9
  S_mci_nc=squeeze(T_all(:,:,:,step));
  
   for sub=1:size(S_mci_nc,1)
%       temp= corrcoef(squeeze(S_mci_nc(sub,:,:)),squeeze(ShortPath(sub,:,:)));
%       corr_data(sub,step)=temp(1,2);
%         temp=squeeze(S_mci_nc(sub,:,:));
%         mu = nanmean(temp(:));
%         sigma = nanstd(temp(:),0);
%         z = bsxfun(@minus,temp(:), mu);
%         z = bsxfun(@rdivide, z, sigma);
%         S_mci_nc(sub,:,:)=reshape(z,416,416);
        
          D(sub,step) = mean(CalcHammingDist(squeeze(S_mci_nc(sub,:,:)),squeeze(P_all(sub,:,:)),2));
        
   end
 
end







%%

for i=1:size(S_all,4)
    for j=1:size(S_all,1)
        temp=squeeze(S_all(j,:,:,i));
        max_time(j,i)=sum(max(temp'))./416;
        
        
    end
end


Y=combat(mean_time',Regressors_allnewTemp(:,5),[],1);
Y=Y';

for j=1:size(S_all,4)
 DependentVariable=mean_time(Regressors_allnewTemp(:,5)==1,j);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end


for j=1:size(S_all,4)
 DependentVariable=Y(:,j);
        [b,r,SSE,SSR_1, ~, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp(:,[1,2,3,4,7]),Contrastnew(1:5),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end

for j=1:size(S_all,4)
 DependentVariable=Y(:,j);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp(:,[1,2,3,5,7]),Contrastnew(1:5),'T'); %YAN Chao-Gan 170714, Added Cohen's f squared (Effect Size) %[b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
end


for i=73:169
    i
    for j=1:490
       corr_value(i,j) = corr2(squeeze(S_mci_nc(i,:,:)),squeeze(S_mci_nc_ADNI(j,:,:)));
        
        
    end
end
for i=1:490
    for j=1:490
       corr_value(i,j) = corr2(squeeze(S_mci_nc_ADNI(i,:,:)),squeeze(S_mci_nc_ADNI(j,:,:)));
        
        
    end
end


for j=1:416
    for k=1:416

    if j~=k
        DependentVariable=S_mci_nc(Regressors_allnewTemp(:,5)==1,j,k);
        [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
            Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,4,7]),Contrastnew(1:5),'T'); 
         
    else
 
         TF_ForContrast(j,k)=0;

    end
end

end
T=TF_ForContrast;
T(isnan(T))=0;

mask=T;
mask(mask~=0)=1;

Ttemp=T(mask~=0);
df=size( Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,4,7]),1)-...
    size( Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,4,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
%     P   = 0.00001;
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(416,416);
Thresholdedmask(mask~=0)=Thresholded;

T_fdr=T.*Thresholdedmask;  




for i=1:416
    DependentVariable=times_in_Temp(Regressors_allnewTemp(:,5)==1,i);
     [b,r,SSE,SSR_1, T, TF_ForContrast_time_outTemp(i), Cohen_f2] =...
         y_regress_ss(DependentVariable,...
          Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,6,7]),Contrastnew(1:5),'T'); 
    
end
T=TF_ForContrast_time_outTemp';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;
Ttemp=T(mask~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,6,7]),1)-....
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,5)==1,[1,2,3,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr=T'.*Thresholdedmask;
find(T_fdr~=0);

%%
ELMCIexcel=xlsread('J:\zhengzihao_ELMCI_20241120\zhengzihao_ELMCI_20241120_3_20_2025.xlsx','B2:J2');

for i=1:length(ELMCIexcel)
    s1=ELMCIexcel{i,1};
    s2=ELMCIexcel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
%     if  str2num(s2{2})<10;
%         s2{2}=['0',s2{2}];
%     end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    ELMCIexcel{i,1}=[s1,'_',s2new];
    ELMCI_COV{i,1}=[s1,'_',s2new];
    if  strcmp(ELMCIexcel{i,2},'EMCI')
        ELMCI_COV{i,2}=-1;
    else
        ELMCI_COV{i,2}=1;
    end
    
    if ELMCIexcel{i,3}=='M'
        ELMCI_COV{i,3}=1;
    else
        ELMCI_COV{i,3}=2;
    end
    ELMCI_COV{i,4}=ELMCIexcel{i,4};
    
    
    
end
path='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);

[a,b1]=ismember(unique(fmri_name),ELMCI_COV(:,1));
b1(b1==0)=[];
Regressor_ELMCIsub1=ELMCI_COV([b1],:);
MCI_head_path='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_nii\preprocess\headmotion';
MCI_head_files=dir(MCI_head_path);
MCI_head_files(1:2)=[];
for i=1:length(MCI_head_files)

    cd([MCI_head_path,'\',MCI_head_files(i).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

    Regressor_ELMCIsub1{i,5}=mean(FD);
    
    cd(['J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open\',MCI_head_files(i).name]);
    
     files=dir('*.dcm');
    info=dicominfo([files(3).folder,'\',files(3).name]);    
    
    if strcmp(info.Manufacturer,'SIEMENS')
         Regressor_ELMCIsub1{i,6}=1;
    else
         Regressor_ELMCIsub1{i,6}=2;
    end

end
T1excel='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_T1_nii\Tissue_Volumes.csv';
T1_name=cellfun(@(name) strsplit(name,'/'),T1excel,'UniformOutput',false);
T1_name=cellfun(@(name) char(name(10)),T1_name,'UniformOutput',false);
T1_name=cellfun(@(name) name(1:end-11),T1_name,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name);
b1(b1==0)=[];
fmri_name=fmri_name(a);
T1volume= T1volume(b1);

[a,b2]=ismember(fmri_name,Regressor_ELMCIsub1(:,1));
Regressor_ELMCIsub1=Regressor_ELMCIsub1(b2,:);


files=dir('I:\ADNI\ELMCI\fmri_Eyes_Open_nii\Spreading\*.mat');
files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);
b4(b4==0)=[];
for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_elmci(i,:,:)=squeeze(T(:,:,7));
    
end

for i=1:length(b4)
    i
    load([files(i).folder,'\',files(b4(i)).name])
    S_elmci_all(i,:,:,:)=T(:,:,:);
    
end


 Regressor_ELMCIsub1=cell2mat(Regressor_ELMCIsub1(:,[2:6]));
 Regressor_ELMCIsub1=[Regressor_ELMCIsub1,T1volume,ones(135,1)];



path='J:\zhengzihao_ELMCI_20241120\resting_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);

T1excel='J:\zhengzihao_ELMCI_20241120\resting_nii_T1_nii\Tissue_Volumes.csv';
T1_name=cellfun(@(name) strsplit(name,'/'),T1excel,'UniformOutput',false);
T1_name=cellfun(@(name) char(name(10)),T1_name,'UniformOutput',false);
T1_name=cellfun(@(name) name(1:end-11),T1_name,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name);
b1(b1==0)=[];
fmri_name=fmri_name(a);
T1volume= T1volume(b1);

times=cellfun(@(name) name(1:10),fmri_name,'UniformOutput',false);
counts = zeros(size(times));             % 初始化结果数组
counterMap = containers.Map();                % 创建映射来计数
for i = 1:length(times)
    key = times{i};                      % 当前字符串
    if isKey(counterMap, key)
        counterMap(key) = counterMap(key) + 1; % 已出现过，计数加一
    else
        counterMap(key) = 1;                   % 首次出现，初始化为1
    end
    counts(i) = counterMap(key);               % 存储当前是第几次出现
end


[a,b2]=ismember(unique(fmri_name),ELMCI_COV(:,1));
Regressor_ELMCIsub2=ELMCI_COV([b2],:);

MCI_head_path='J:\zhengzihao_ELMCI_20241120\resting_nii\preprocess\headmotion';
MCI_head_files=dir(MCI_head_path);
MCI_head_files(1:2)=[];
MCI_head_name=cellfun(@(name) name(1:end-11),{MCI_head_files.name},'UniformOutput',false);
[a,b3]=ismember(fmri_name,MCI_head_name);
for i=1:length(fmri_name)

    cd([MCI_head_path,'\',MCI_head_files(b3(i)).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

    Regressor_ELMCIsub2{i,5}=mean(FD);

end
Regressor_ELMCIsub2=cell2mat(Regressor_ELMCIsub2(:,[2:5]));
Regressor_ELMCIsub2=[Regressor_ELMCIsub2,T1volume,ones(286,1)];
 

files=dir('I:\ADNI\ELMCI\resting_nii\Spreading\*.mat');

files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);

for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_mci_nc(i,:,:)=squeeze(T(:,:,7));
    
end
for i=1:length(b4)
    i
    load([files(i).folder,'\',files(b4(i)).name])
    S_all(i,:,:,:)=T;
    
end

files=dir('I:\ADNI\ELMCI\fmri_Eyes_Open_nii\Spreading\*.mat');
S_mci_nc=[];
for i=1:length(files)
    load([files(i).folder,'\',files(i).name])
    S_mci_nc(i,:,:)=squeeze(T(:,:,7));
    
end




clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for step =1:9
S_mci_nc=squeeze(T_all(:,:,:,step));
S_mci_nc=S_mci_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);


for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==32)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end



% for sub=1:size(S_mci_nc,1)
%     temp=squeeze(S_mci_nc(sub,:,:));
% 

%     maxV=max(temp(:));
%     minV=min(temp(:));
%     S_mci_nc(sub,:,:)=(temp-minV)./(maxV-minV); 
% end

new=S_mci_nc_temp(Regressor_ELMCIsub2(:,4)<0.5,:,:);
new(isnan(new))=0;
cInd=find(Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,1)==-1);
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data = myzscore(data,0,2);
data = zScoreToSubset(data,cInd);
data=reshape(data,size(data,1),18,18);

% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,18,18);
% end

for i=1:18
for j=1:18
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),[1,0,0,0,0,0],'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
    df=size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),1)-...
        size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;

T_all(:,:,step)=T_fdr;
end





new=S_mci_nc_temp;
new(isnan(new))=0;
T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(Regressor_ELMCIsub2(:,4)<0.5,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),[1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),1)-...
        size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
for i=1:18
    
    try
    plot_hemispherewithboundarynew_zhengnew(T_fdr(yeo_17labels,i),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
    print(gcf,['I:\ADNI\yeonet_',num2str(i),'.png'],'-dpng','-r600');
    close all;
    
    catch
    end
end
T2map=[];
for net=1:18
    network=new(:,:,yeoindex17==net);

    network=mean(network,3);
    network=squeeze(network);
   for i=1:416
        DependentVariable=network(Regressor_ELMCIsub2(:,4)<0.5,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),[1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),1)-...
        size( Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end

for j=1:416
for k=1:416

if j~=k
    DependentVariable=new(Regressor_ELMCIsub2(:,4)<0.5,j,k);
    [b,r,SSE,SSR_1, T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),[1,0,0,0,0,0],'T'); 

else

     TF_ForContrast(j,k)=0;

end
end
end





new=S_mci_nc_temp;
new(isnan(new))=0;

for i=1:18
    
    try
    plot_hemispherewithboundarynew_zhengnew(T_fdr(yeo_17labels_17_to_7,i),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
%     print(gcf,['I:\ADNI\Combatyeonet_',num2str(i),'.png'],'-dpng','-r600');
%     close all;
    
    catch
    end
end

net=18;
network=new(:,yeoindex17==net,:);
network=mean(network,2);
network=squeeze(network);
DependentVariable=network(Regressor_ELMCIsub2(:,4)<0.5,:);
[a,b]=sortrows(Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5,:),[1,3]);

idx=[];
for i=1:18
    
    idx=[idx;find(yeoindex17==i)];
end

DependentVariable=DependentVariable(b,idx);
imagesc(DependentVariable');



NCexcel=xlsread('J:\zhengzihao_ELMCI_20241120\NC\zhengzihao_NC_resting_5_16_2025.csv','B2:J2');

for i=1:length(NCexcel)
    s1=NCexcel{i,1};
    s2=NCexcel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
%     if  str2num(s2{2})<10;
%         s2{2}=['0',s2{2}];
%     end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    NCexcel{i,1}=[s1,'_',s2new];
    NCexcel_COV{i,1}=[s1,'_',s2new];
    if  strcmp(NCexcel{i,2},'EMCI')
        NCexcel_COV{i,2}=-1;
    else
        NCexcel_COV{i,2}=1;
    end
    
    if NCexcel{i,3}=='M'
        NCexcel_COV{i,3}=1;
    else
        NCexcel_COV{i,3}=2;
    end
    NCexcel_COV{i,4}=NCexcel{i,4};
    
    
    
end

path='J:\zhengzihao_ELMCI_20241120\NC\resting_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);


T1excel_NC='J:\zhengzihao_ELMCI_20241120\resting_nii_T1_nii\Tissue_Volumes.csv';

T1_name_NC=cellfun(@(name) strsplit(name,'/'),T1excel_NC,'UniformOutput',false);
T1_name_NC=cellfun(@(name) char(name(10)),T1_name_NC,'UniformOutput',false);
T1_name_NC=cellfun(@(name) name(1:end-11),T1_name_NC,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name_NC);
b1(b1==0)=[];
fmri_name=fmri_name(a);
T1volume_NC= T1volume_NC(b1);
times=cellfun(@(name) name(1:10),fmri_name,'UniformOutput',false);
countsNC = zeros(size(times));             % 初始化结果数组
counterMap = containers.Map();                % 创建映射来计数
for i = 1:length(times)
    key = times{i};                      % 当前字符串
    if isKey(counterMap, key)
        counterMap(key) = counterMap(key) + 1; % 已出现过，计数加一
    else
        counterMap(key) = 1;                   % 首次出现，初始化为1
    end
    countsNC(i) = counterMap(key);               % 存储当前是第几次出现
end





[a,b2]=ismember(unique(fmri_name),NCexcel_COV(:,1));
Regressor_NC=NCexcel_COV([b2],:);

NC_head_path='J:\zhengzihao_ELMCI_20241120\NC\resting_nii\preprocess\headmotion';
NC_head_files=dir(NC_head_path);
NC_head_files(1:2)=[];
NC_head_name=cellfun(@(name) name(1:end-11),{NC_head_files.name},'UniformOutput',false);
[a,b3]=ismember(fmri_name,NC_head_name);
for i=1:length(fmri_name)

    cd([NC_head_path,'\',NC_head_files(b3(i)).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

   Regressor_NC{i,5}=mean(FD);

end

 Regressor_NC=cell2mat(Regressor_NC(:,[2:5]));
 Regressor_NC=[ Regressor_NC,T1volume_NC,ones(152,1)];
 
files=dir('J:\zhengzihao_ELMCI_20241120\NC\Spreading\*.mat');

files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);

for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_nc(i,:,:)=squeeze(T(:,:,7));
    
end
for i=1:length(files_name)
    load([files(i).folder,'\',files(i).name])
    S_nc(i,:,:)=squeeze(T(:,:,7));
    
end
num=286;
for i=1:length(files_name)
    load([files(i).folder,'\',files(i).name])
    num=num+1;
    S_all(num,:,:,:)=T(:,:,:);
    
end


S_nc=S_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);

for i=1:size(S_nc,1)
    temp=squeeze(S_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_nc_temp(i,:,:)=ZScoreMatrix(temp);
end
 
new=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=2 ,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5,:,:));
new(isnan(new))=0;

Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=2) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=2)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];

    for j=1:416
        for k=1:416

        if j~=k
            DependentVariable=new(:,j,k);
         [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(j,k), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
        else

             TF_ForContrast(j,k)=0;

        end
    end

    end




% cInd=find(Regressor(:,1)==-1);
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data = myzscore(data,0,2);
data = zScoreToSubset(data,cInd);
data=reshape(data,size(data,1),18,18);

% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,18,18);
% end

for i=1:18
for j=1:18
DependentVariable=data(:,i,j);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
  Regressor,[1,0,0,0,0,0],'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
    df=size(  Regressor,1)-...
        size(  Regressor,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;


T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(:,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end

%%
NC='I:\ADNI\MCI\BOLD_Smooth_new\*.mat';

MCI='I:\ADNI\NC_All_BOLD_Smooth_new\*.mat';


NCdir=dir(NC);
MCIdir=dir(MCI);
enrichedFC=[];
receptors_name={'5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2','CB1', 'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5','MOR', 'NET', 'NMDA', 'VAChT'};
E_name={ '5HT2a', '5HT4', '5HT6', 'A4B2', 'D1', 'M1', 'mGluR5', 'NMDA'};
I_name={'5HT1a', '5HT1b','CB1', 'D2', 'GABAa', 'H3','MOR'};
[~,E_index]=ismember(E_name,receptors_name);
[~,I_index]=ismember(I_name,receptors_name);
receptors_E=mean(zscore(receptor(:,E_index)),2);
receptors_I=mean(zscore(receptor(:,I_index)),2);
for type=1:19
   type
    num=0;
for i=1:length(MCIdir)
    num=num+1;
    load([MCIdir(i).folder,'\',MCIdir(i).name])
    
        BOLD=theROITimeCoursesTotal;
        beta1=[];
            for time=1:size(BOLD,1)

        %         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
                [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(time,:))',[zscore(receptor(:,type)),ones(416,1)],[1,0],'T');
                beta1(time)=b(1);

            end
        beta2=[];
           for roi=1:416

        %         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
                [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(:,roi)),[zscore(beta1'),ones(size(BOLD,1),1)],[1,0],'T');
                beta2(roi)=b(1);

            end  
    enrichedFC(num,:,type)=zscore(beta2);
    
end
for i=1:length(NCdir)
    num=num+1;
    load([NCdir(i).folder,'\',NCdir(i).name])
    
        BOLD=theROITimeCoursesTotal;
        beta1=[];
            for time=1:size(BOLD,1)

        %         lm=fitglm(zscore(receptor(:,3)),zscore(BOLD(time,:))');
                [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(time,:))',[zscore(receptor(:,type)),ones(416,1)],[1,0],'T');
                beta1(time)=b(1);

            end
        beta2=[];
           for roi=1:416

        %         lm=fitglm(zscore(ASL(i,:)'),zscore(BOLD(time,:)));
                [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(zscore(BOLD(:,roi)),[zscore(beta1'),ones(size(BOLD,1),1)],[1,0],'T');
                beta2(roi)=b(1);

            end  
    enrichedFC(num,:,type)=zscore(beta2);
    
end
end
% enrichedFC(208,:,:)=[];
data_harmonized = combat(squeeze(enrichedFC(:,:,6))', Regressors_allnewTemp(:,5)',[], 0);
for i=1:416
    DependentVariable=data_harmonized(i,Regressors_allnewTemp(:,6)<0.5)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
end

Ttemp=TF_ForContrast';
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdrenriched=Ttemp.*Thresholded;

for type=1:19
for i=1:416
    DependentVariable=squeeze(enrichedFC(Regressors_allnewTemp(:,5)==1 & Regressors_allnewTemp(:,6)<0.5,i,type));
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regress_Simens,Contrastnew(1:6),'T'); 
end
Ttemp=TF_ForContrast';
df=size(Regress_Simens,1)-...
size(Regress_Simens,2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdrenriched(:,type)=Ttemp.*Thresholded;


end


for type=1:19
    
data_harmonized = combat(squeeze(enrichedFC(:,:,type))', Regressors_allnewTemp(:,5)',[], 0);
for i=1:416
    DependentVariable=data_harmonized(i,Regressors_allnewTemp(:,6)<0.5)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
end
T_fdrenrichedmap(:,type)=TF_ForContrast';
Ttemp=TF_ForContrast';
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));

if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdrenriched(:,type)=Ttemp.*Thresholded;


end

h=[];
h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.3,1]);
allw=T_fdrenrichedmap;
vmax=5;
    for j=1:19
         vertices_l=[left_surface.vertices];
        vertices_l=double(vertices_l);
        faces_l=double([left_surface.faces]);
        
        vertices_r=[right_surface.vertices];
        vertices_r=double(vertices_r);
        faces_r=double([right_surface.faces]);
        global BOUNDARY_L BOUNDARY_R
           T=allw(:,j);
%            T=T(yeo_17labels_17_to_7);
            T(isnan(T))=0;
            T_pos=T;
            T_pos(T<=0)=0;
            T_neg=T;
            T_neg(T>=0)=0;

            resultMap=zeros(size(brainmask,1),1);
            for i=1:400
                resultMap(brainmask==i)=T(i);
            end
            resultsubMap=zeros(16,1);
            for i=401:416
                resultsubMap(i-400)=T(i);
            end
               brainmask_new=brainmask;
                   brainmask_new(resultMap==0)=0;
                    
                  if size(unique(brainmask_new(1:32492)),1)==1
                      BOUNDARY_L={};
                  else
                       [BOUNDARY_L,BOUNDARY_ROI_ID_L] = findROIboundaries(vertices_l,faces_l,brainmask_new(1:32492));
                  end
                  if size(unique(brainmask_new(32493:end)),1)==1
                      BOUNDARY_R={};
                  else
                       [BOUNDARY_R,BOUNDARY_ROI_ID_R] = findROIboundaries(vertices_r,faces_r, brainmask_new(32493:end));
                  end


                    maxpos=vmax;
                    minpos=0;
                    minneg=-vmax;
                    maxneg=0;
            resultcol=zeros(size(resultMap,1),3);
            if  maxpos==minpos
                if find(T==maxpos)<=400
                    flagpos=1;
                else
                    flagpos=2;
                end
                idx_pos=100;
                pos_step=1;
    %                         poscolor=(othercolor('YlOrRd9',100)); 
                resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

            else
                flagpos=3;
                pos_step=((maxpos)-(minpos))/100;
                idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                idx_pos(idx_pos<1)=1;
                idx_pos(idx_pos>100)=100;
    %                         poscolor=(othercolor('YlOrRd9',100));  
                resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
            end

           if  maxneg==minneg
                if find(T==maxneg)<=400
                    flagneg=1;
                else
                    flagneg=2;
                end
                idx_neg=1;
                neg_step=1;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

            else
                flagneg=3;
                neg_step=((-minneg)-(-maxneg))/100;
                idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                idx_neg(idx_neg<1)=1;
                idx_neg(idx_neg>100)=100;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
           end
            resultcol(resultMap==0,:)=repmat([0.75 0.75 0.75],sum(resultMap==0),1);


            S = convert_surface(surface);
            data=resultcol;
            D={};
            for ii = 1:2
            D{1,ii} = data(1:size(S{1}.coord,2),:);
            if numel(S) == 2
                D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
            end
            end
            jj=1;
                for ii = 1:numel(S)*2
                    idx = ceil(ii/2);
                    h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-(1/20)*j-0.02 .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                    h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                ,'FaceColor', 'interp','EdgeColor', 'none');
%                                       if ii==1|ii==2
%                             for bounday_i = 1:length(BOUNDARY_L)
%                                        hold on;
%                               plot3(BOUNDARY_L{bounday_i}(:,1), BOUNDARY_L{bounday_i}(:,2), BOUNDARY_L{bounday_i}(:,3), 'Color', 'w', 'LineWidth',0.5,'Clipping','off');
%                             end
%                              hold off;
%                        else
%                              for bounday_i = 1:length(BOUNDARY_R)
%                                        hold on;
%                               plot3(BOUNDARY_R{bounday_i}(:,1), BOUNDARY_R{bounday_i}(:,2), BOUNDARY_R{bounday_i}(:,3), 'Color', 'w', 'LineWidth',0.5,'Clipping','off');
%                              end
%                              hold off;
%                                
%                         end
            
   
                    material dull; lighting phong;
                end
                set(h.axes(:,jj)                    , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]      )
            set(h.axes(1,:),'View',[-90 0]);
            set(h.axes(2,:),'View',[90 0]);
            set(h.axes(3,:),'View',[-90 0]);
            set(h.axes(4,:),'View',[90 0]);


            for ii = 1:4
                axes(h.axes(ii));
                h.camlight(ii) = camlight();
            end

           idx_subcluster = find(resultsubMap ~= 0);
            idx_subcluster_out=setdiff([1:16],idx_subcluster);

            for ii=5:6
                for jj=1
                    h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-(1/20)*j-0.02 .1 .1]);

                    for i=1:length(idx_subcluster)
                        thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                        temp1=thala.dat;
                        temp1(temp1~=idx_subcluster(i))=0;
                        thala.dat=temp1;
                        thalaregion=region(thala);
                        if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                            subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                            subidx_pos(subidx_pos<1)=1;
                            subidx_pos(subidx_pos>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                            subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                            subidx_neg(subidx_neg<1)=1;
                            subidx_neg(subidx_neg>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                        end


                    end
                for i=1:length(idx_subcluster_out)
                    thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                    temp=thala.dat;
                    temp(temp~=idx_subcluster_out(i))=0;
                    thala.dat=temp;
                    thalaregion=region(thala);
                    p = imageCluster('cluster',region2struct(thalaregion),'color',[0.75,0.75,0.75],'alpha',1);
                end

                  material dull; lighting phong;
                set(h.axes(:,jj)               , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  ); 
                end
            end
%             set(h.axes(5,1),'View',[-90 0]);
%             set(h.axes(6,1),'View',[90 0]);   
            
            set(h.axes(5,1),'View',[45 0]);
            set(h.axes(6,1),'View',[-45 0]);     
       
%             set(h.axes(7,1),'View',[45 0]);    
%             set(h.axes(8,1),'View',[-45 0]);  
            material dull; lighting phong;
                set(gca                 , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  );
            set(gca,'color','w');
            set(gcf,'color','w');



    end


    
h=[];
h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.3,1]);
allw=T_fdrenriched;
    for j=1:19
         vertices_l=[left_surface.vertices];
        vertices_l=double(vertices_l);
        faces_l=double([left_surface.faces]);
        
        vertices_r=[right_surface.vertices];
        vertices_r=double(vertices_r);
        faces_r=double([right_surface.faces]);
        global BOUNDARY_L BOUNDARY_R
           T=allw(:,j);
%            T=T(yeo_17labels_17_to_7);
            T(isnan(T))=0;
            T_pos=T;
            T_pos(T<=0)=0;
            T_neg=T;
            T_neg(T>=0)=0;

            resultMap=zeros(size(brainmask,1),1);
            for i=1:400
                resultMap(brainmask==i)=T(i);
            end
            resultsubMap=zeros(16,1);
            for i=401:416
                resultsubMap(i-400)=T(i);
            end
               brainmask_new=brainmask;
                   brainmask_new(resultMap==0)=0;
                    
                  if size(unique(brainmask_new(1:32492)),1)==1
                      BOUNDARY_L={};
                  else
                       [BOUNDARY_L,BOUNDARY_ROI_ID_L] = findROIboundaries(vertices_l,faces_l,brainmask_new(1:32492));
                  end
                  if size(unique(brainmask_new(32493:end)),1)==1
                      BOUNDARY_R={};
                  else
                       [BOUNDARY_R,BOUNDARY_ROI_ID_R] = findROIboundaries(vertices_r,faces_r, brainmask_new(32493:end));
                  end


                    maxpos=max(T_pos);
                    minpos=min(T_pos(T_pos~=0));
                    minneg=min(T_neg);
                    maxneg=max(T_neg(T_neg~=0));
            resultcol=zeros(size(resultMap,1),3);
            if  maxpos==minpos
                if find(T==maxpos)<=400
                    flagpos=1;
                else
                    flagpos=2;
                end
                idx_pos=100;
                pos_step=1;
    %                         poscolor=(othercolor('YlOrRd9',100)); 
                resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

            else
                flagpos=3;
                pos_step=((maxpos)-(minpos))/100;
                idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                idx_pos(idx_pos<1)=1;
                idx_pos(idx_pos>100)=100;
    %                         poscolor=(othercolor('YlOrRd9',100));  
                resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
            end

           if  maxneg==minneg
                if find(T==maxneg)<=400
                    flagneg=1;
                else
                    flagneg=2;
                end
                idx_neg=1;
                neg_step=1;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

            else
                flagneg=3;
                neg_step=((-minneg)-(-maxneg))/100;
                idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                idx_neg(idx_neg<1)=1;
                idx_neg(idx_neg>100)=100;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
           end
            resultcol(resultMap==0,:)=repmat([0.75 0.75 0.75],sum(resultMap==0),1);


            S = convert_surface(surface);
            data=resultcol;
            D={};
            for ii = 1:2
            D{1,ii} = data(1:size(S{1}.coord,2),:);
            if numel(S) == 2
                D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
            end
            end
            jj=1;
                for ii = 1:numel(S)*2
                    idx = ceil(ii/2);
                    h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-(1/20)*j-0.02 .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                    h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                ,'FaceColor', 'interp','EdgeColor', 'none');
%                                       if ii==1|ii==2
%                             for bounday_i = 1:length(BOUNDARY_L)
%                                        hold on;
%                               plot3(BOUNDARY_L{bounday_i}(:,1), BOUNDARY_L{bounday_i}(:,2), BOUNDARY_L{bounday_i}(:,3), 'Color', 'w', 'LineWidth',0.5,'Clipping','off');
%                             end
%                              hold off;
%                        else
%                              for bounday_i = 1:length(BOUNDARY_R)
%                                        hold on;
%                               plot3(BOUNDARY_R{bounday_i}(:,1), BOUNDARY_R{bounday_i}(:,2), BOUNDARY_R{bounday_i}(:,3), 'Color', 'w', 'LineWidth',0.5,'Clipping','off');
%                              end
%                              hold off;
%                                
%                         end
            
   
                    material dull; lighting phong;
                end
                set(h.axes(:,jj)                    , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]      )
            set(h.axes(1,:),'View',[-90 0]);
            set(h.axes(2,:),'View',[90 0]);
            set(h.axes(3,:),'View',[-90 0]);
            set(h.axes(4,:),'View',[90 0]);


            for ii = 1:4
                axes(h.axes(ii));
                h.camlight(ii) = camlight();
            end

           idx_subcluster = find(resultsubMap ~= 0);
            idx_subcluster_out=setdiff([1:16],idx_subcluster);

            for ii=5:8
                for jj=1
                    h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-(1/20)*j-0.02 .1 .1]);

                    for i=1:length(idx_subcluster)
                        thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                        temp1=thala.dat;
                        temp1(temp1~=idx_subcluster(i))=0;
                        thala.dat=temp1;
                        thalaregion=region(thala);
                        if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                            subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                            subidx_pos(subidx_pos<1)=1;
                            subidx_pos(subidx_pos>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                            subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                            subidx_neg(subidx_neg<1)=1;
                            subidx_neg(subidx_neg>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                        end


                    end
                for i=1:length(idx_subcluster_out)
                    thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                    temp=thala.dat;
                    temp(temp~=idx_subcluster_out(i))=0;
                    thala.dat=temp;
                    thalaregion=region(thala);
                    p = imageCluster('cluster',region2struct(thalaregion),'color',[0.75,0.75,0.75],'alpha',1);
                end

                  material dull; lighting phong;
                set(h.axes(:,jj)               , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  ); 
                end
            end
            set(h.axes(5,1),'View',[-90 0]);
            set(h.axes(6,1),'View',[90 0]);     
            set(h.axes(7,1),'View',[45 0]);    
            set(h.axes(8,1),'View',[-45 0]);  
            material dull; lighting phong;
                set(gca                 , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  );
            set(gca,'color','w');
            set(gcf,'color','w');



    end
    
    


T_fdrenriched=Ttemp.*Thresholded;
%%

load('I:\ADNI\T_all.mat')
S_mci_nc=squeeze(T_all(:,:,:,7));
clear T_all

load('I:\ADNI\NCI_NC_COV_NEW_489.mat')
Regressors_allnewTemp=cell2mat(MCI_NC_COV_new(:,2:end));
Regressors_allnewTemp=[Regressors_allnewTemp,ones(size(Regressors_allnewTemp,1),1)];

S_mci_nc=S_mci_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for i=1:size(S_mci_nc,1)
    i
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end
S_nc=S_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for i=1:size(S_nc,1)
    temp=squeeze(S_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_nc_temp(i,:,:)=ZScoreMatrix(temp);
end
S_


S_ELMCI=S_ELMCI(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for i=1:size(S_ELMCI,1)
    temp=squeeze(S_ELMCI(i,:,:));
    temp(temp==1 | temp==32)=NaN;
    S_ELMCI_temp(i,:,:)=ZScoreMatrix(temp);
end
T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
    network1=S_mci_nc_temp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
    
    
    network2=S_ELMCI_temp(:,yeoindex17==net,:);
    network2=mean(network2,2);
    network2=squeeze(network2);
    

    network3=S_nc_temp(:,yeoindex17==net,:);
    network3=mean(network3,2);
    network3=squeeze(network3);
    
    batch=[Regressors_allnewTemp(:,5)',3*ones(1,size(Regressor_ELMCIsub2,1)),3*ones(1,size(Regressor_NC,1))];
    dataall=[network1;network2;network3];
    datanew=combat(dataall',batch,[], 0);
    
    datanew_mci_nc=datanew(:,batch==1 | batch==2);
    datanew_elmci=datanew(:,490:490+285);
    datanew_nc=datanew(:,490+286:end);
    
      for i=1:416
        DependentVariable=datanew_mci_nc(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_2=[T2map_2;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_2(:,net)=Ttemp.*Thresholded;
    
    
    
     new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_3=[T2map_3;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_3(:,net)=Ttemp.*Thresholded;
    
end
real_data1=T2map_1(:,yeo_17labels_17_to_7);
real_data2=T2map_2(:,yeo_17labels_17_to_7);
real_data3=T2map_3(:,yeo_17labels_17_to_7);
[U1,S1,V1]=svd(real_data1);
[U2,S2,V2]=svd(real_data2);
[U3,S3,V3]=svd(real_data3);



T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
    network1=S_mci_nc_temp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
    
    
    network2=S_ELMCI_temp(:,yeoindex17==net,:);
    network2=mean(network2,2);
    network2=squeeze(network2);
    

    network3=S_nc_temp(:,yeoindex17==net,:);
    network3=mean(network3,2);
    network3=squeeze(network3);
    
    batch=[Regressors_allnewTemp(:,5)',3*ones(1,size(Regressor_ELMCIsub2,1)),3*ones(1,size(Regressor_NC,1))];
    dataall=[network1;network2;network3];
    datanew=combat(dataall',batch,[], 0);
    
    datanew_mci_nc=datanew(:,batch==1 | batch==2);
    datanew_elmci=datanew(:,490:490+285);
    datanew_nc=datanew(:,490+286:end);
    
      for i=1:416
        DependentVariable=datanew_mci_nc(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_2=[T2map_2;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_2(:,net)=Ttemp.*Thresholded;
    
    
    
     new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_3=[T2map_3;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_3(:,net)=Ttemp.*Thresholded;
    
end
real_data1=T2map_1(:,yeo_17labels_17_to_7);
real_data2=T2map_2(:,yeo_17labels_17_to_7);
real_data3=T2map_3(:,yeo_17labels_17_to_7);
[coeff, score1, latent, tsquared, explained, mu] = pca(real_data1');
[coeff, score2, latent, tsquared, explained, mu] = pca(real_data2');
[coeff, score3, latent, tsquared, explained, mu] = pca(real_data3');





T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
    network1=S_mci_nc_temp(:,:,yeoindex17==net);

    network1=mean(network1,3);
    network1=squeeze(network1);
    
    
    network2=S_ELMCI_temp(:,:,yeoindex17==net);
    network2=mean(network2,3);
    network2=squeeze(network2);
    

    network3=S_nc_temp(:,:,yeoindex17==net);
    network3=mean(network3,3);
    network3=squeeze(network3);
    
    batch=[Regressors_allnewTemp(:,5)',3*ones(1,size(Regressor_ELMCIsub2,1)),3*ones(1,size(Regressor_NC,1))];
    dataall=[network1;network2;network3];
    datanew=combat(dataall',batch,[], 0);
    
    datanew_mci_nc=datanew(:,batch==1 | batch==2);
    datanew_elmci=datanew(:,490:490+285);
    datanew_nc=datanew(:,490+286:end);
    
      for i=1:416
        DependentVariable=datanew_mci_nc(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_2=[T2map_2;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_2(:,net)=Ttemp.*Thresholded;
    
    
    
     new=cat(2,datanew_elmci(:,Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1),...
    datanew_nc(:,Regressor_NC(:,4)<0.5));
    new(isnan(new))=0;

    Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1) ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5,:)];

    Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==-1)),1);...
   -1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
    
     for i=1:416
        DependentVariable=new(i,:)';
    [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
     end
    Ttemp=TF_ForContrast';
        T2map_3=[T2map_3;TF_ForContrast];
    df=size(  Regressor,1)-...
        size(  Regressor,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_3(:,net)=Ttemp.*Thresholded;
    
end
real_data1=T2map_1(:,yeo_17labels_17_to_7);
real_data2=T2map_2(:,yeo_17labels_17_to_7);
real_data3=T2map_3(:,yeo_17labels_17_to_7);
[coeff, score1, latent, tsquared, explained, mu] = pca(real_data1');
[coeff, score2, latent, tsquared, explained, mu] = pca(real_data2');
[coeff, score3, latent, tsquared, explained, mu] = pca(real_data3');



dataall=[reshape(S_mci_nc_temp,489,[]);reshape(S_ELMCI_temp,286,[]);reshape(S_nc_temp,152,[])];
dataallnew=dataall;
dataallnew(:,find(eye(416)))=[];
datanew=combat(dataallnew', [Regressors_allnewTemp(:,5)',3*ones(1,size(Regressor_ELMCIsub2,1)),3*ones(1,size(Regressor_NC,1))],[], 0);


%%
Regressor_ELMCIsub1=Regressor_ELMCIsub1(:,[1,2,3,6,5,4,7]);
Regressor_ELMCIsub1(:,1)=ones(size(Regressor_ELMCIsub1,1),1);

Regressors_allnewTempNC=Regressors_allnewTemp(Regressors_allnewTemp(:,1)==-1,:);

S_mci_nc_tempNC=S_mci_nc_temp(Regressors_allnewTemp(:,1)==-1,:,:);

S_elmci=S_elmci(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for i=1:size(S_elmci,1)
temp=squeeze(S_elmci(i,:,:));
temp(temp==1 | temp==32)=NaN;
S_elmci_temp(i,:,:)=ZScoreMatrix(temp);
end
Datatemp=cat(1,S_elmci_temp,S_mci_nc_tempNC);
Regressorstemp=[Regressor_ELMCIsub1;Regressors_allnewTempNC];

T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
network1=Datatemp(:,yeoindex17==net,:);

    network1=mean(network1,2);
    network1=squeeze(network1);
    
    data_harmonized = combat(network1', Regressorstemp(:,5)',[], 0);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    

    
end





%%
T2map_1=[];
T2map_2=[];
T2map_3=[];
age_gt= 50
age_lt= 80
Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:,:));

RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:)];
RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt)),1);...
    -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt)),1)];

T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
    
    
    
    
    network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
    
    
 
    
      for i=1:416
        DependentVariable=network1(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
       RegressorsTemp,[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(RegressorsTemp,1)-...
        size(RegressorsTemp,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    
end


T2map_1=[];
T2map_2=[];
T2map_3=[];
age_gt= 50
age_lt= 80
Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & Regressor_ELMCIsub2(:,1)==1,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:,:));

RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & Regressor_ELMCIsub2(:,1)==1,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:)];
RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & Regressor_ELMCIsub2(:,1)==1)),1);...
    -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt)),1)];

for net=1:18
    

    network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);

    
      for i=1:416
        DependentVariable=network1(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
       RegressorsTemp,[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(RegressorsTemp,1)-...
        size(RegressorsTemp,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    
end
real_data1=T2map_1(:,yeo_17labels_17_to_7);
[coeff, score1, latent, tsquared, explained, mu] = pca(real_data1');
plot_hemisphere_zhengnew(score1(:,1),{left_surface,right_surface},brainmask,1,poscolor,negcolor);
plot_hemisphere_zhengnew(mean(real_data1),{left_surface,right_surface},brainmask,1,poscolor,negcolor);





%%
close all;
h=[];
h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.3,1]);
allw=T_fdr;
    for j=1:18

           T=allw(:,j);
           T=T(yeo_17labels_17_to_7);
            T(isnan(T))=0;
            T_pos=T;
            T_pos(T<=0)=0;
            T_neg=T;
            T_neg(T>=0)=0;

            resultMap=zeros(size(brainmask,1),1);
            for i=1:400
                resultMap(brainmask==i)=T(i);
            end
            resultsubMap=zeros(16,1);
            for i=401:416
                resultsubMap(i-400)=T(i);
            end

            maxpos=max(T_pos);
            minpos=min(T_pos(T_pos~=0));
            minneg=min(T_neg);
            maxneg=max(T_neg(T_neg~=0));
            resultcol=zeros(size(resultMap,1),3);
            if  maxpos==minpos
                if find(T==maxpos)<=400
                    flagpos=1;
                else
                    flagpos=2;
                end
                idx_pos=100;
                pos_step=1;
    %                         poscolor=(othercolor('YlOrRd9',100)); 
                resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

            else
                flagpos=3;
                pos_step=((maxpos)-(minpos))/100;
                idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                idx_pos(idx_pos<1)=1;
                idx_pos(idx_pos>100)=100;
    %                         poscolor=(othercolor('YlOrRd9',100));  
                resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
            end

           if  maxneg==minneg
                if find(T==maxneg)<=400
                    flagneg=1;
                else
                    flagneg=2;
                end
                idx_neg=1;
                neg_step=1;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

            else
                flagneg=3;
                neg_step=((-minneg)-(-maxneg))/100;
                idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                idx_neg(idx_neg<1)=1;
                idx_neg(idx_neg>100)=100;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
           end
            resultcol(resultMap==0,:)=repmat([0.75 0.75 0.75],sum(resultMap==0),1);


            S = convert_surface(surface);
            data=resultcol;
            D={};
            for ii = 1:2
            D{1,ii} = data(1:size(S{1}.coord,2),:);
            if numel(S) == 2
                D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
            end
            end
            jj=1;
                for ii = 1:numel(S)*2
                    idx = ceil(ii/2);
                    h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-(1/19)*j-0.02 .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                    h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                ,'FaceColor', 'interp','EdgeColor', 'none');

                    material dull; lighting phong;
                end
                set(h.axes(:,jj)                    , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]      )
            set(h.axes(1,:),'View',[-90 0]);
            set(h.axes(2,:),'View',[90 0]);
            set(h.axes(3,:),'View',[-90 0]);
            set(h.axes(4,:),'View',[90 0]);


            for ii = 1:4
                axes(h.axes(ii));
                h.camlight(ii) = camlight();
            end

           idx_subcluster = find(resultsubMap ~= 0);
            idx_subcluster_out=setdiff([1:16],idx_subcluster);

            for ii=5:8
                for jj=1
                    h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-(1/19)*j-0.02 .1 .1]);

                    for i=1:length(idx_subcluster)
                        thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                        temp1=thala.dat;
                        temp1(temp1~=idx_subcluster(i))=0;
                        thala.dat=temp1;
                        thalaregion=region(thala);
                        if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                            subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                            subidx_pos(subidx_pos<1)=1;
                            subidx_pos(subidx_pos>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                            subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                            subidx_neg(subidx_neg<1)=1;
                            subidx_neg(subidx_neg>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                        end


                    end
                   for i=1:length(idx_subcluster_out)
                    thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                    temp=thala.dat;
                    temp(temp~=idx_subcluster_out(i))=0;
                    thala.dat=temp;
                    thalaregion=region(thala);
                    p = imageCluster('cluster',region2struct(thalaregion),'color',[0.75,0.75,0.75],'alpha',1);
                   end

                  material dull; lighting phong;
                set(h.axes(:,jj)               , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  ); 
                end
            end
            set(h.axes(5,1),'View',[-90 0]);
            set(h.axes(6,1),'View',[90 0]);     
            set(h.axes(7,1),'View',[45 0]);    
            set(h.axes(8,1),'View',[-45 0]);  
            material dull; lighting phong;
                set(gca                 , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  );
            set(gca,'color','w');
            set(gcf,'color','w');



    end
    
    
    
    
    close all;
h=[];
h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.3,1]);
allw=T_fdr;
    for j=1:18

           T=allw(:,j);
           T=T(yeo_17labels_17_to_7);
            T(isnan(T))=0;
            T_pos=T;
            T_pos(T<=0)=0;
            T_neg=T;
            T_neg(T>=0)=0;

            resultMap=zeros(size(brainmask,1),1);
            for i=1:400
                resultMap(brainmask==i)=T(i);
            end
            resultsubMap=zeros(16,1);
            for i=401:416
                resultsubMap(i-400)=T(i);
            end

                maxpos=vmax;
                minpos=0;
                minneg=-vmax;
                maxneg=0;
            resultcol=zeros(size(resultMap,1),3);
            if  maxpos==minpos
                if find(T==maxpos)<=400
                    flagpos=1;
                else
                    flagpos=2;
                end
                idx_pos=100;
                pos_step=1;
    %                         poscolor=(othercolor('YlOrRd9',100)); 
                resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

            else
                flagpos=3;
                pos_step=((maxpos)-(minpos))/100;
                idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                idx_pos(idx_pos<1)=1;
                idx_pos(idx_pos>100)=100;
    %                         poscolor=(othercolor('YlOrRd9',100));  
                resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
            end

           if  maxneg==minneg
                if find(T==maxneg)<=400
                    flagneg=1;
                else
                    flagneg=2;
                end
                idx_neg=1;
                neg_step=1;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

            else
                flagneg=3;
                neg_step=((-minneg)-(-maxneg))/100;
                idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                idx_neg(idx_neg<1)=1;
                idx_neg(idx_neg>100)=100;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
           end
            resultcol(resultMap==0,:)=repmat([0.75 0.75 0.75],sum(resultMap==0),1);


            S = convert_surface(surface);
            data=resultcol;
            D={};
            for ii = 1:2
            D{1,ii} = data(1:size(S{1}.coord,2),:);
            if numel(S) == 2
                D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
            end
            end
            jj=1;
                for ii = 1:numel(S)*2
                    idx = ceil(ii/2);
                    h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-(1/19)*j-0.02 .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                    h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                ,'FaceColor', 'interp','EdgeColor', 'none');

                    material dull; lighting phong;
                end
                set(h.axes(:,jj)                    , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]      )
            set(h.axes(1,:),'View',[-90 0]);
            set(h.axes(2,:),'View',[90 0]);
            set(h.axes(3,:),'View',[-90 0]);
            set(h.axes(4,:),'View',[90 0]);


            for ii = 1:4
                axes(h.axes(ii));
                h.camlight(ii) = camlight();
            end

           idx_subcluster = find(resultsubMap ~= 0);
            idx_subcluster_out=setdiff([1:16],idx_subcluster);

            for ii=5:8
                for jj=1
                    h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-(1/19)*j-0.02 .1 .1]);

                    for i=1:length(idx_subcluster)
                        thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                        temp1=thala.dat;
                        temp1(temp1~=idx_subcluster(i))=0;
                        thala.dat=temp1;
                        thalaregion=region(thala);
                        if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                            subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                            subidx_pos(subidx_pos<1)=1;
                            subidx_pos(subidx_pos>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                            subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                            subidx_neg(subidx_neg<1)=1;
                            subidx_neg(subidx_neg>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                        end


                    end
                 for i=1:length(idx_subcluster_out)
                    thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                    temp=thala.dat;
                    temp(temp~=idx_subcluster_out(i))=0;
                    thala.dat=temp;
                    thalaregion=region(thala);
                    p = imageCluster('cluster',region2struct(thalaregion),'color',[0.75,0.75,0.75],'alpha',1);
                end

                  material dull; lighting phong;
                set(h.axes(:,jj)               , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  ); 
                end
            end
            set(h.axes(5,1),'View',[-90 0]);
            set(h.axes(6,1),'View',[90 0]);     
            set(h.axes(7,1),'View',[45 0]);    
            set(h.axes(8,1),'View',[-45 0]);  
            material dull; lighting phong;
                set(gca                 , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  );
            set(gca,'color','w');
            set(gcf,'color','w');



    end
    
close all;
h=[];
h.figure = figure('Color','white','Units','normalized','Position',[0 0 0.5,1]);
allw=selection_freq_all;
    for j=1:8

           T=allw(:,j);
  
            T(isnan(T))=0;
            T_pos=T;
            T_pos(T<=0)=0;
            T_neg=T;
            T_neg(T>=0)=0;

            resultMap=zeros(size(brainmask,1),1);
            for i=1:400
                resultMap(brainmask==i)=T(i);
            end
            resultsubMap=zeros(16,1);
            for i=401:416
                resultsubMap(i-400)=T(i);
            end

                maxpos=vmax;
                minpos=0;
                minneg=-vmax;
                maxneg=0;
            resultcol=zeros(size(resultMap,1),3);
            if  maxpos==minpos
                if find(T==maxpos)<=400
                    flagpos=1;
                else
                    flagpos=2;
                end
                idx_pos=100;
                pos_step=1;
    %                         poscolor=(othercolor('YlOrRd9',100)); 
                resultcol(resultMap>0,:)=repmat(poscolor(round(idx_pos),:),sum(resultMap>0),1);

            else
                flagpos=3;
                pos_step=((maxpos)-(minpos))/100;
                idx_pos=(maxpos-resultMap(resultMap>0))./pos_step;
                idx_pos(idx_pos<1)=1;
                idx_pos(idx_pos>100)=100;
    %                         poscolor=(othercolor('YlOrRd9',100));  
                resultcol(resultMap>0,:)=poscolor(round(idx_pos),:);
            end

           if  maxneg==minneg
                if find(T==maxneg)<=400
                    flagneg=1;
                else
                    flagneg=2;
                end
                idx_neg=1;
                neg_step=1;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=repmat(negcolor(round(idx_neg),:),sum(resultMap<0),1);

            else
                flagneg=3;
                neg_step=((-minneg)-(-maxneg))/100;
                idx_neg=(-minneg-(-resultMap(resultMap<0)))./neg_step;
                idx_neg(idx_neg<1)=1;
                idx_neg(idx_neg>100)=100;
                 %negcolor=(othercolor('YlGnBu9',100));
                resultcol(resultMap<0,:)=negcolor(round(idx_neg),:);
           end
            resultcol(resultMap==0,:)=repmat([0.75 0.75 0.75],sum(resultMap==0),1);


            S = convert_surface(surface);
            data=resultcol;
            D={};
            for ii = 1:2
            D{1,ii} = data(1:size(S{1}.coord,2),:);
            if numel(S) == 2
                D{2,ii} = data(size(S{1}.coord,2)+1:end,:);
            end
            end
            jj=1;
                for ii = 1:numel(S)*2
                    idx = ceil(ii/2);
                    h.axes(ii,jj) = axes('Position',[-.09+ii*.1 1-(1/9)*j-0.02 .1 .1]); %[-.1+ii*.133 .9-.2*jj .1 .1]
                    h.patchs(ii,jj) = patch('Vertices',S{idx}.coord','Faces',S{idx}.tri,'FaceVertexCData',D{idx,jj}...
                ,'FaceColor', 'interp','EdgeColor', 'none');

                    material dull; lighting phong;
                end
                set(h.axes(:,jj)                    , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]      )
            set(h.axes(1,:),'View',[-90 0]);
            set(h.axes(2,:),'View',[90 0]);
            set(h.axes(3,:),'View',[-90 0]);
            set(h.axes(4,:),'View',[90 0]);


            for ii = 1:4
                axes(h.axes(ii));
                h.camlight(ii) = camlight();
            end

           idx_subcluster = find(resultsubMap ~= 0);
            idx_subcluster_out=setdiff([1:16],idx_subcluster);

            for ii=5:8
                for jj=1
                    h.axes(ii,jj) = axes('Position',[-.09+(ii)*.1 1-(1/9)*j-0.02 .1 .1]);

                    for i=1:length(idx_subcluster)
                        thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                        temp1=thala.dat;
                        temp1(temp1~=idx_subcluster(i))=0;
                        thala.dat=temp1;
                        thalaregion=region(thala);
                        if resultsubMap(idx_subcluster(i))>0 && flagpos==3                                   

                            subidx_pos=(maxpos-resultsubMap(idx_subcluster(i)))./pos_step;
                            subidx_pos(subidx_pos<1)=1;
                            subidx_pos(subidx_pos>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(round(subidx_pos),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))>0 && flagpos==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',poscolor(100,:),'alpha',1);

                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==3

                            subidx_neg=(-minneg-(-resultsubMap(idx_subcluster(i))))./neg_step;
                            subidx_neg(subidx_neg<1)=1;
                            subidx_neg(subidx_neg>100)=100;
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(round(subidx_neg),:),'alpha',1);
                        elseif resultsubMap(idx_subcluster(i))<0 && flagneg==2
                            p = imageCluster('cluster',region2struct(thalaregion),'color',negcolor(1,:),'alpha',1);

                        end


                    end
                 for i=1:length(idx_subcluster_out)
                    thala=fmri_data('D:\subcortex-master\Group-Parcellation\3T\Subcortex-Only\Tian_Subcortex_S1_3T.nii');
                    temp=thala.dat;
                    temp(temp~=idx_subcluster_out(i))=0;
                    thala.dat=temp;
                    thalaregion=region(thala);
                    p = imageCluster('cluster',region2struct(thalaregion),'color',[0.75,0.75,0.75],'alpha',1);
                end

                  material dull; lighting phong;
                set(h.axes(:,jj)               , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  ); 
                end
            end
            set(h.axes(5,1),'View',[-90 0]);
            set(h.axes(6,1),'View',[90 0]);     
            set(h.axes(7,1),'View',[45 0]);    
            set(h.axes(8,1),'View',[-45 0]);  
            material dull; lighting phong;
                set(gca                 , ...
                'Visible'           , 'off'         , ...
                'DataAspectRatio'   , [1 1 1]       , ...
                'PlotBoxAspectRatio', [1 1 1]  );
            set(gca,'color','w');
            set(gcf,'color','w');



    end   
    
    
    
    
    %% one sample
    
    
T2map_1=[];
T2map_2=[];
T2map_3=[];
age_gt= 30
age_lt= 100
% Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & counts'==1,:,:),...
%     S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt  & countsNC'==1,:,:));
% 
% RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & counts'==1,:);...
%     Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt & countsNC'==1,:)];
% RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt & counts'==1)),1);...
%     -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt & countsNC'==1)),1)];
% 
Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt ,:,:));

RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt ,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:)];
RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt)),1);...
    -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt )),1)];
RegressorsTemp(:,[2,3,4,5])=RegressorsTemp(:,[2,3,4,5])-repmat(mean(RegressorsTemp(:,[2,3,4,5]),1),[length(RegressorsTemp(:,[2,3,4,5])),1]);

T2map_1=[];
T2map_2=[];
T2map_3=[];
for net=1:18
    
    network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
   
      for i=1:416
        DependentVariable=network1(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
       RegressorsTemp,[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(RegressorsTemp,1)-...
        size(RegressorsTemp,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    
end





T2map_1=[];
T2map_2=[];
T2map_3=[];
Cov=RegressorsTemp(RegressorsTemp(:,1)==-1,[2,3,4,5]);
Cov= Cov- repmat(mean(Cov,1),[length(Cov),1]);
Cov=[ones(length(Cov),1),Cov];
for net=1:18
    
    network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
   
      for i=1:416
        DependentVariable=network1(RegressorsTemp(:,1)==-1,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
       Cov,[1,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(RegressorsTemp,1)-...
        size(RegressorsTemp,2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    
end




new=S_mci_nc_temp(:,:,:);
new(isnan(new))=0;
T2map=[];
Regressors_allnewTemp2=Regressors_allnewTemp(:,[1:7]);
age_gt=50
age_lt=95
Cov=Regressors_allnewTemp2(Regressors_allnewTemp2(:,6)<0.5 & ...
    Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt,[1,2,3,4,6,7]);
Cov(:,[2,3,4,5])=Cov(:,[2,3,4,5])-repmat(mean(Cov(:,[2,3,4,5]),1),[length(Cov(:,[2,3,4,5])),1]);

for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp2(:,5)',[], 0);
    for i=1:416
        DependentVariable=data_harmonized(i,Regressors_allnewTemp2(:,6)<0.5 & Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Cov,Contrastnew(1:6),'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size(Regressors_allnewTemp2(Regressors_allnewTemp2(:,6)<0.5 & Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp2(Regressors_allnewTemp2(:,6)<0.5 & Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end

new=S_mci_nc_temp(:,:,:);
new(isnan(new))=0;
T2map=[];
Regressors_allnewTemp2=Regressors_allnewTemp(:,[1:7]);
age_gt=50
age_lt=95
Cov=Regressors_allnewTemp2(Regressors_allnewTemp2(:,1)==-1 & Regressors_allnewTemp2(:,6)<0.5 & Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt,[2,3,4,6]);

Cov= Cov- repmat(mean(Cov,1),[length(Cov),1]);
Cov=[ones(length(Cov),1),Cov];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp2(:,5)',[], 0);
    for i=1:416
        DependentVariable=data_harmonized(i,Regressors_allnewTemp2(:,1)==-1 & Regressors_allnewTemp2(:,6)<0.5 & Regressors_allnewTemp2(:,3)>=age_gt & Regressors_allnewTemp2(:,3)<=age_lt)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Cov,[1,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
%     df=size(Cov,1)-...
%         size(Cov,2)
%     PMap=2*(1-tcdf(abs(Ttemp), df));
%     qThreshold=0.05;%设置p值
%     SortP=sort(PMap);
%     V=length(SortP);
%     
%     I=(1:V)';
%     cVID = 1;
%     cVN  = sum(1./(1:V));
%     P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
%     Thresholded=zeros(size(PMap));
%     if ~isempty(P)
%     Thresholded(find(PMap<=P))=1;
%     end
%     T_fdr(:,net)=Ttemp.*Thresholded;

    
end




%%
Regressor_ELMCIsub1=Regressor_ELMCIsub1(:,[1,2,3,6,5,4,7]);
% Regressor_ELMCIsub1(:,1)=ones(size(Regressor_ELMCIsub1,1),1);
Regressor_ELMCIsub2new(:,[1,2,3,4,6,7])=Regressor_ELMCIsub2(:,[1,2,3,5,4,6]);
Regressor_ELMCIsub2new(:,5)=3*ones(size(Regressor_ELMCIsub2new,1),1);
Regressor_NCnew(:,[1,2,3,4,6,7])=Regressor_NC(:,[1,2,3,5,4,6]);
Regressor_NCnew(:,5)=3*ones(size(Regressor_NCnew,1),1);
ELMCI_MCindex=[Regressor_ELMCIsub1(:,1);Regressor_ELMCIsub2(:,1);-2*Regressors_allnewTempNC(:,1);2*Regressor_NCnew(:,1)];

Datatemp=cat(1,S_elmci_temp,S_ELMCI_temp,S_mci_nc_tempNC,S_nc_temp);

Regressorstemp=[Regressor_ELMCIsub1;Regressor_ELMCIsub2new;Regressors_allnewTempNC;Regressor_NCnew];

Regressorstemp(:,1)=[1*ones(size(S_elmci_temp,1),1);1*ones(size(S_ELMCI_temp,1),1);...
    -1*ones(size(S_mci_nc_tempNC,1),1);-1*ones(size(S_nc_temp,1),1)];

age_gt=65;
age_lt=80;
Cov=Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]);
Cov(:,[2,3,4,5])= Cov(:,[2,3,4,5])- repmat(mean(Cov(:,[2,3,4,5]),1),[length(Cov(:,[2,3,4,5])),1]);
T2map_1=[];
for net=1:18
network1=Datatemp(:,:,yeoindex17==net);

    network1=mean(network1,3);
    network1=squeeze(network1);
    
%     data_harmonized = combat(network1', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
      data_harmonized = combat(network1', Regressorstemp(:,5)',[], 0);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Cov,[1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    

    
end

network=squeeze(mean(Datatemp,2));
data_harmonized = combat(network', Regressorstemp(:,5)', Regressorstemp(:,[2,3,4,6]), 0);
%     data_harmonized = combat(network', Regressors_allnewTemp(:,5)',Regressors_allnewTemp(:,[2,3,4,6]), 0);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Cov,[1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;





    
new=Datatemp;
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end

data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
% data = myzscore(data,0,2);
% data = zScoreToSubset(data,cInd);
% data=reshape(data,size(data,1),18,18);

% for sub=1:size(data_mat_mean_all,1)
%     temp=squeeze(data_mat_mean_all(sub,:,:));
%     mu = nanmean(temp(:));
%     sigma = nanstd(temp(:),0);
%     z = bsxfun(@minus,temp(:), mu);
%     z = bsxfun(@rdivide, z, sigma);
%     data_mat_mean_all(sub,:,:)=reshape(z,18,18);
% end
data_harmonized = combat(data',  Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
data=reshape(data_harmonized',size(data_harmonized',1),18,18);
for i=1:18
for j=1:18
        DependentVariable=data(Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,i,j);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
        Cov,[1,0,0,0,0,0],'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;


age_gt=55;
age_lt=95;

T2map_1=[];
for net=1:18
network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
%     data_harmonized = combat(network1', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
      data_harmonized = combat(network1', Regressorstemp(:,5)',[], 0);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
            (ELMCI_MCindex==-1 | ELMCI_MCindex==2))';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==-1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
        [1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==-1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==-1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    

    
end


age_gt=55;
age_lt=95;

T2map_1=[];
index=1;
for net=1:18
network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);
%     data_harmonized = combat(network1', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
      data_harmonized = combat(network1', Regressorstemp(:,5)',[], 0);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
            (ELMCI_MCindex==index | ELMCI_MCindex==2))';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
        [1,0,0,0,0,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(Regressorstemp(Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index| ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    

    
end





age_gt=55;
age_lt=95;
index=1
network=squeeze(mean(Datatemp,2));
% data_harmonized = combat(network', Regressorstemp(:,5)',[], 0);
    data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);

          for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
            (ELMCI_MCindex==index | ELMCI_MCindex==2))';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,5,6,7]),...
        [1,0,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
   df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,5,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,5,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));

    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;
    
    
age_gt=55;
age_lt=95;

network=squeeze(mean(Datatemp,2));
% data_harmonized = combat(network', Regressorstemp(:,5)',[], 0);
    data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);

          for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
            (ELMCI_MCindex==index | ELMCI_MCindex==2))';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
        [1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
   df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));

    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;
       
    
    
    
    
    
    
    
age_gt=55;
age_lt=95;
index=1
network=squeeze(mean(Datatemp,2));
data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 1);
      for i=1:416
        DependentVariable=data_harmonized(i,Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
            (ELMCI_MCindex==index | ELMCI_MCindex==2))';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
        Regressorstemp(Regressorstemp(:,5)==3 &  Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
        [1,0,0,0,0,0],'T'); 
    end
        Ttemp=TF_ForContrast';
   df=size(Regressorstemp(Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,5)==3 & Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
%     P=0.01
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;


Datatemp=cat(1,S_elmci_temp,S_ELMCI_temp,S_mci_nc_tempNC);
Regressorstemp=[Regressor_ELMCIsub1;Regressor_ELMCIsub2new;Regressors_allnewTempNC];
Regressorstemp(:,1)=[1*ones(size(S_elmci_temp,1),1);1*ones(size(S_ELMCI_temp,1),1);...
    -1*ones(size(S_mci_nc_tempNC,1),1)];
ELMCI_MCindex=[Regressor_ELMCIsub1(:,1);Regressor_ELMCIsub2(:,1);-2*Regressors_allnewTempNC(:,1)];
age_gt=60;
age_lt=95;
index=1
network=squeeze(mean(Datatemp,3));
data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
for i=1:416
DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
    (ELMCI_MCindex==index | ELMCI_MCindex==2))';
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressorstemp( Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
[1,0,0,0,0,0],'T'); 
end
Ttemp=TF_ForContrast';
df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
%     P=0.01
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=Ttemp.*Thresholded;

Datatemp=cat(1,S_ELMCI_temp,S_mci_nc_tempNC);
Regressorstemp=[Regressor_ELMCIsub2new;Regressors_allnewTempNC];
Regressorstemp(:,1)=[1*ones(size(S_ELMCI_temp,1),1);...
    -1*ones(size(S_mci_nc_tempNC,1),1)];
ELMCI_MCindex=[Regressor_ELMCIsub2(:,1);-2*Regressors_allnewTempNC(:,1)];
age_gt=60;
age_lt=95;index=-1
network=squeeze(mean(Datatemp,2));
data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
for i=1:416
DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
    (ELMCI_MCindex==index | ELMCI_MCindex==2))';
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressorstemp( Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
[1,0,0,0,0,0],'T'); 
end
Ttemp=TF_ForContrast';
df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
%     P=0.01
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=Ttemp.*Thresholded;





Datatemp=cat(1,S_nc_temp,S_mci_nc_tempNC);
Regressorstemp=[Regressor_NCnew;Regressors_allnewTempNC];
Regressorstemp(:,1)=[1*ones(size(S_nc_temp,1),1);...
    -1*ones(size(S_mci_nc_tempNC,1),1)];
ELMCI_MCindex=[Regressor_NCnew(:,1);-2*Regressors_allnewTempNC(:,1)];
age_gt=60;
age_lt=95;
index=1
network=squeeze(mean(Datatemp,2));
data_harmonized = combat(network', Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
for i=1:416
DependentVariable=data_harmonized(i,Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt &...
    (ELMCI_MCindex==index | ELMCI_MCindex==2))';
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressorstemp( Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),...
[1,0,0,0,0,0],'T'); 
end
Ttemp=TF_ForContrast';
df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
(ELMCI_MCindex==index | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
%     P=0.01
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=Ttemp.*Thresholded;






data_harmonized = combat(data',  Regressorstemp(:,5)',Regressorstemp(:,[2,3,4,6]), 0);
data=reshape(data_harmonized',size(data_harmonized',1),18,18);
for i=1:18
for j=1:18
        DependentVariable=data(Regressorstemp(:,6)<0.5  & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==1 | ELMCI_MCindex==2),i,j);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
         Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
    df=size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt& ...
        (ELMCI_MCindex==1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),1)-...
        size(Regressorstemp(Regressorstemp(:,6)<0.5 & Regressorstemp(:,3)>=age_gt & Regressorstemp(:,3)<=age_lt&...
        (ELMCI_MCindex==1 | ELMCI_MCindex==2),[1,2,3,4,6,7]),2)
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;



%%
option.method = 3;
option.num_boot = 5000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata =X;%term

exp{1} = Y;%gene
result = pls_analysis(exp,416, 1, option); 
lv=1;
[B,I] = sort(result.boot_result.orig_corr(:,lv),'descend');  % sort loadings
npos = length(find(result.boot_result.orig_corr(:,lv) > 0)); % number of positive loadings
nneg = length(find(result.boot_result.orig_corr(:,lv) < 0)); % number of negative loadings
tpos = I(1:floor(0.25*npos));                                % get top 25% of positive loadings
tneg = I(end-floor(0.25*nneg):end);                          % get top 25% of negative loadings
pos_terms = t(tpos);                                         % these are the positive terms contributing most
neg_terms = t(tneg);   
load('H:\TSA_testCamcan\FC\covariance\groupnew.mat')
roiindex=[];
yeoindex=[groupnew;ones(16,1)*8];
yeoindex_new=[];
for i=1:8
[ind1,ind2]=find(yeoindex==i);
roiindex=[roiindex;ind1];
yeoindex_new=[yeoindex_new;ones(size(ind1,1),1)*i];
end

map = [
    219 2 10;
    231 95 27;
    238 146 43;
    246 191 65;
    246 236 84;
    202 222 169;
    147 205 137;
    100 190 120;   % ← 插入的新颜色
    76 177 99
]/ 255;

index_sequence = [1:9];     % 哪几行
repeat_list = [2,2,2,2,2,3,3,1,1];        % 每行分别重复几次

% 构建 colormap
resultmap = [];
for i = 1:length(index_sequence)
    resultmap = [resultmap; repmat(map(index_sequence(i),:), repeat_list(i), 1)];
end


rsn=yeoindex17;
figure;               % distribute scores (as boxplots) in 7 networks
subplot(2,1,1)
boxplot(-1*result.usc(yeo_17labels_7_to_17,1),rsn,'Colors',resultmap,...
    'BoxStyle','filled','Whisker',0,'OutlierSize',2.5,'Symbol','.','MedianStyle','target','Widths',1) % temp scores

h = findobj(gca, 'Tag', 'Box');
set(h,'LineWidth',8);
% 
% h = findobj(gca, 'Tag', 'Median');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h = findobj(gca, 'Tag', 'Lower Adjacent Value');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% h = findobj(gca, 'Tag', 'Upper Adjacent Value') ;
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h = findobj(gca, 'Tag', 'Lower Whisker');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% h = findobj(gca, 'Tag', 'Upper Whisker');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h=findobj(gca,'tag','Outliers'); 
% set(h,'Marker','o','MarkerFaceColor', [0.5,0.5,0.5], 'MarkerSize',6,....
%     'Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5]); 

name={'VisC','VisP','SomMotA','SomMotB','DorsAttnA','DorsAttnB','SalVentAttnA','SalVentAttnB','LimbicA','LimbicB','ContA','ContB','ContC','DefaultA','DefaultB','DefaultC','TempPar','Subcortex'};
set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','bold')
set(gca,'XTickLabels',name);
set(gca,'XTickLabelRotation',35);
set(gca,'linewidth',1);
box off

rsn=yeoindex;            % distribute scores (as boxplots) in 7 networks
subplot(2,1,2)
boxplot(-1*result.usc(:,1),rsn,'Colors',map([1,2,3,4,5,6,7,9],:),'BoxStyle',...
    'filled','Whisker',0,'OutlierSize',2.5,'Symbol','.','MedianStyle','target') % gene scores

h = findobj(gca, 'Tag', 'Box');
set(h,'LineWidth',8);
% h = findobj(gca, 'Tag', 'Box');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h = findobj(gca, 'Tag', 'Median');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h = findobj(gca, 'Tag', 'Lower Adjacent Value');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% h = findobj(gca, 'Tag', 'Upper Adjacent Value') ;
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% 
% h = findobj(gca, 'Tag', 'Lower Whisker');
% set(h,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
% h = findobj(gca, 'Tag', 'Upper Whisker');
% set(h,'Color',[0.5,5,0.5],'LineWidth',1.5);
% 
% h=findobj(gca,'tag','Outliers'); 
% set(h,'Marker','o','MarkerFaceColor', [0.5,0.5,0.5], 'MarkerSize',6,....
%     'Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5]); 
name={'VN','SMN','DAN','VAN','Limbic','FPN','DMN','Subcortex'};
set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','bold')
set(gca,'XTickLabels',name);
set(gca,'XTickLabelRotation',35);
set(gca,'linewidth',1);
box off
set(gcf,'color','w');
%% camcan
files=dir([pwd,'/*.mat']);
for i=1:length(files)
    i
    load(files(i).name)
    T_all(i,:,:)=S(:,:,5);
    
end
S_mci_nc=T_all(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);


for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end

S_mci_nc=T_all;
for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    S_mci_nc_camcan_temp(i,:,:)=ZScoreMatrix(temp,1);
end


new =S_mci_nc_camcan_temp;
new(isnan(new))=0;
Regressorsnew=Regressors;
[a,b]=find(isnan(Regressors));
new(a,:,:)=[];
Regressorsnew(a,:)=[];


    network=mean(new,2);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(Regressorsnew(:,4)<0.5,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],[1,0,0,0,0],'T'); 
    end

        Ttemp=TF_ForContrast';
    df=size([Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],1)-...
        size([Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;


T2map=[];
for net=1:18
    network=new(:,:,yeoindex17==net);

    network=mean(network,3);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network(Regressorsnew(:,4)<0.5,i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],[1,0,0,0,0],'T'); 
    end
    T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size([Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],1)-...
        size([Regressorsnew(Regressorsnew(:,4)<0.5,:),ones(size(Regressorsnew(Regressorsnew(:,4)<0.5,:),1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end






young_index=find(Regressorsnew(:,1)>20 & Regressorsnew(:,1)<=40 & Regressorsnew(:,4)<=0.5);

middle_index=find(Regressorsnew(:,1)>40 & Regressorsnew(:,1)<=60 & Regressorsnew(:,4)<=0.5);

old_index=find(Regressorsnew(:,1)>60  & Regressorsnew(:,4)<=0.5);

Regressors_temp=[Regressorsnew(old_index,:);Regressorsnew(young_index,:)];
Regressors_temp=[Regressors_temp(:,[2,4]),zscore(Regressors_temp(:,3)),[ones(size(old_index,1),1);-1*ones(size(young_index,1),1)]];


Regressors_temp=[Regressorsnew(middle_index,:);Regressorsnew(young_index,:)];
Regressors_temp=[Regressors_temp(:,[2,4]),zscore(Regressors_temp(:,3)),[ones(size(middle_index,1),1);-1*ones(size(young_index,1),1)]];

T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network([middle_index;young_index],i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressors_temp,ones(size(Regressors_temp,1),1)],[0,0,0,1,0],'T'); 
    end
    T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size([Regressors_temp,ones(size(Regressors_temp,1),1)],1)-...
        size([Regressors_temp,ones(size(Regressors_temp,1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
[coeff, score1, latent, tsquared, explained, mu] = pca(T2map(:,yeo_17labels_17_to_7)');
plot_hemisphere_zhengnew(score1(:,1),{left_surface,right_surface},brainmask,1,poscolor,negcolor);


T2map=[];
for net=1:18
    network=new(:,:,yeoindex17==net);

    network=mean(network,3);
    network=squeeze(network);
    for i=1:416
        DependentVariable=network([middle_index;young_index],i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressors_temp,ones(size(Regressors_temp,1),1)],[0,0,0,1,0],'T'); 
    end
    T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size([Regressors_temp,ones(size(Regressors_temp,1),1)],1)-...
        size([Regressors_temp,ones(size(Regressors_temp,1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end
[coeff, score1, latent, tsquared, explained, mu] = pca(T2map(:,yeo_17labels_17_to_7)');
plot_hemisphere_zhengnew(score1(:,1),{left_surface,right_surface},brainmask,1,poscolor,negcolor);


new=enrichedFC;
new(isnan(new))=0;
Regressorsnew=Regressors;
[a,b]=find(isnan(Regressors));
new(a,:,:)=[]; 
Regressorsnew(a,:)=[];

T2map_enriched=[];
for net=1:19
   
    for i=1:416
        DependentVariable=enrichedFC([old_index;young_index],i,net);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressors_temp,ones(size(Regressors_temp,1),1)],[0,0,0,1,0],'T'); 
    end
    T2map_enriched=[T2map_enriched;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size([Regressors_temp,ones(size(Regressors_temp,1),1)],1)-...
        size([Regressors_temp,ones(size(Regressors_temp,1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_enriched(:,net)=Ttemp.*Thresholded;

    
end
new=S_mci_nc_temp;
new(isnan(new))=0;
Regressorsnew=Regressors;
[a,b]=find(isnan(Regressors));
new(a,:,:)=[];
Regressorsnew(a,:)=[];
network=squeeze(mean(new,2));
    for i=1:416
        DependentVariable=network([middle_index;young_index],i);
[b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   [Regressors_temp,ones(size(Regressors_temp,1),1)],[0,0,0,1,0],'T'); 
    end
%     T2map=[T2map;TF_ForContrast];
        Ttemp=TF_ForContrast';
    df=size([Regressors_temp,ones(size(Regressors_temp,1),1)],1)-...
        size([Regressors_temp,ones(size(Regressors_temp,1),1)],2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;
    
    
    
name='RdBu10';
cmp = load('colorData.mat',name);
cmp_pos=cmp.(name)(1:5,:);
cmp_neg=cmp.(name)(6:end,:);

cmp_pos = interp1(linspace(0,1,size(cmp_pos,1)),cmp_pos,linspace(0,1,100),'cubic');
cmp_pos(cmp_pos<0) = 0;
cmp_pos(cmp_pos>1) = 1;

cmp_neg = interp1(linspace(0,1,size(cmp_neg ,1)),cmp_neg ,linspace(0,1,100),'cubic');
cmp_neg(cmp_neg<0) = 0;
cmp_neg(cmp_neg>1) = 1;

poscolor=cmp_pos;

negcolor=flipud(cmp_neg);

















index=8;
MCI_data=S_mci_nc_temp(1:248,:,:);
NC_data=S_mci_nc_temp(249:end,:,:);
Regressors_allnewTemp_MCI=Regressors_allnewTemp(1:248,:);
Regressors_allnewTemp_NC=Regressors_allnewTemp(249:end,:);



new=cat(1,MCI_data(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5,:,:),...
    NC_data(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5,:,:));



Regressors_Temp=[Regressors_allnewTemp_MCI(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5,:);...
    Regressors_allnewTemp_NC(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5,:)];


clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data', Regressors_Temp(:,5)',[], 1);

data=reshape(data_harmonized',size(data_harmonized',1),18,18);
for i=1:18
for j=1:18
        DependentVariable=data(:,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_Temp(:,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_Temp(:,[1,2,3,4,6,7]),1)-...
    size(Regressors_Temp(:,[1,2,3,4,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;


tempdir=dir('I:\ADNI\MCI\Spreading_Smooth');
tempdir(1:2)=[];
tempdir(208)=[];
name_MCI=cellfun(@(folder) folder(1:10),{tempdir.name},'UniformOutput',false);


tempdir=dir('I:\ADNI\NC_all_Spreading_Smooth');
tempdir(1:2)=[];
name_nc=cellfun(@(folder) folder(1:10),{tempdir.name},'UniformOutput',false);

subname=cat(2,name_MCI,name_nc)';
subname=subname(Regressors_allnewTemp(:,6)<0.5);
group=Regressors_Temp(:,1);
age=Regressors_Temp(:,3);
sex=Regressors_Temp(:,2);
tiv=Regressors_Temp(:,4);
fd=Regressors_Temp(:,6);
times_new=times(Regressors_allnewTemp(:,6)<0.5);


for i=1:18
for j=1:18
        DependentVariable=data(:,i,j);
    T = table(subname, group, times_new, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Time', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group + Time + Age + Sex + Motion + Tiv +(1|Subject)');
TF_ForContrast_brain_mean3(i,j)=mdl.Coefficients.tStat(2);
p_ForContrast_brain_mean3(i,j)=mdl.Coefficients.pValue(2);
end
end

PMap=p_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;



new=cat(1,MCI_data(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5 ,:,:),...
    NC_data(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5 ,:,:));
Regressors_Temp=[Regressors_allnewTemp_MCI(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5 ,:);...
    Regressors_allnewTemp_NC(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5,:)];
subname=cat(2,name_MCI,name_nc)';
subname=subname(Regressors_allnewTemp(:,6)<0.5);
group=Regressors_Temp(:,1);
age=Regressors_Temp(:,3);
sex=Regressors_Temp(:,2);
tiv=Regressors_Temp(:,4);
fd=Regressors_Temp(:,6);
times_new=times(Regressors_allnewTemp(:,6)<0.5);
for i=1:416
    i
for j=1:416
    if i~=j
        DependentVariable=new(:,i,j);
    T = table(subname, group, times_new, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Time', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group + Time + Age + Sex + Motion + Tiv +(1|Subject)');
TF_ForContrast_brain_mean3(i,j)=mdl.Coefficients.tStat(2);
p_ForContrast_brain_mean3(i,j)=mdl.Coefficients.pValue(2);
    else
        TF_ForContrast_brain_mean3(i,j)=0;
        p_ForContrast_brain_mean3(i,j)=mdl.Coefficients.pValue(2);
    end
end
end


new=cat(1,MCI_data(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5,:,:),...
    NC_data(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5,:,:));
Regressors_Temp=[Regressors_allnewTemp_MCI(times_MCI'~=index & Regressors_allnewTemp_MCI(:,6)<0.5,:);...
    Regressors_allnewTemp_NC(times_NC'~=index & Regressors_allnewTemp_NC(:,6)<0.5,:)];
subname=cat(2,name_MCI,name_NC)';
subname=subname(Regressors_allnewTemp(:,6)<0.5  & Regressors_allnewTemp(:,5)==1);
group=Regressors_Temp(:,1);
age=Regressors_Temp(:,3);
sex=Regressors_Temp(:,2);
tiv=Regressors_Temp(:,4);
fd=Regressors_Temp(:,6);
times_new=times(Regressors_allnewTemp(:,6)<0.5);



network=squeeze(mean(new(:,yeo_17labels_17_to_7,yeo_17labels_17_to_7),3));%2是入

Z_scores=Z_scores(:,yeo_17labels_17_to_7);
network=Z_scores(Regressors_allnewTemp(:,6)<0.5,:);
network = combat(network', Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,5)',[], 0);
    for i=1:416
        DependentVariable=network(i,:)';
  T = table(subname, group, times_new, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Time', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group + Time + Age + Sex + Motion + Tiv +(1|Subject)');
TF_ForContrast(i)=mdl.Coefficients.tStat(2);
p_ForContrast(i)=mdl.Coefficients.pValue(2);
    end
        Ttemp=TF_ForContrast';

 
    PMap=p_ForContrast';
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr=Ttemp.*Thresholded;

S_mci_nc_tempNC=S_mci_nc_temp(Regressors_allnewTemp(:,1)==-1,:,:);




net=18
T_fdr_pos=T_fdr(:,net);
T_fdr_pos(T_fdr_pos<=0)=0;
T_fdr_neg=T_fdr(:,net);
T_fdr_neg(T_fdr_neg>=0)=0;

index=sum(logical(T_fdr_neg),2);
dataTemp=squeeze(Z_scoresMCI(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2))


Z_scoresNC=Z_scores(:,Regressors_allnewTemp(:,1)==1 &Regressors_allnewTemp(:,6)<0.5 ,:);
dataTemp=squeeze(Z_scoresNC(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))

Z_scoresNC=Z_scores(:,Regressors_allnewTemp(:,1)==-1 &Regressors_allnewTemp(:,6)<0.5 ,:);
dataTemp=squeeze(Z_scoresNC(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))



% dataTemp=squeeze(Z_scoresMCI2(net,:,yeo_17labels_17_to_7));
% dataTemp=dataTemp(:,find(index));
% mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==-1,:),2))
% mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==1,:),2))


dataTemp=squeeze(Z_scoresTMS(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(1:17,:),2))
mean(mean(dataTemp(18:end,:),2))



net=18
T_fdr_pos=T_fdr(:,net);
T_fdr_pos(T_fdr_pos<=0)=0;
T_fdr_neg=T_fdr(:,net);
T_fdr_neg(T_fdr_neg>=0)=0;

index=sum(logical(T_fdr_neg),2);
dataTemp=squeeze(abs(Z_scoresMCI(net,:,yeo_17labels_17_to_7)));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2))


Z_scoresMCIall=Z_scores(:,Regressors_allnewTemp(:,1)==1 &Regressors_allnewTemp(:,6)<0.5 ,:);
dataTemp=squeeze(abs(Z_scoresNC(net,:,yeo_17labels_17_to_7)));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))

Z_scoresNC=Z_scores(:,Regressors_allnewTemp(:,1)==-1 &Regressors_allnewTemp(:,6)<0.5 ,:);
dataTemp=squeeze(abs(Z_scoresNC(net,:,yeo_17labels_17_to_7)));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))



% dataTemp=squeeze(Z_scoresMCI2(net,:,yeo_17labels_17_to_7));
% dataTemp=dataTemp(:,find(index));
% mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==-1,:),2))
% mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==1,:),2))


dataTemp=squeeze(abs(Z_scoresTMS(net,:,yeo_17labels_17_to_7)));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(1:17,:),2))
mean(mean(dataTemp(18:end,:),2))



for i=1:18
    for j=1:416
        Z_scoresMCIall(i,:,j)=Z_scoresMCIall(i,:,j)-mean(Z_scoresNC(i,:,j));
     
    end
end
for i=1:18
    for j=1:416
        Z_scoresMCI(i,:,j)=Z_scoresMCI(i,:,j)-mean(Z_scoresNC(i,:,j));
     
    end
end
for i=1:18
    for j=1:416
        Z_scoresMCI2(i,:,j)=Z_scoresMCI2(i,:,j)-mean(Z_scoresNC(i,:,j));
     
    end
end

for i=1:18
    for j=1:416
        Z_scoresMCI2(i,:,j)=Z_scoresMCI2(i,:,j)-mean(Z_scoresNC2(i,:,j));
     
    end
end



for i=1:18
    for j=1:416
        Z_scoresTMS(i,:,j)=Z_scoresTMS(i,:,j)-mean(Z_scoresNC(i,:,j));
     
    end
end


for i=1:18
    for j=1:416
        Z_scoresNC(i,:,j)=Z_scoresNC(i,:,j)-mean(Z_scoresNC(i,:,j));
     
    end
end

for i=1:18
    for j=1:416
        Z_scoresNC2(i,:,j)=Z_scoresNC2(i,:,j)-mean(Z_scoresNC2(i,:,j));
     
    end
end





net=18
T_fdr_pos=T_fdr(:,net);
T_fdr_pos(T_fdr_pos<=0)=0;
T_fdr_neg=T_fdr(:,net);
T_fdr_neg(T_fdr_neg>=0)=0;

index=sum(logical(T_fdr_neg),2);
dataTemp=squeeze(Z_scoresMCI(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2))

dataTemp=squeeze(Z_scoresMCIall(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))

dataTemp=squeeze(Z_scoresNC(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))

dataTemp=squeeze(Z_scoresMCI2(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==1,:),2))

dataTemp=squeeze(Z_scoresTMS(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(1:17,:),2))
mean(mean(dataTemp(18:end,:),2))

T_fdr_neg=T_fdr;
T_fdr_neg(T_fdr_neg>0)=0;

T_fdr_pos=T_fdr;
T_fdr_pos(T_fdr_pos<0)=0;
sum(logical(T_fdr),2)


index=sum(logical(T_fdr_pos),2);
dataTemp=squeeze(Z_scoresMCI(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2))

dataTemp=squeeze(Z_scoresMCIall(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))

dataTemp=squeeze(Z_scoresNC(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))



net=18
T_fdr_pos=T_fdr(:,net);
T_fdr_pos(T_fdr_pos<=0)=0;
T_fdr_neg=T_fdr(:,net);
T_fdr_neg(T_fdr_neg>=0)=0;


index=sum(logical(T_fdr_pos),2);
dataTemp=squeeze(Z_scoresMCI(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2))

dataall{1}=mean(dataTemp(Regressor_ELMCIsub1(:,1)==-1,:),2);
dataall{2}=mean(dataTemp(Regressor_ELMCIsub1(:,1)==1,:),2);


dataTemp=squeeze(Z_scoresMCIall(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))
dataall{3}=mean(dataTemp(:,:),2);

dataTemp=squeeze(Z_scoresNC(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))
dataall{4}=mean(dataTemp(:,:),2);


dataTemp=squeeze(Z_scoresNC2(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(:,:),2))
dataall{5}=mean(dataTemp(:,:),2);


dataTemp=squeeze(Z_scoresMCI2(net,:,yeo_17labels_17_to_7));
dataTemp=dataTemp(:,find(index));
mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==-1,:),2))
mean(mean(dataTemp(Regressor_ELMCIsub2(:,1)==1,:),2))

dataall{6}=mean(dataTemp(Regressor_ELMCIsub2(:,1)==-1,:),2);
dataall{7}=mean(dataTemp(Regressor_ELMCIsub2(:,1)==1,:),2);









for i=1:size(mci_mean,2)
             
    sss=mci_mean(:,i);
%     mci_Qmed(i)=median(sss);
    mci_Qmed(i)=mean(sss);
    mci_Q25e(i)=mci_Qmed(i)-quantile(sss,.25);
    mci_Q75e(i)=quantile(sss,.75)-mci_Qmed(i);
end

for i=1:size(nc_mean,2)
             
    sss=nc_mean(:,i);
%     nc_Qmed(i)=median(sss);
     nc_Qmed(i)=mean(sss);
    nc_Q25e(i)=nc_Qmed(i)-quantile(sss,.25);
    nc_Q75e(i)=quantile(sss,.75)-nc_Qmed(i);
end
hold on;
hl=errorbar(0.9:size(mci_mean,2)-1+0.9,mci_Qmed,mci_Q25e,mci_Q75e,'LineStyle', 'None');
hr=errorbar(1.1:size(nc_mean,2)-1+1.1,nc_Qmed,nc_Q25e,nc_Q75e,'LineStyle', 'None');
hl.Marker='o';
hl.MarkerEdgeColor= [.2 .2 .2];
hl.MarkerFaceColor= [208,79,88]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hl.MarkerSize=10;  
hl.LineWidth=1;
hl.Color='k';

hr.Marker='o';
hr.MarkerEdgeColor= [.2 .2 .2];
hr.MarkerFaceColor= [23,112,159]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hr.MarkerSize=10;  
hr.LineWidth=1;
hr.Color='k';








for i=1:size(dataall,2)
             
    sss=dataall{i};
%     mci_Qmed(i)=median(sss);
    mci_Qmed(i)=mean(sss);
    mci_Q25e(i)=mci_Qmed(i)-quantile(sss,.25);
    mci_Q75e(i)=quantile(sss,.75)-mci_Qmed(i);
end
hold on;
hl=errorbar(2.5,mci_Qmed(3),mci_Q25e(3),mci_Q75e(3),'LineStyle', 'None');
hr=errorbar(.5,mci_Qmed(4),mci_Q25e(4),mci_Q75e(4),'LineStyle', 'None');
hl1=errorbar(1.5,mci_Qmed(1),mci_Q25e(1),mci_Q75e(1),'LineStyle', 'None');
hr1=errorbar(3.5,mci_Qmed(2),mci_Q25e(2),mci_Q75e(2),'LineStyle', 'None');
[23,112,159]

hl.Marker='o';
hl.MarkerEdgeColor= [.2 .2 .2];
hl.MarkerFaceColor= [23,112,159]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hl.MarkerSize=10;  
hl.LineWidth=2;
hl.Color='k';

hr.Marker='o';
hr.MarkerEdgeColor= [.2 .2 .2];
hr.MarkerFaceColor= [208,79,88]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hr.MarkerSize=10;  
hr.LineWidth=2;
hr.Color='k';




hl1.Marker='o';
hl1.MarkerEdgeColor= [.2 .2 .2];
hl1.MarkerFaceColor= [ 35, 146, 199]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hl1.MarkerSize=10;  
hl1.LineWidth=2;
hl1.Color='k';

hr1.Marker='o';
hr1.MarkerEdgeColor= [.2 .2 .2];
hr1.MarkerFaceColor= [18, 92, 132 ]./255;%light blue (default) %dark blue  [0 0.4470 0.7410] (default chosen as first color)
hr1.MarkerSize=10;  
hr1.LineWidth=2;
hr1.Color='k';




xlim([0,5])





for net=1:18
Z_score=squeeze(Z_scores(net,:,yeo_17labels_17_to_7));
% Z_score=Z_score-repmat(squeeze(mean(Z_scores(net,Regressors_allnewTemp(:,1)==-1,yeo_17labels_17_to_7),2))',489,1);
network=Z_score(Regressors_allnewTemp(:,6)<0.5,:);
network = combat(network', Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,5)',[], 0);
for i=1:416
DependentVariable=network(i,:)';
[b,r,SSE_1,SSR_1, T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T');
end
T=TF_ForContrast';
T(isnan(T))=0;
mask=T;
mask(mask~=0)=1;
Ttemp=T(mask~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)...
-size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
Thresholdedmask=zeros(1,416);
Thresholdedmask(mask~=0)=Thresholded;
T_fdr(:,net)=T'.*Thresholdedmask;
end



%%
load('I:\ADNI\T_allnew.mat')
S_mci_nc=squeeze(T_all(:,:,:,5));
clear T_all


load('I:\ADNI\S_allnew.mat')
S_mci_nc=squeeze(S_all(:,:,:,5));
clear S_all

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==22)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1)=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end


for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end










%%
ELMCIexcel=xlsread('J:\zhengzihao_ELMCI_20241120\zhengzihao_ELMCI_20241120_3_20_2025.xlsx','B2:J2');

for i=1:length(ELMCIexcel)
    s1=ELMCIexcel{i,1};
    s2=ELMCIexcel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
%     if  str2num(s2{2})<10;
%         s2{2}=['0',s2{2}];
%     end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    ELMCIexcel{i,1}=[s1,'_',s2new];
    ELMCI_COV{i,1}=[s1,'_',s2new];
    if  strcmp(ELMCIexcel{i,2},'EMCI')
        ELMCI_COV{i,2}=-1;
    else
        ELMCI_COV{i,2}=1;
    end
    
    if ELMCIexcel{i,3}=='M'
        ELMCI_COV{i,3}=1;
    else
        ELMCI_COV{i,3}=2;
    end
    ELMCI_COV{i,4}=ELMCIexcel{i,4};
    
    
    
end
path='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);

[a,b1]=ismember(unique(fmri_name),ELMCI_COV(:,1));
b1(b1==0)=[];
Regressor_ELMCIsub1=ELMCI_COV([b1],:);
MCI_head_path='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_nii\preprocess\headmotion';
MCI_head_files=dir(MCI_head_path);
MCI_head_files(1:2)=[];
for i=1:length(MCI_head_files)

    cd([MCI_head_path,'\',MCI_head_files(i).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

    Regressor_ELMCIsub1{i,5}=mean(FD);
    
    cd(['J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open\',MCI_head_files(i).name]);
    
     files=dir('*.dcm');
    info=dicominfo([files(3).folder,'\',files(3).name]);    
    
    if strcmp(info.Manufacturer,'SIEMENS')
         Regressor_ELMCIsub1{i,6}=1;
    else
         Regressor_ELMCIsub1{i,6}=2;
    end

end
% T1excel1='J:\zhengzihao_ELMCI_20241120\fmri_Eyes_Open_T1_nii\Tissue_Volumes.csv';
T1_name=cellfun(@(name) strsplit(name,'/'),T1excel1,'UniformOutput',false);
T1_name=cellfun(@(name) char(name(10)),T1_name,'UniformOutput',false);
T1_name=cellfun(@(name) name(1:end-11),T1_name,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name);
b1(b1==0)=[];
fmri_name=fmri_name(a);
subname1=cellfun(@(folder) folder(1:10),fmri_name,'UniformOutput',false)';
T1volume1= T1volume1(b1);

[a,b2]=ismember(fmri_name,Regressor_ELMCIsub1(:,1));
Regressor_ELMCIsub1=Regressor_ELMCIsub1(b2,:);


files=dir('I:\ADNI\ELMCI\fmri_Eyes_Open_nii\Spreading_smooth_new\*.mat');
files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);
b4(b4==0)=[];
for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_elmci(i,:,:)=squeeze(T(:,:,7));
    
end


for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_elmci(i,:,:)=squeeze(S(:,:,5));
    
end
 Regressor_ELMCIsub1=cell2mat(Regressor_ELMCIsub1(:,[2:6]));
 Regressor_ELMCIsub1=[Regressor_ELMCIsub1,T1volume1,ones(135,1)];
 
 
 
 

 

 
 
 

 
 

 
 
 

%%
path='J:\zhengzihao_ELMCI_20241120\resting_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);

% T1excel2='J:\zhengzihao_ELMCI_20241120\resting_nii_T1_nii\Tissue_Volumes.csv';
T1_name=cellfun(@(name) strsplit(name,'/'),T1excel2,'UniformOutput',false);
T1_name=cellfun(@(name) char(name(10)),T1_name,'UniformOutput',false);
T1_name=cellfun(@(name) name(1:end-11),T1_name,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name);
b1(b1==0)=[];
fmri_name=fmri_name(a);
subname2=cellfun(@(folder) folder(1:10),fmri_name,'UniformOutput',false)';
T1volume2= T1volume2(b1);
[a,b2]=ismember(unique(fmri_name),ELMCI_COV(:,1));
Regressor_ELMCIsub2=ELMCI_COV([b2],:);

MCI_head_path='J:\zhengzihao_ELMCI_20241120\resting_nii\preprocess\headmotion';
MCI_head_files=dir(MCI_head_path);
MCI_head_files(1:2)=[];
MCI_head_name=cellfun(@(name) name(1:end-11),{MCI_head_files.name},'UniformOutput',false);
[a,b3]=ismember(fmri_name,MCI_head_name);
for i=1:length(fmri_name)

    cd([MCI_head_path,'\',MCI_head_files(b3(i)).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

    Regressor_ELMCIsub2{i,5}=mean(FD);

end
Regressor_ELMCIsub2=cell2mat(Regressor_ELMCIsub2(:,[2:5]));
Regressor_ELMCIsub2=[Regressor_ELMCIsub2,T1volume2,ones(286,1)];


files=dir('I:\ADNI\ELMCI\resting_nii\Spreading_smooth_new\*.mat');

files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);

for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_ELMCI(i,:,:)=squeeze(T(:,:,7));
    
end


for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_ELMCI(i,:,:)=squeeze(S(:,:,5));
    
end

files=dir('J:\zhengzihao_ELMCI_20241120\resting_nii\preprocess\BOLD_smooth\*.mat');
files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);
b4(b4==0)=[];
for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    FC_elmci2(i,:,:)=atanh(corr(theROITimeCoursesTotal));
    
end




%%
% NCexcel=xlsread('J:\zhengzihao_ELMCI_20241120\NC\zhengzihao_NC_resting_5_16_2025.csv','B2:J2');
for i=1:length(NCexcel)
    s1=NCexcel{i,1};
    s2=NCexcel{i,9};
    s2=s2(2:end);
    s2=strsplit(s2,'-');
    if str2num(s2{1})<10;
        s2{1}=['0',s2{1}];
    end
    
%     if  str2num(s2{2})<10;
%         s2{2}=['0',s2{2}];
%     end
    s2new=[s2{3},'-',s2{1},'-',s2{2}];
    NCexcel{i,1}=[s1,'_',s2new];
    NCexcel_COV{i,1}=[s1,'_',s2new];
    if  strcmp(NCexcel{i,2},'EMCI')
        NCexcel_COV{i,2}=-1;
    else
        NCexcel_COV{i,2}=1;
    end
    
    if NCexcel{i,3}=='M'
        NCexcel_COV{i,3}=1;
    else
        NCexcel_COV{i,3}=2;
    end
    NCexcel_COV{i,4}=NCexcel{i,4};
    
    
    
end

path='J:\zhengzihao_ELMCI_20241120\NC\resting_nii\preprocess\smooth_regress';
files=dir(path);
files(1:2)=[];
fmri_name=cellfun(@(name) name(1:end-11),{files.name},'UniformOutput',false);


% T1excel_NC='J:\zhengzihao_ELMCI_20241120\NC\Tissue_Volumes.csv';

T1_name_NC=cellfun(@(name) strsplit(name,'/'),T1excel_NC,'UniformOutput',false);
T1_name_NC=cellfun(@(name) char(name(10)),T1_name_NC,'UniformOutput',false);
T1_name_NC=cellfun(@(name) name(1:end-11),T1_name_NC,'UniformOutput',false);
[a,b1]=ismember(unique(fmri_name),T1_name_NC);
b1(b1==0)=[];
fmri_name=fmri_name(a);
subname3=cellfun(@(folder) folder(1:10),fmri_name,'UniformOutput',false)';
T1volume_NC= T1volume_NC(b1);

[a,b2]=ismember(unique(fmri_name),NCexcel_COV(:,1));
Regressor_NC=NCexcel_COV([b2],:);

NC_head_path='J:\zhengzihao_ELMCI_20241120\NC\resting_nii\preprocess\headmotion';
NC_head_files=dir(NC_head_path);
NC_head_files(1:2)=[];
NC_head_name=cellfun(@(name) name(1:end-11),{NC_head_files.name},'UniformOutput',false);
[a,b3]=ismember(fmri_name,NC_head_name);
for i=1:length(fmri_name)

    cd([NC_head_path,'\',NC_head_files(b3(i)).name]);
    head_file=dir(['rp*']);
    rp1=importdata([head_file.name]);

    rp = rp1;
    rp(:,4:6) = rp1(:,4:6)*180/pi;
    rpMax = max(abs(rp));
    rpDiff = diff(rp1);
    rpDiff(:,4:6) = rpDiff(:,4:6)*50;
    FD = [0;sum(abs(rpDiff),2)];

   Regressor_NC{i,5}=mean(FD);

end

 Regressor_NC=cell2mat(Regressor_NC(:,[2:5]));
 Regressor_NC=[ Regressor_NC,T1volume_NC,ones(152,1)];
 
files=dir('J:\zhengzihao_ELMCI_20241120\NC\Spreading_smooth_new\*.mat');
files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);
for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_nc(i,:,:)=squeeze(T(:,:,7));
    
end
for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    S_nc(i,:,:)=squeeze(S(:,:,5));
    
end

files=dir('J:\zhengzihao_ELMCI_20241120\NC\resting_nii\preprocess\BOLD_smooth\*.mat');
files_name=cellfun(@(name) name(1:end-32),{files.name},'UniformOutput',false);
[a,b4]=ismember(fmri_name,files_name);

for i=1:length(b4)
    load([files(i).folder,'\',files(b4(i)).name])
    FC_NC2(i,:,:)=atanh(corr(theROITimeCoursesTotal));
    
end




for i=1:size(S_elmci,1)
    temp=squeeze(S_elmci(i,:,:));
    temp(temp==1 | temp==max(S_elmci(:)))=NaN;
    S_elmci_temp(i,:,:)=ZScoreMatrix(temp);
end

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    temp(temp==1 | temp==max(S_mci_nc(:)))=NaN;
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp,1);
end


for i=1:size(S_ELMCI,1)
    temp=squeeze(S_ELMCI(i,:,:));
    temp(temp==1 | temp==max(S_ELMCI(:)))=NaN;
    S_ELMCI_temp(i,:,:)=ZScoreMatrix(temp);
end
for i=1:size(S_nc,1)
    temp=squeeze(S_nc(i,:,:));
    temp(temp==1 | temp==max(S_nc(:)))=NaN;
    S_nc_temp(i,:,:)=ZScoreMatrix(temp);
end




for i=1:size(S_elmci,1)
    temp=squeeze(S_elmci(i,:,:));
    S_elmci_temp(i,:,:)=ZScoreMatrix(temp);
end

for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp,1);
end


for i=1:size(S_ELMCI,1)
    temp=squeeze(S_ELMCI(i,:,:));
    S_ELMCI_temp(i,:,:)=ZScoreMatrix(temp);
end
for i=1:size(S_nc,1)
    temp=squeeze(S_nc(i,:,:));
    S_nc_temp(i,:,:)=ZScoreMatrix(temp);
end




new=S_mci_nc_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp(:,5)',Regressors_allnewTemp(:,[2,3,4,6]), 0);
    for i=1:416
        DependentVariable=data_harmonized(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end








T2map_1=[];
T2map_2=[];
T2map_3=[];
age_gt= 20
age_lt= 100
Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,1)==-1 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:,:));

RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,1)==-1 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:)];
RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,1)==-1 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt)),1);...
    -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt)),1)];


Datatemp=cat(1,S_ELMCI_temp(Regressor_ELMCIsub2(:,1)~=2 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:,:),...
    S_nc_temp(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:,:));

RegressorsTemp=[Regressor_ELMCIsub2(Regressor_ELMCIsub2(:,1)~=2 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt,:);...
    Regressor_NC(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt,:)];
RegressorsTemp(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,1)~=2 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,3)>=age_gt & Regressor_ELMCIsub2(:,3)<=age_lt)),1);...
    -1*ones(length(find(Regressor_NC(:,4)<0.5 & Regressor_NC(:,3)>=age_gt & Regressor_NC(:,3)<=age_lt)),1)];


T2map_1=[];
T2map_2=[];
T2map_3=[];
Datatemp=Datatemp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for net=1:18
    
    
    
    
    network1=Datatemp(:,yeoindex17==net,:);
    network1=mean(network1,2);
    network1=squeeze(network1);

      for i=1:416
        DependentVariable=network1(:,i);
        [b,r,SSE_1(i),SSR_1(i), T, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
       RegressorsTemp(:,[1,6]),[1,0],'T'); 
    end
     Ttemp=TF_ForContrast';
        T2map_1=[T2map_1;TF_ForContrast];
    df=size(RegressorsTemp(:,[1,6]),1)-...
        size(RegressorsTemp(:,[1,6]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr_1(:,net)=Ttemp.*Thresholded;
    
    
end


plot_hemisphere_zhengnew2new(squeeze(mean(mean(S_mci_nc_temp(Regressors_allnewTemp(:,1)==1,yeoindex==8,:),2),1)),{left_surface,right_surface},brainmask,1,poscolor,negcolor,1);
plot_hemisphere_zhengnew2new(squeeze(mean(mean(S_nc_temp(:,yeoindex==8,:),2),1)),{left_surface,right_surface},brainmask,1,poscolor,negcolor,1);
net=18;
corr(squeeze(mean(Z_scoresNC(net,:,yeo_17labels_17_to_7),2)),squeeze(mean(Z_scoresNC2(net,:,yeo_17labels_17_to_7),2)))







%%
Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,3)>65 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1) ,:);...
Regressor_NC(Regressor_NC(:,4)<0.5,:)];
Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,3)>65 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1)),1);...
-1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
new=cat(1,FC_elmci2((Regressor_ELMCIsub2(:,3)>65 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1),:,:),FC_NC2(Regressor_NC(:,4)<0.5,:,:));
subnamenew=cat(1,subname2((Regressor_ELMCIsub2(:,3)>65 & Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)==1)),subname3(Regressor_NC(:,4)<0.5));
group=Regressor(:,1);
age=Regressor(:,3);
sex=Regressor(:,2);
tiv=Regressor(:,5);
fd=Regressor(:,4);
T_all1=zeros(416,416);
for i=1:416
    for j=i+1:416
     DependentVariable=new(:,i,j);
    [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressor,[1,0,0,0,0,0],'T'); 
    T_all1(i,j)=TF_ForContrast;
    T_all1(j,i)=TF_ForContrast;
      
    end
end

T_all1=zeros(416,416);
for i=1:416
    i
    for j=i+1:416
     DependentVariable=new(:,i,j);
    T = table(subnamenew, group, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group  + Age + Sex + Motion + Tiv +(1|Subject)');
T_all1(i,j)=mdl.Coefficients.tStat(2);
p_all1(i,j)=mdl.Coefficients.pValue(2);
T_all1(j,i)=mdl.Coefficients.tStat(2);
p_all1(j,i)=mdl.Coefficients.pValue(2);
      
    end
end









Regressor_ELMCIsub1=Regressor_ELMCIsub1(:,[1,2,3,6,5,4,7]);
Regressor_ELMCIsub1(:,1)=ones(size(Regressor_ELMCIsub1,1),1);

new=cat(1,FC_elmci,fc_all(Regressors_allnewTemp(:,1)==-1,:,:));
T_all2=zeros(416,416);
Regressors_allnewTempNC=Regressors_allnewTemp(Regressors_allnewTemp(:,1)==-1,:);
Regressorstemp=[Regressor_ELMCIsub1;Regressors_allnewTempNC];

subnamenew=cat(1,subname1(( Regressor_ELMCIsub1(:,6)<0.5 )),...
    subname(Regressors_allnewTemp(:,1)==-1&Regressors_allnewTemp(:,6)<0.5));
group=Regressorstemp(Regressorstemp(:,6)<0.5,1);
age=Regressorstemp(Regressorstemp(:,6)<0.5,3);
sex=Regressorstemp(Regressorstemp(:,6)<0.5,2);
tiv=Regressorstemp(Regressorstemp(:,6)<0.5,4);
fd=Regressorstemp(Regressorstemp(:,6)<0.5,6);

for i=1:416
    for j=i+1:416
     DependentVariable=new(:,i,j);
%      data_harmonized = combat(DependentVariable', Regressorstemp(:,5)',[], 1);
    [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable(Regressorstemp(:,6)<0.5),...
   Regressorstemp(Regressorstemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    T_all2(i,j)=TF_ForContrast;
    T_all2(j,i)=TF_ForContrast;
      
    end
end
T_all1=zeros(416,416);
for i=1:416
    for j=i+1:416
     DependentVariable=new(Regressorstemp(:,6)<0.5,i,j);
    T = table(subnamenew, group, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group  + Age + Sex + Motion + Tiv +(1|Subject)');
T_all2(i,j)=mdl.Coefficients.tStat(2);
p_all2(i,j)=mdl.Coefficients.pValue(2);
T_all2(j,i)=mdl.Coefficients.tStat(2);
p_all2(j,i)=mdl.Coefficients.pValue(2);
      
    end
end



subname=subname(Regressors_allnewTemp(:,6)<0.5);
group=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,1);
age=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,3);
sex=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,2);
tiv=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,4);
fd=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,6);

T_all3=zeros(416,416);
for i=1:416
    for j=i+1:416
     DependentVariable=fc_all(Regressors_allnewTemp(:,6)<0.5,i,j);
%      data_harmonized = combat(DependentVariable', Regressorstemp(:,5)',[], 1);
    [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),[1,0,0,0,0,0],'T'); 
    T_all3(i,j)=TF_ForContrast;
    T_all3(j,i)=TF_ForContrast;
      
    end
end

T_all3=zeros(416,416);
for i=1:416
    for j=i+1:416
     DependentVariable=fc_all(Regressors_allnewTemp(:,6)<0.5,i,j);
    T = table(subnamenew, group, DependentVariable, age, sex, fd,tiv, ...
    'VariableNames', {'Subject', 'Group', 'Value', 'Age', 'Sex', 'Motion','Tiv'});
    mdl = fitlme(T, 'Value ~ Group  + Age + Sex + Motion + Tiv +(1|Subject)');
T_all2(i,j)=mdl.Coefficients.tStat(2);
p_all2(i,j)=mdl.Coefficients.pValue(2);
T_all2(j,i)=mdl.Coefficients.tStat(2);
p_all2(j,i)=mdl.Coefficients.pValue(2);
      Cohen_d(i)=mdl.Coefficients.Estimate(2)./sqrt(mdl.MSE);
    end
end

clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
Regressor=[Regressor_ELMCIsub2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=3) ,:);...
Regressor_NC(Regressor_NC(:,4)<0.5,:)];

Regressor(:,1)=[ones(length(find(Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=3)),1);...
-1*ones(length(find(Regressor_NC(:,4)<0.5)),1)];
new=cat(1,FC_elmci2((Regressor_ELMCIsub2(:,4)<0.5 & Regressor_ELMCIsub2(:,1)~=3),:,:),FC_NC2(Regressor_NC(:,4)<0.5,:,:));


new=new(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end


temp1=squeeze(mean(fc_all(Regressors_allnewTemp(:,1)==1 & Regressors_allnewTemp(:,6)<0.5,:,:),1));

temp1=squeeze(mean(FC_NC2,1));
temp2=squeeze(mean(fc_all(Regressors_allnewTemp(:,1)==-1 & Regressors_allnewTemp(:,6)<0.5,:,:),1));
temp1(isinf(temp1))=0;
temp2(isinf(temp2))=0;
corr2(temp1,temp2)

temp1=squeeze(mean(FC_elmci,1));
temp2=squeeze(mean(FC_elmci,1));
temp1(isinf(temp1))=0;
temp2(isinf(temp2))=0;
corr2(temp1,temp2)







%%
T=T_all(:,:,:,5);
for i=1:489
    maxV(i)=max(max(max(T(i,:,:))));
end

for sub=1:489
    sub
    temp1=squeeze(T(sub,:,:));
    temp2=squeeze(S_mci_nc(sub,:,:));
    maxT=max(temp1(:));
    if maxT==22
        num=0;
        nets=[];
        Unique=unique(temp1(:));
        for j=1:Unique(end-1)
            if j==1
                num=num+1;
                net=zeros(416,416);
                nets(:,:,num)=net;
                
            else
                num=num+1;
                net=zeros(416,416);
                net(temp1<=j)=temp2(temp1<=j);
                nets(:,:,num)=net;
            end
%             S{sub}=nets;
        end
            [dx,dy,dt]=gradient(nets,1,1,1);
            temp=mean(dt,3);
            gradient1(:,:,sub)=temp;
            
            [dx,dy,dt]=gradient(dt,1,1,1);
            temp2=mean(dt,3);
            gradient2(:,:,sub)=temp2;
            
        
    else
        num=0;
        nets=[];
        for j=1:maxT
            if j==1
                num=num+1;
                net=zeros(416,416);
                nets(:,:,num)=net;
            else
                num=num+1;
                net=zeros(416,416);
                net(temp1<=j)=temp2(temp1<=j);
                nets(:,:,num)=net;
            end

%             S{sub}=nets;
        end
%         gradient(nets,3);
            [dx,dy,dt]=gradient(nets,1,1,1);
            temp=mean(dt,3);
            gradient1(:,:,sub)=temp;
            
            [dx,dy,dt]=gradient(dt,1,1,1);
            temp2=mean(dt,3);
            gradient2(:,:,sub)=temp2;
    end
    
end
gradient1=permute(gradient1,[3,1,2]);
gradient2=permute(gradient2,[3,1,2]);

for i=1:size(gradient1,1)
    temp=squeeze(gradient1(i,:,:));
    gradient1_temp(i,:,:)=ZScoreMatrix(temp,1);
end
for i=1:size(gradient2,1)
    temp=squeeze(gradient2(i,:,:));
    gradient2_temp(i,:,:)=ZScoreMatrix(temp);
end

new=gradient2_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
T2map=[];
for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp(:,5)',Regressors_allnewTemp(:,[2,3,4,6]), 0);
    for i=1:416
        DependentVariable=data_harmonized(i,Regressors_allnewTemp(:,6)<0.5)';
        [b,r,SSE_1(i),SSR_1(i), T1, TF_ForContrast(i), Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
    end
        Ttemp=TF_ForContrast';
        T2map=[T2map;TF_ForContrast];
    df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),1)-...
        size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),2)
    PMap=2*(1-tcdf(abs(Ttemp), df));
    qThreshold=0.05;%设置p值
    SortP=sort(PMap);
    V=length(SortP);
    I=(1:V)';
    cVID = 1;
    cVN  = sum(1./(1:V));
    P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
    Thresholded=zeros(size(PMap));
    if ~isempty(P)
    Thresholded(find(PMap<=P))=1;
    end
    T_fdr(:,net)=Ttemp.*Thresholded;

    
end



[~,index]=sort(Regressors_allnewTemp(:,3));
for i=1:size(S_mci_nc,1)
    temp=squeeze(S_mci_nc(i,:,:));
    S_mci_nc_temp(i,:,:)=ZScoreMatrix(temp);
end
for i=1:size(S_elmci,1)
    temp=squeeze(S_elmci(i,:,:));
    S_elmci_temp(i,:,:)=ZScoreMatrix(temp);
end

path='I:\ADNI\NCsetdiff_Spreading_Smooth_new\*.mat';
file=dir(path);
for i=1:length(file)
   
    load([file(i).folder,'\',file(i).name]);
     S_ncsetdiff(i,:,:)=S(:,:,5);
end
for i=1:size(S_ncsetdiff,1)
    temp=squeeze(S_ncsetdiff(i,:,:));
    S_ncsetdiff_temp(i,:,:)=ZScoreMatrix(temp);
end






new=S_mci_nc(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new=S_mci_nc_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data', Regressors_allnewTemp(:,5)',[], 1);
data=reshape(data_harmonized',size(data_harmonized',1),18,18);



new=S_elmci_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
data_elmci=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data_elmci',Regressor_ELMCIsub1(:,5)',[], 1);
data_elmci=reshape(data_harmonized',size(data_harmonized',1),18,18);





new=S_mci_nc_temp;
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:8];
    mask=yeoindex;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
data=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data', Regressors_allnewTemp(:,5)',[], 1);
data=reshape(data_harmonized',size(data_harmonized',1),8,8);


[dt,dx,dy]=gradient(data,1,1,1);
[dt,dx,dy]=gradient(data(index,:,:),1,1,1);
[dt_mci,dx,dy]=gradient(data(Regressors_allnewTemp(:,1)==1,:,:),1,1,1);
[dt_nc,dx,dy]=gradient(data(Regressors_allnewTemp(:,1)==-1,:,:),1,1,1);

[dt,dx,dy]=gradient(data(index,:,:),1,1,1);
[dt2_mci,dx,dy]=gradient(gradient(data(Regressors_allnewTemp(:,1)==1,:,:),1,1,1),1,1,1);
[dt2_nc,dx,dy]=gradient(gradient(data(Regressors_allnewTemp(:,1)==-1,:,:),1,1,1),1,1,1);


gradient_magnitude=squeeze(sqrt(sum(dt.^2,1)));
figure,imagesc(gradient_magnitude);

Regressors_allnewTempSort=Regressors_allnewTemp(index,:);



dt=cat(1,dt_mci,dt_nc);
for i=1:8
for j=1:8
        DependentVariable=dt(Regressors_allnewTemp(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end


for i=1:8
for j=1:8
        DependentVariable=dt(Regressors_allnewTemp(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end

for i=1:8
for j=1:8
        DependentVariable=dt(Regressors_allnewTemp(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end


for i=1:18
for j=1:18
        DependentVariable=data(Regressors_allnewTemp(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end




[~,index_mci]=sort(Regressors_allnewTemp(Regressors_allnewTemp(:,1)==1,3));
[dt_mci,dx,dy]=gradient(data(index_mci,:,:),1,1,1);

[~,index_nc]=sort(Regressors_allnewTemp(Regressors_allnewTemp(:,1)==-1,3));
[dt_nc,dx,dy]=gradient(data(index_nc+248,:,:),1,1,1);

Regressors_allnewTempSort=Regressors_allnewTemp([index_mci;index_nc+248],:);

dt=cat(1,dt_mci,dt_nc);
for i=1:18
for j=1:18
        DependentVariable=dt(Regressors_allnewTempSort(:,6)<0.5,i,j);
        [b,r,SSE_1,SSR_1, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,...
   Regressors_allnewTempSort(Regressors_allnewTempSort(:,6)<0.5,[1,2,3,4,6,7]),Contrastnew(1:6),'T'); 
TF_ForContrast_brain_mean3(i,j)=TF_ForContrast;


end
end
boxplot(dt(Regressors_allnewTempSort(:,6)<0.5,18,17),Regressors_allnewTempSort(Regressors_allnewTempSort(:,6)<0.5,1));



Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),1)-...
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(8,8);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;

Ttemp=TF_ForContrast_brain_mean3(TF_ForContrast_brain_mean3~=0);
df=size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),1)-...
    size(Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[1,2,3,6,7]),2);
PMap=2*(1-tcdf(abs(Ttemp), df));
qThreshold=0.05;%设置p值
SortP=sort(PMap);
V=length(SortP);
I=(1:V)';
cVID = 1;
cVN  = sum(1./(1:V));
P   = SortP(find(SortP <= I/V*qThreshold/cVID, 1, 'last' ));
Thresholded=zeros(size(PMap));
if ~isempty(P)
Thresholded(find(PMap<=P))=1;
end
T_fdr=zeros(18,18);
T_fdr(TF_ForContrast_brain_mean3~=0)=Ttemp.*Thresholded;

G=digraph(T_fdr);
w=G.Edges.Weight; 
LWidths=w;
h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,18);
h.ArrowSize=15;
h.NodeColor=resultmap;
axis equal




G=digraph(T_fdr);
w=G.Edges.Weight; 
LWidths=w;
figure,h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,8);
h.ArrowSize=15;
h.NodeColor=map([1,2,3,4,5,6,7,9],:);
axis equal
colorbar;
caxis([-5,5])
ylim([-2,2])
xlim([-2,2])

h.NodeLabel='';
axis off
set(gcf,'color','w');
% 归一化到 [0.5, 0.9]
% abs_weights = abs(w);  % 取绝对值
% min_alpha = 0.2;
% max_alpha = 0.7;
% norm_alpha = (abs_weights - min(abs_weights)) / (max(abs_weights) - min(abs_weights));
% edge_alphas = norm_alpha * (max_alpha - min_alpha) + min_alpha;
% for i = 1:numedges(G)
%     h.EdgeAlpha(i) = edge_alphas(i);
% end

%%
data=data(Regressors_allnewTemp(:,6)<0.5,:,:);
labels=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,1);
labels(labels==-1)=0;

data=dt(Regressors_allnewTemp(:,6)<0.5,:,:);
labels=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,1);
labels(labels==-1)=0;


data=dt(Regressor(:,4)<0.5,:,:);
data=data(1:end,:,:);
labels= Regressor(Regressor(:,4)<0.5,1);
labels(labels==-1)=0;
labels=labels(1:end,:);
Regressor= Regressor(Regressor(:,4)<0.5,:);
Regressor=Regressor(1:end,:);


[M, N, ~] = size(data);
edge_mask = ~eye(N); % 所有非对角线的边

edge_mask = logical(logical(T_fdr)-diag(diag(logical(T_fdr))));

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);

X_all = X_all - [Regressor(:,[2,3,4,5]), ones(size(X_all,1),1)]...
    *([Regressor(:,[2,3,4,5]), ones(size(X_all,1),1)]\X_all);

%%

cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);        

C_list = 2.^(-20:1:20);
gamma_list = 2.^(-20:1:20);
best_acc = 0;
for c = C_list
    for g = gamma_list
        param = sprintf('-s 0 -t 2 -c %f -g %f -v 5 -q', c, g);  % 5折交叉验证
        acc = svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), param);
        if acc > best_acc
            best_acc = acc;
            best_c = c;
            best_g = g;
        end
    end
end

fprintf('最佳参数: C = %f, gamma = %f, 交叉验证准确率 = %.2f%%\n', best_c, best_g, best_acc);




param = sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
model=svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2),param);
[~, acc, prob_estimates] = svmpredict(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), model, '-q -b 1');  % 获取类别1的概率
[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTrain), prob_estimates(:, 1),1);
 figure,plot(X_roc, Y_roc);
 
 
 [~, acc, prob_estimates] = svmpredict(labels(idxTest), zscore(X_all(idxTest,:),[],2), model, '-q -b 1');  % 获取类别1的概率
[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob_estimates(:, 1),1);
 figure,plot(X_roc, Y_roc);
 
 %%

libsvm_options=sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
% === Step 1: 用完整数据训练一次模型 ===
X_all_z=zscore(X_all,[],2);
% model = svmtrain(labels, X_all_z, libsvm_options);
[y_pred_baseline,acc,~]=svmpredict(labels,X_all_z,model);

[~, acc, prob_estimates] = svmpredict(labels, zscore(X_all,[],2), model, '-q -b 1');  % 获取类别1的概率
[X_roc, Y_roc, ~, AUC] = perfcurve(labels, prob_estimates(:, 1),1);
figure,plot(X_roc, Y_roc);



% === Step 2: bootstrap扰动边，使用固定模型做预测 ===
n_runs=100;
n_bootstrap_samples=400;
% === 每条边扰动评估 ===
edge_importances = zeros(N, N);
for k = 1:n_edges

    for r = 1:n_runs
        idx = randsample(M, n_bootstrap_samples, true);
        X_sample = X_all_z(idx, :);
        y_sample = labels(idx);

        % 删除边
        X_sample(:, k) = 0;

        % 使用训练好的模型预测
        [y_pred_dropped, acc, ~] = svmpredict(y_sample, X_sample, model, '-q -b 1');



       y_pred_orig= y_pred_baseline(idx);
        diff=sum(y_pred_orig~=y_pred_dropped);
        diff_list(r)=diff;
    end

     [i, j] = ind2sub([N, N], edge_indices(k));
    edge_importances(i, j) = mean(diff_list);
end
edge_importances_z=zeros(8,8);
edge_importances_z=(edge_importances-(mean(edge_importances(edge_indices))))./(std(edge_importances(edge_indices)));
edge_importances_z(abs(edge_importances_z)<1.96)=0;

%%
[M, N, ~] = size(data);
edge_mask = ~eye(N); % 所有非对角线的边

edge_mask = logical(logical(T_fdr)-diag(diag(logical(T_fdr))));

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);

X_all = X_all - [Regressor(:,[2,3,4,5]), ones(size(X_all,1),1)]...
    *([Regressor(:,[2,3,4,5]), ones(size(X_all,1),1)]\X_all);

%%
cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);    
X_roc_all={};
Y_roc_all={};
AUC_all=[];

for net=1:8
    
[M, N, ~] = size(data);
edge_mask = zeros(N,N);
edge_mask(net,:)=1;
edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
    
    
C_list = 2.^(-20:1:20);
gamma_list = 2.^(-20:1:20);
best_acc = 0;
for c = C_list
    for g = gamma_list
        param = sprintf('-s 0 -t 2 -c %f -g %f -v 5 -q', c, g);  % 5折交叉验证
        acc = svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), param);
        if acc > best_acc
            best_acc = acc;
            best_c = c;
            best_g = g;
        end
    end
end

fprintf('最佳参数: C = %f, gamma = %f, 交叉验证准确率 = %.2f%%\n', best_c, best_g, best_acc);




param = sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
model=svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2),param);
% [~, acc, prob_estimates] = svmpredict(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTrain), prob_estimates(:, 1),1);
%  figure,plot(X_roc, Y_roc);
 
 
 [~, acc, prob_estimates] = svmpredict(labels(idxTest), zscore(X_all(idxTest,:),[],2), model, '-q -b 1');  % 获取类别1的概率
acc_final(net)=acc(1);
[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob_estimates(:, 1),1);
X_roc_all{net}=X_roc;
Y_roc_all{net}=Y_roc;
AUC_all(net)=AUC;

 figure,plot(X_roc, Y_roc);
end
 
cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);    
X_roc_all={};
Y_roc_all={};
AUC_all=[];

for net=1:8
    
[M, N, ~] = size(data);
edge_mask = zeros(N,N);
edge_mask(net,:)=1;
edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
  X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  
    
% C_list = 2.^(-20:1:20);
% gamma_list = 2.^(-20:1:20);
% best_acc = 0;
% for c = C_list
%     for g = gamma_list
%         param = sprintf('-s 0 -t 2 -c %f -g %f -v 5 -q', c, g);  % 5折交叉验证
%         acc = svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), param);
%         if acc > best_acc
%             best_acc = acc;
%             best_c = c;
%             best_g = g;
%         end
%     end
% end
% 
% fprintf('最佳参数: C = %f, gamma = %f, 交叉验证准确率 = %.2f%%\n', best_c, best_g, best_acc);

svm_model = fitcsvm(zscore(X_all(idxTrain,:),[],2), labels(idxTrain), ...
    'KernelFunction', 'rbf', ...
    'Standardize', false, ...
    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ShowPlots', false, ...
        'KFold', 5));


% param = sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
% model=svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2),param);
% [~, acc, prob_estimates] = svmpredict(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTrain), prob_estimates(:, 1),1);
%  figure,plot(X_roc, Y_roc);
 
%  
%  [~, acc, prob_estimates] = svmpredict(labels(idxTest), zscore(X_all(idxTest,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% acc_final(net)=acc(1);

% 对训练好的 SVM 模型做概率拟合
svm_model = fitPosterior(svm_model);

% 用测试集做预测，并输出概率
[label_pred, prob_estimates] = predict(svm_model, zscore(X_all(idxTest,:), [], 2));
true_labels = labels(idxTest);  % 真实标签
acc(net) = mean(label_pred == true_labels) * 100;  % 百分比形式（可去掉 *100）

[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob_estimates(:, 2),1);
X_roc_all{net}=X_roc;
Y_roc_all{net}=Y_roc;
AUC_all(net)=AUC;

 figure,plot(X_roc, Y_roc);
end


cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);    
X_roc_all={};
Y_roc_all={};
AUC_all=[];

for net=1:18
    
[M, N, ~] = size(data);
edge_mask = zeros(N,N);
edge_mask(:,net)=1;
edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
%   X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
%     *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  
    

svm_model = fitcsvm(X_all(idxTrain,:), labels(idxTrain), ...
    'KernelFunction', 'rbf', ...
    'Standardize', false, ...
    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ShowPlots', false, ...
        'KFold', 5));



% 对训练好的 SVM 模型做概率拟合
svm_model = fitPosterior(svm_model);

% 用测试集做预测，并输出概率
[label_pred, prob_estimates] = predict(svm_model, X_all(idxTest,:));
true_labels = labels(idxTest);  % 真实标签
acc(net) = mean(label_pred == true_labels) * 100;  % 百分比形式（可去掉 *100）

[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob_estimates(:, 2),1);
X_roc_all{net}=X_roc;
Y_roc_all{net}=Y_roc;
AUC_all(net)=AUC;

%  figure,plot(X_roc, Y_roc);

% nb_model = fitcnb(X_all(idxTrain,:), labels(idxTrain));
% cv_nb_model = crossval(nb_model, 'KFold', 5);
% cvLoss = kfoldLoss(cv_nb_model);
% fprintf('5-fold Cross-Validated Classification Error (Naive Bayes): %.2f%%\n', cvLoss*100);
% 
% 
% [label_pred, prob_estimates] = predict(nb_model, X_all(idxTest,:));
% true_labels = labels(idxTest);
% acc(net) = mean(label_pred == true_labels) * 100; 
% [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob_estimates(:,2), 1);
% X_roc_all{net}=X_roc;
% Y_roc_all{net}=Y_roc;
% AUC_all(net)=AUC;




% rf_model = fitcensemble(X_all(idxTrain,:), labels(idxTrain), ...
%     'Method', 'Bag', ...
%     'NumLearningCycles', 100, ...
%     'Learners', templateTree('MaxNumSplits', 10));
% 
% % K折交叉验证
% cv_rf_model = crossval(rf_model, 'KFold', 5);
% cvLoss = kfoldLoss(cv_rf_model);
% fprintf('5-fold Cross-Validated Classification Error (Random Forest): %.2f%%\n', cvLoss * 100);
% 
% % 测试集预测
% [label_pred, prob_estimates] = predict(rf_model, X_all(idxTest,:));
% 
% % 准确率
% true_labels = labels(idxTest);
% acc(net) = mean(label_pred == true_labels) * 100;
% 
% % AUC 和 ROC 曲线
% [X_roc, Y_roc, ~, AUC] = perfcurve(true_labels, prob_estimates(:,2), 1);
% X_roc_all{net} = X_roc;
% Y_roc_all{net} = Y_roc;
% AUC_all(net) = AUC;

end

% cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
% idxTrain = training(cv);   
% idxTest = test(cv);   

new=cat(1,S_elmci_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17),S_ncsetdiff_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17));
new(isnan(new))=0;
clear data_mat_mean_all TF_ForContrast_brain_mean3 T_all data_mat_mean
Regressors=cat(1,Regressor_ELMCIsub1(:,1:6),NCnewCov);
Regressors(:,1)=[ones(size(Regressor_ELMCIsub1(:,1:6),1),1);zeros(size(NCnewCov,1),1)];
for sub=1:size(new,1)
    data_mat=squeeze(new(sub,:,:));
    order=[1:18];
    mask=yeoindex17;
    t=[];
    start = 1;
    lines = 1;
    for i = 1:length(order)
        add = find(mask==order(i));
        t = [t;add];
        start = start + length(add);
        lines(i+1) = start;
    end
    data_reorder = data_mat(t,t);
    idx_begin = lines(1:end-1);
    idx_end = lines(2:end)-1;
    for i = 1:length(lines)-1
        for j = 1:length(lines)-1
        data_temp = data_reorder(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j));
%         if i == j
%         data_temp = convet_matrix_to_vector(data_temp);
%         end
        data_temp = data_temp(:);
        data_mat_mean(i,j) = mean(data_temp);
    
        end
    end
    data_mat_mean_all(sub,:,:)=data_mat_mean;
end
data_Test=reshape(data_mat_mean_all,size(data_mat_mean_all,1),[]);
data_harmonized = combat(data_Test',Regressors(:,5)',[], 1);
data_Test=reshape(data_harmonized',size(data_harmonized',1),18,18);
data_Test=data_mat_mean_all;

indexNC=randperm(374 - 136 + 1, 135) + 135;
index=[1:135];
ttest(Regressors(index,3),Regressors(indexNC,3))
for net=1:18
    
[M, N, ~] = size(data);
edge_mask = zeros(N,N);
edge_mask(net,:)=1;
edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end

[M1, N1, ~] = size(data_Test);
X_test=zeros(M1, n_edges);
for i = 1:M1
    mat = squeeze(data_Test(i, :, :));
    X_test(i, :) = mat(edge_mask);
end

%   X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
%     *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  
    
% C_list = 2.^(-20:1:20);
% gamma_list = 2.^(-20:1:20);
% best_acc = 0;
% for c = C_list
%     for g = gamma_list
%         param = sprintf('-s 0 -t 2 -c %f -g %f -v 5 -q', c, g);  % 5折交叉验证
%         acc = svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), param);
%         if acc > best_acc
%             best_acc = acc;
%             best_c = c;
%             best_g = g;
%         end
%     end
% end
% 
% fprintf('最佳参数: C = %f, gamma = %f, 交叉验证准确率 = %.2f%%\n', best_c, best_g, best_acc);

svm_model = fitcsvm(X_all(:,:), labels(:), ...
    'KernelFunction', 'rbf', ...
    'Standardize', true, ...
    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ShowPlots', false, ...
        'KFold', 5));


% param = sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
% model=svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2),param);
% [~, acc, prob_estimates] = svmpredict(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTrain), prob_estimates(:, 1),1);
%  figure,plot(X_roc, Y_roc);
 
%  
%  [~, acc, prob_estimates] = svmpredict(labels(idxTest), zscore(X_all(idxTest,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% acc_final(net)=acc(1);

% 对训练好的 SVM 模型做概率拟合
svm_model = fitPosterior(svm_model);

% 用测试集做预测，并输出概率
[label_pred, prob_estimates] = predict(svm_model, X_test([index;indexNC],:));
true_labels = Regressors([index;indexNC],1);  % 真实标签
acc(net) = mean(label_pred == true_labels) * 100;  % 百分比形式（可去掉 *100）

[X_roc, Y_roc, ~, AUC] = perfcurve(Regressors([index;indexNC],1), prob_estimates(:, 2),1);
X_roc_all{net}=X_roc;
Y_roc_all{net}=Y_roc;
AUC_all(net)=AUC;

%  figure,plot(X_roc, Y_roc);
end
 figure,plot(X_roc_all{14}, Y_roc_all{14});

indexNC=randperm(374 - 136 + 1, 135) + 135;
index=[1:135];
ttest(Regressors(index,3),Regressors(indexNC,3))
[M, N, ~] = size(data);
edge_mask = ~eye(N);

edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end

[M1, N1, ~] = size(data_Test);
X_test=zeros(M1, n_edges);
for i = 1:M1
    mat = squeeze(data_Test(i, :, :));
    X_test(i, :) = mat(edge_mask);
end

svm_model = fitcsvm(X_all(:,:), labels(:), ...
    'KernelFunction', 'rbf', ...
    'Standardize', true, ...
    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ShowPlots', false, ...
        'KFold', 5));



% 对训练好的 SVM 模型做概率拟合
svm_model = fitPosterior(svm_model);

% 用测试集做预测，并输出概率
[label_pred, prob_estimates] = predict(svm_model, X_test([index;indexNC],:));
true_labels = Regressors([index;indexNC],1);  % 真实标签
acc = mean(label_pred == true_labels) * 100;  % 百分比形式（可去掉 *100）

[X_roc, Y_roc, ~, AUC] = perfcurve(Regressors([index;indexNC],1), prob_estimates(:, 2),1);
figure,plot(X_roc, Y_roc);




X_roc_all={};
Y_roc_all={};
AUC_all=[];
acc=[];
temp=Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,:);
index=randperm(461 - 232 + 1, 135) + 231;
dataTest=cat(1,data_elmci,[data(index,:,:)]);
Regressors=cat(1,Regressor_ELMCIsub1(:,[1,2,3,6,4]),temp(index,[1,2,3,4,6]));
Regressors(:,1)=[ones(size(Regressor_ELMCIsub1(:,[1,2,3,6,4]),1),1);zeros(size(temp(index,[1,2,3,4,6]),1),1)];
for net=1:18
    
[M, N, ~] = size(data);
edge_mask = zeros(N,N);
edge_mask(net,:)=1;
edge_mask=logical(edge_mask);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end

[M1, N1, ~] = size(dataTest);
X_test=zeros(M1, n_edges);
for i = 1:M1
    mat = squeeze(dataTest(i, :, :));
    X_test(i, :) = mat(edge_mask);
end

%   X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]...
%     *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  
    
% C_list = 2.^(-20:1:20);
% gamma_list = 2.^(-20:1:20);
% best_acc = 0;
% for c = C_list
%     for g = gamma_list
%         param = sprintf('-s 0 -t 2 -c %f -g %f -v 5 -q', c, g);  % 5折交叉验证
%         acc = svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), param);
%         if acc > best_acc
%             best_acc = acc;
%             best_c = c;
%             best_g = g;
%         end
%     end
% end
% 
% fprintf('最佳参数: C = %f, gamma = %f, 交叉验证准确率 = %.2f%%\n', best_c, best_g, best_acc);

svm_model = fitcsvm(X_all(:,:), labels(:), ...
    'KernelFunction', 'rbf', ...
    'Standardize', false, ...
    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ShowPlots', false, ...
        'KFold', 5));


% param = sprintf('-s 0 -t 2 -c %f -g %f -b 1 -q', best_c, best_g);
% model=svmtrain(labels(idxTrain), zscore(X_all(idxTrain,:),[],2),param);
% [~, acc, prob_estimates] = svmpredict(labels(idxTrain), zscore(X_all(idxTrain,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTrain), prob_estimates(:, 1),1);
%  figure,plot(X_roc, Y_roc);
 
%  
%  [~, acc, prob_estimates] = svmpredict(labels(idxTest), zscore(X_all(idxTest,:),[],2), model, '-q -b 1');  % 获取类别1的概率
% acc_final(net)=acc(1);

% 对训练好的 SVM 模型做概率拟合
svm_model = fitPosterior(svm_model);

% 用测试集做预测，并输出概率
[label_pred, prob_estimates] = predict(svm_model, X_test);
true_labels = Regressors(:,1);  % 真实标签
acc(net) = mean(label_pred == true_labels) * 100;  % 百分比形式（可去掉 *100）

[X_roc, Y_roc, ~, AUC] = perfcurve(Regressors(:,1), prob_estimates(:, 1),1);
X_roc_all{net}=X_roc;
Y_roc_all{net}=Y_roc;
AUC_all(net)=AUC;

%  figure,plot(X_roc, Y_roc);
end
 figure,plot(X_roc_all{14}, Y_roc_all{14});

 
 
%%
cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);    
    
[M, N, ~] = size(data);
edge_mask = ~eye(N)

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
  X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5 ,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  


[B, FitInfo] = lassoglm(X_all(idxTrain,:), labels(idxTrain), 'binomial', 'CV', 5);
bestB = B(:, FitInfo.IndexMinDeviance);
idxImportant = find(bestB ~= 0);  % 找到非零的连接（重要特征）

X_train_sel = X_all(idxTrain, idxImportant);
X_test_sel = X_all(idxTest, idxImportant);

model = fitglm(X_train_sel, labels(idxTrain), 'Distribution', 'binomial');

prob = predict(model, X_test_sel);
y_pred = double(prob >= 0.5);

accuracy = mean(y_pred == labels(idxTest));
fprintf('测试集准确率: %.2f%%\n', accuracy * 100);

[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob, 1);
fprintf('AUC值: %.3f\n', AUC);



idxImportantnew=edge_indices(idxImportant);

[x,y]=ind2sub([18,18],idxImportantnew);
Weights = bestB(idxImportant);
nets=zeros(18,18);
for i=1:length(x)
   
        
        nets(x(i),y(i))=Weights(i);

end
G=digraph(nets);
w=G.Edges.Weight;
LWidths=5*w./max(w);
figure,h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,18);
h.ArrowSize=15;
axis equal
 






cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv);    
    
[M, N, ~] = size(data);
edge_mask =logical(ones(N,N));

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
  X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5 ,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  


[B, FitInfo] = lassoglm(X_all(idxTrain,:), labels(idxTrain), 'binomial', 'CV', 5);
bestB = B(:, FitInfo.IndexMinDeviance);
idxImportant = find(bestB ~= 0);  % 找到非零的连接（重要特征）

X_train_sel = X_all(idxTrain, idxImportant);
X_test_sel = X_all(idxTest, idxImportant);

model = fitglm(X_train_sel, labels(idxTrain), 'Distribution', 'binomial');

prob = predict(model, X_test_sel);
y_pred = double(prob >= 0.5);

accuracy = mean(y_pred == labels(idxTest));
fprintf('测试集准确率: %.2f%%\n', accuracy * 100);

[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob, 1);
fprintf('AUC值: %.3f\n', AUC);

[x,y]=ind2sub([8,8],idxImportant);
Weights = bestB(idxImportant);
nets=zeros(8,8);
for i=1:length(x)
   
        
        nets(x(i),y(i))=Weights(i);

end
G=digraph(nets);
w=G.Edges.Weight;
LWidths=5*w./max(w);
figure,h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,8);
h.ArrowSize=15;
axis equal





%%
n_bootstrap = 100;        % 重复次数
sample_fraction = 0.8;    % 每次抽取70%样本
  
    
[M, N, ~] = size(data);
edge_mask = ~eye(N)

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
  X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5 ,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  

selection_count = zeros(n_edges, 1);
% Stability Selection 主循环
for k = 1:n_bootstrap
    k
    % Bootstrap 抽样
    idx_boot = randsample(M, round(sample_fraction * M), true);
    X_boot = X_all(idx_boot, :);
    y_boot = labels(idx_boot);

    % Lasso
    [B, FitInfo] = lassoglm(X_boot, y_boot, 'binomial', 'CV', 5);

    % 最优lambda下的系数
    bestB = B(:, FitInfo.IndexMinDeviance);
    selection_count = selection_count + (bestB ~= 0);
end
% 计算每个特征被选中的频率
selection_freq = selection_count / n_bootstrap;

% 阈值选择（例如频率 > 0.6）
threshold = 0.6;
idxImportant = find(selection_freq >= threshold);


cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv); 

model = fitclinear(X_all(idxTrain, idxImportant),labels(idxTrain), 'Learner', 'logistic');
[y_pred,prob] = predict(model, X_all(idxTest, idxImportant));
acc = mean(y_pred == labels(idxTest));
fprintf('用稳定连接训练模型的测试准确率: %.2f%%\n', acc * 100);


[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob(:,2), 1);
fprintf('AUC值: %.3f\n', AUC);


idxImportantnew=edge_indices(idxImportant);

[x,y]=ind2sub([18,18],idxImportantnew);
Weights = selection_freq(idxImportant);
nets=zeros(18,18);
for i=1:length(x)
   
        
        nets(x(i),y(i))=Weights(i);

end
G=digraph(nets);
w=G.Edges.Weight;
LWidths=5*w./max(w);
figure,h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,18);
h.ArrowSize=15;
axis equal




G=digraph(T_fdr);
w=G.Edges.Weight;
LWidths=5*w./max(w);
figure,h=plot(G,'LineWidth',abs(LWidths),'Layout','circle','EdgeCData', w);
load('cmap.mat','cmap')
colormap(cmap)
h.MarkerSize=10*ones(1,18);
h.ArrowSize=15;
axis equal





data=S_mci_nc_temp(Regressors_allnewTemp(:,6)<0.5,:,:);
n_bootstrap = 100;        % 重复次数
sample_fraction = 0.8;    % 每次抽取70%样本
  
    
[M, N, ~] = size(data);
edge_mask = ~eye(N);

edge_indices = find(edge_mask);
n_edges = numel(edge_indices);

X_all = zeros(M, n_edges);
for i = 1:M
    mat = squeeze(data(i, :, :));
    X_all(i, :) = mat(edge_mask);
end
  X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5 ,[2,3,4,6]), ones(size(X_all,1),1)]...
    *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  

selection_count = zeros(n_edges, 1);
% Stability Selection 主循环
for k = 1:n_bootstrap
    k
    % Bootstrap 抽样
    idx_boot = randsample(M, round(sample_fraction * M), true);
    X_boot = X_all(idx_boot, :);
    y_boot = labels(idx_boot);

    % Lasso
    [B, FitInfo] = lassoglm(X_boot, y_boot, 'binomial', 'CV', 5);

    % 最优lambda下的系数
    bestB = B(:, FitInfo.IndexMinDeviance);
    selection_count = selection_count + (bestB ~= 0);
end
% 计算每个特征被选中的频率
selection_freq = selection_count / n_bootstrap;

% 阈值选择（例如频率 > 0.6）
threshold = 0.6;
idxImportant = find(selection_freq >= threshold);


cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
idxTrain = training(cv);   
idxTest = test(cv); 

model = fitclinear(X_all(idxTrain, idxImportant),labels(idxTrain), 'Learner', 'logistic');
[y_pred,prob] = predict(model, X_all(idxTest, idxImportant));
acc = mean(y_pred == labels(idxTest));
fprintf('用稳定连接训练模型的测试准确率: %.2f%%\n', acc * 100);


[X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob(:,2), 1);
fprintf('AUC值: %.3f\n', AUC);


idxImportantnew=edge_indices(idxImportant);

[x,y]=ind2sub([18,18],idxImportantnew);
Weights = selection_freq(idxImportant);
nets=zeros(18,18);
for i=1:length(x)
   
        
        nets(x(i),y(i))=Weights(i);

end

%%
new=S_mci_nc_temp(:,yeo_17labels_7_to_17,yeo_17labels_7_to_17);
new(isnan(new))=0;

for net=1:18
    network=new(:,yeoindex17==net,:);

    network=mean(network,2);
    network=squeeze(network);
    data_harmonized = combat(network', Regressors_allnewTemp(:,5)',Regressors_allnewTemp(:,[2,3,4,6]), 0);
       data=data_harmonized';
       data=data(Regressors_allnewTemp(:,6)<0.5,:,:);
     
        n_bootstrap = 100;        % 重复次数
        sample_fraction = 0.8;    % 每次抽取70%样本


        [M, N, ~] = size(data);
          X_all=data;
          X_all = X_all - [Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5 ,[2,3,4,6]), ones(size(X_all,1),1)]...
            *([Regressors_allnewTemp(Regressors_allnewTemp(:,6)<0.5,[2,3,4,6]), ones(size(X_all,1),1)]\X_all);  
           % 阈值选择（例如频率 > 0.6）
        idxImportant = find(T2map(:,net));  
    cv = cvpartition(size(data,1), 'HoldOut', 0.3); 
    idxTrain = training(cv);   
    idxTest = test(cv); 

    model = fitclinear(X_all(idxTrain, idxImportant),labels(idxTrain), 'Learner', 'logistic');
    [y_pred,prob] = predict(model, X_all(idxTest, idxImportant));
    acc = mean(y_pred == labels(idxTest));
    fprintf('用稳定连接训练模型的测试准确率: %.2f%%\n', acc * 100);


    [X_roc, Y_roc, ~, AUC] = perfcurve(labels(idxTest), prob(:,2), 1);
    fprintf('AUC值: %.3f\n', AUC);     


    
end

colors=tab20c(20);
colors=colors([1:4:16],:);
h.figure = figure('Color','white','Units','normalized','Position',[0 0 1,1]);
nRows = 2;
nCols = 9;
total_axes = nRows * nCols;
% 间距参数（相对单位）
hSpacing = 0.025;  % 水平间距
vSpacing = 0.04;  % 垂直间距
marginL = 0.02;
marginB = 0.05;
% 每个 axes 的宽高
axWidth = (1 - 2*marginL - (nCols - 1) * hSpacing) / nCols;
axHeight = (1 - 2*marginB - (nRows - 1) * vSpacing) / nRows;

h.axes = gobjects(total_axes,1);  % 初始化 axes 句柄数组
for k = 1:total_axes
    row = floor((k - 1) / nCols);        
    col = mod((k - 1), nCols);           


    posX = marginL + col * (axWidth + hSpacing);
    posY = 1 - marginB - (row + 1) * axHeight - row * vSpacing;

    h.axes(k) = axes('Parent', h.figure, ...
                     'Units', 'normalized', ...
                     'Position', [posX, posY, axWidth, axHeight]);
    axes(h.axes(k)); hold on;
    plot(X_roc_o_y{k}, Y_roc_o_y{k}, 'Color', colors(1,:), 'LineWidth', 2);
    plot(X_roc_o_m{k}, Y_roc_o_m{k}, 'Color', colors(2,:), 'LineWidth', 2);
    plot(X_roc_m_y{k}, Y_roc_m_y{k}, 'Color', colors(3,:), 'LineWidth', 2);
    plot(X_roc_all{k}, Y_roc_all{k}, 'Color', colors(4,:), 'LineWidth', 2);
    axis([0 1 0 1]);
%     grid on;

    set(h.axes(k), ...
        'FontSize', 14, ...
        'LineWidth',2, ...
        'XTick', 0:0.2:1, ...
        'YTick', 0:0.2:1,...
        'FontName','Times New Roman',...
        'FontWeight','bold');

end

h.figure = figure('Color','white', 'Units','normalized', 'Position',[0.1 0.1 0.8 0.7]);
nRows = 2;
nCols = 4;
total_axes = nRows * nCols;
% 间距和边距设置
hSpacing = 0.03;
vSpacing = 0.05;
marginL = 0.05;
marginR = 0.05;
marginT = 0.05;
marginB = 0.07;
axWidth = (1 - marginL - marginR - (nCols - 1) * hSpacing) / nCols;
axHeight = (1 - marginT - marginB - (nRows - 1) * vSpacing) / nRows;


h.axes = gobjects(total_axes,1);  % 初始化 axes 句柄数组
for k = 1:total_axes
    row = floor((k - 1) / nCols);        
    col = mod((k - 1), nCols);           


    posX = marginL + col * (axWidth + hSpacing);
    posY = 1 - marginB - (row + 1) * axHeight - row * vSpacing;

    h.axes(k) = axes('Parent', h.figure, ...
                     'Units', 'normalized', ...
                     'Position', [posX, posY, axWidth, axHeight]);
    axes(h.axes(k)); hold on;
    plot(X_roc_o_y{k}, Y_roc_o_y{k}, 'Color', colors(1,:), 'LineWidth', 2);
    plot(X_roc_o_m{k}, Y_roc_o_m{k}, 'Color', colors(2,:), 'LineWidth', 2);
    plot(X_roc_m_y{k}, Y_roc_m_y{k}, 'Color', colors(3,:), 'LineWidth', 2);
    plot(X_roc_all{k}, Y_roc_all{k}, 'Color', colors(4,:), 'LineWidth', 2);
    axis([0 1 0 1]);
%     grid on;

    set(h.axes(k), ...
        'FontSize', 14, ...
        'LineWidth',2, ...
        'XTick', 0:0.2:1, ...
        'YTick', 0:0.2:1,...
        'FontName','Times New Roman',...
        'FontWeight','bold');

end

