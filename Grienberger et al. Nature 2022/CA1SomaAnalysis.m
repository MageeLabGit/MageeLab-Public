%% run section by section

%input variables
S.parentfolder=''; 
 
if ~isfolder(S.parentfolder)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', S.parentfolder);
        uiwait(warndlg(errorMessage));
        return;  
end

cd (S.parentfolder);
d=dir('*.h5'); S.fileNames={d(:).name}; 
wsCell=cell(1,length(d)); 

for i=1:length(d)
    data=ws.loadDataFile(S.fileNames{i});
    S.acq=data.header.AcquisitionSampleRate;
    names=fieldnames(data);
    ws1=[]; ws1=data.(names{2}).analogScans;
   
    idxl = ws1(:,1)>=2;
    idxl(1) = 0;
    idx = find(idxl);
    yest = ws1(idx-1,1)<2; 
    startlevels=idx(yest)-1; 
       
    if i>1
        wsCell{1,i}=ws1(:,1:7);
    else        
        wsCell{1,i}=ws1(startlevels(1):length(ws1),1:7);
    end
    
    S.wsALL=vertcat(S.wsALL, (wsCell{1,i})); 
end


idxl = S.wsALL(:,5)>=2;
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,5)<2; 
S.frametiming=idx(yest);
length(S.frametiming)

idxl = S.wsALL(:,1)>=2;
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,1)<2; 
S.startlevels=idx(yest)-1;
S.startlevels(length(S.startlevels)+1,1)=length(S.wsALL(:,1))+1; 

% divide up into individual laps and calculate position signal from analog running signal 

S.wsL=cell(1,size(S.wsALL,2));

test=cell(1,length(S.startlevels)-1); 

for k=1:size(S.wsALL,2)
    for i=1:(length(S.startlevels)-1)
     test{:,i}=S.wsALL(S.startlevels(i):S.startlevels(i+1)-1,k); 
    end
    S.wsL{1,k}=test; 
end

S.maxpos=zeros(length(S.startlevels),1); 
position=S.wsL{1,2};
for i=1:(length(S.startlevels)-1)
    position{1,i}=cumtrapz(position{1,i}); 
    S.maxpos(i,1)=max(position{1,i}); 
end

S.wsL{1,size(S.wsALL,2)+1}=position; 

interim=[]; 
for k=1:length(S.startlevels)-1
    interim=vertcat(interim, (position{1,k})); 
end

S.wsALL=horzcat(S.wsALL, interim); 

cd (S.parentfolder);
d=dir('*.h5'); S.fileNames={d(:).name}; 
wsCell=cell(1,length(d)); 
interim=[];
for i=1:length(d)
    data=ws.loadDataFile(S.fileNames{i});
    names=fieldnames(data);
    ws1=data.(names{2}).analogScans;

    idxl = ws1(:,1)>=2;
    idxl(1) = 0;
    idx = find(idxl);
    yest = ws1(idx-1,1)<2; 
    startlevels=idx(yest)-1; 
    
    
    if i>1
        if size(ws1,2)>7
            wsCell{1,i}=ws1(:,8);
        else
            wsCell{1,i}=ws1(:,2);
            wsCell{1,i}(:)=0;
        end    
    else
        if size(ws1,2)>7
            wsCell{1,i}=ws1(startlevels(1):length(ws1),8);
        else
            wsCell{1,i}=ws1(startlevels(1):length(ws1),2);
            wsCell{1,i}(:)=0;
        end
    end

    interim=vertcat(interim, (wsCell{1,i}));
end

S.wsALL(:,9)=interim;
test=cell(1,length(S.startlevels)-1); 

for k=9
    for i=1:(length(S.startlevels)-1)
        test{:,i}=S.wsALL(S.startlevels(i):S.startlevels(i+1)-1,k); 
    end
    S.wsL{1,k}=test; 
end

S.lightposition=zeros((length(S.startlevels)-2),2); S.lightposition(:,:)=nan; 
for k=1:(length(S.startlevels)-2)
    lap=S.wsL{1,8}{1,k};
    a=find(S.wsL{1,9}{1,k}>2,1); 
    if (a>0)
        S.lightposition(k,1)=floor(lap(a)/((max(lap)/50))); 
        if S.lightposition(k,1)<1
            S.lightposition(k,1)=1;
        end
    end
    a=find(S.wsL{1,9}{1,k}>2,1,'last'); 
    if (a>0)
        S.lightposition(k,2)=ceil(lap(a)/((max(lap)/50)));
    end
end

control=[S.maxpos S.startlevels]; 

figure; plot(control(:,1)); 

S.frametimingOriginal=S.frametiming;
for k=1:length(S.frametiming)
    if S.wsALL(S.frametiming(k,1),2)<=0.05
        c=find(S.wsALL(S.frametiming(k,1):S.framelength(k,2),2)>0.05, 1,'first');
            if isempty(c)==0 && c<(0.032*S.acq) %length of 1 frame
                S.frametiming(k,1)=S.frametiming(k,1)+c-1;
            end
    end
end

clearvars -except S control

%check variable 'control' for problems with startmarkers


%% ROIs extracted via SIMA package; remove shutter artifact
% load dataset.txt

S.lapsused=1:(length(S.startlevels)-2);
S.medfiltC=500; 
S.polynom=3; 
S.medfiltSm=5;

S.datasetorig=dataset;

startframe=1:400:size(dataset,1);
datasetSm=dataset; 
for i=1:size(dataset,2)
    datasetSm(:,i)=medfilt1(dataset(:,i),S.medfiltSm); 

    for k=1:length(startframe)
        dataset(startframe(k):startframe(k)+1,i)=NaN; 
        datasetSm(startframe(k):startframe(k)+1,i)=NaN; 
    end

end

S.dataset=dataset;
S.datasetSm=datasetSm; 

S.frametimesperlap=cell(1,(length(S.startlevels)-2));
S.dataperlap=cell(1,size(dataset,2));

for i=2:size(dataset,2)
    test=dataset(:,i); 
    for j=S.lapsused
        C=S.frametiming>S.startlevels(j)&S.frametiming<S.startlevels(j+1)-1; 
        S.frametimesperlap{1,j}=S.frametiming(C); 
        S.frametimesperlap{1,j}=S.frametimesperlap{1,j}-S.startlevels(j); 
        D=find(S.frametiming>y1&S.frametiming<y2); 
        S.dataperlap{1,i}{1,j}=test(D); 
     end
end


%%

close all
C=S.wsALL(:,2)>0.05; 
allrunning=S.wsALL(:,2); S.running.allrunningR=allrunning(C); 
allstart=S.wsALL(:,1); S.running.allstartR=allstart(C); 
S.running.allstartR(1,1)=0;S.running.allstartR(2:10,1)=5;
idxl = S.running.allstartR>=2;
idxl(1) = 0;
idx = find(idxl);
yest = S.running.allstartR(idx-1,1)<2; 
S.startlevelsR=idx(yest)-1;
S.startlevelsR(length(S.startlevelsR)+1,1)=length(S.running.allstartR)+1;

allposition=S.wsALL(:,8); S.running.allpositionR=allposition(C);
allvalve=S.wsALL(:,3); S.running.allvalveR=allvalve(C); 
alllicking=S.wsALL(:,4); S.running.alllickingR=alllicking(C);
alllight=S.wsALL(:,9); S.running.allightR=alllight(C); 

S.wsR=cell(1,6);

test=cell(1,length(S.startlevelsR)-1); 
for i=1:(length(S.startlevelsR)-1)
    test{:,i}=S.running.allstartR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end

S.wsR{1,1}=test; 

test2=cell(1,length(S.startlevelsR)-1); 
for i=1:(length(S.startlevelsR)-1)
    test2{:,i}=S.running.alllickingR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end

S.wsR{1,3}=test2; 

for i=1:(length(S.startlevelsR)-1)
    test{:,i}=S.running.allpositionR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end
S.wsR{1,2}=test; 

test3=cell(1,length(S.startlevelsR)-1); 
for i=1:(length(S.startlevelsR)-1)
    test3{:,i}=S.running.allvalveR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end

S.wsR{1,4}=test3;

test4=cell(1,length(S.startlevelsR)-1); 
for i=1:(length(S.startlevelsR)-1)
    test4{:,i}=S.running.allrunningR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end

S.wsR{1,5}=test4;

test5=cell(1,length(S.startlevelsR)-1); 
for i=1:(length(S.startlevelsR)-1)
    test5{:,i}=S.running.allightR(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
end

S.wsR{1,6}=test5;


S.maxposR=zeros(length(S.startlevelsR),1); 
for i=1:(length(S.startlevelsR)-1)
    cg=test{1,i};
    cg=cg(2:length(test{1,i}),1);
    S.maxposR(i,1)=max(cg); 
end

control=[control S.maxposR S.startlevelsR]; 

%% is number of startlevelsR and startlevels the same? 

S.threshold=zeros(50,length(S.startlevelsR)-2); % must have 51 rows at the end
 
for i=1:(length(S.startlevelsR)-2)
    binsize=(S.maxposR(i,1)/50); 
    cg=test{1,i};
    for binnumber=1:50
        thres=binnumber * binsize; 
        idxl = cg>=thres;
        idxl(1) = 0;
        idx = find(idxl, 1);
        yest = cg(idx-1,1)<thres; 
        if idx(yest)>0
            S.threshold(binnumber,i)=idx(yest);
        else
            S.threshold(binnumber,i)=numel(cg);
        end
    end
end

interim=zeros(1,(length(S.startlevelsR))-2); interim(:,:)=1; 
S.threshold=vertcat(interim, S.threshold);

S.Lrate=zeros(50,(length(S.startlevelsR)-2)); S.Lrate(:,:)=NaN; 
S.Ltick=cell(1,length(S.startlevelsR)-2); 
for j=S.lapsused
    cg=test2{:,j}; 
    for binnumber=1:50
        cgg=cg(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1);
        idxl = cgg>=2; 
        idxl(1) = 0;
        idx = find(idxl);
        yest = cgg(idx-1,1)<2; 
        if numel(idx(yest))>0
            S.Lrate(binnumber,j)=numel(idx(yest))/(length(cgg)/S.acq); 
        else
            S.Lrate(binnumber,j)=0;  
        end
    end
    %licks as a function of position; for each lick determine position; 
    idxl = cg>=2; 
    idxl(1) = 0;
    idx = find(idxl);
    yest = cg(idx-1,1)<2; 
    licktiming=idx(yest);
    cgg=test{:,j};
    S.Ltick{1,j}=cgg(licktiming); 
end

S.MeanLrate=mean(S.Lrate, 2, 'omitnan');


S.V=zeros(50,(length(S.startlevelsR)-2)); S.V(:,:)=nan; 
S.Vtick=cell(1,length(S.startlevelsR)-2); 
for j=S.lapsused
    cg=test3{:,j}; 
    for binnumber=1:50
        cgg=cg(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1);
        idxl = cgg>=2; 
        idxl(1) = 0;
        idx = find(idxl);
        yest = cgg(idx-1,1)<2; 
        if numel(idx(yest))>0
            S.V(binnumber,j)=numel(idx(yest)); 
        else
            S.V(binnumber,j)=0;    
        end
    end
    %reward location as a function of position; for each lick determine position; 
    idxl = cg>=2; 
    idxl(1) = 0;
    idx = find(idxl);
    yest = cg(idx-1,1)<2; 
    valvetiming=idx(yest);
    cgg=test{:,j};
    S.Vtick{1,j}=cgg(valvetiming); 
end

S.SumValve=sum(S.V, 2, 'omitnan');
S.MeanValve=mean(S.V, 2, 'omitnan');


S.Running=zeros(50,(length(S.startlevels)-2)); %50 x number of laps
S.Running(:,:)=nan; 
for j=S.lapsused
    cg=S.wsR{1,5}{:,j}; 
    for binnumber=1:50
     S.Running(binnumber, j)=mean(cg(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1),'omitnan');
    end
end

S.MeanRunning=mean(S.Running,2,'omitnan');

S.light=zeros(50,(length(S.startlevels)-2)); %50 x number of laps
S.light(:,:)=nan; 
for j=S.lapsused
    cg=S.wsR{1,6}{:,j}; 
    for binnumber=1:50
        S.light(binnumber, j)=mean(cg(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1),'omitnan');
    end
end

S.MeanLight=mean(S.light,2,'omitnan');

S.R=zeros(50, (length(S.startlevelsR)-2), size(dataset,2));
S.R(:,:,:)=NaN; 

S.Rmax=zeros(50, (length(S.startlevelsR)-2), size(dataset,2));
S.Rmax(:,:,:)=NaN; 

datasetE=zeros(size(S.wsALL,1), 1); datasetE(:,:)=NaN;
datasetE(S.frametiming(:,1),1)=dataset(:,1); datasetER=datasetE(C); 

for j=S.lapsused
    y1=S.startlevelsR(j); y2=(S.startlevelsR(j+1))-1; 
    cg=datasetER(y1:y2,1);
    for binnumber=1:50
        x1=S.threshold(binnumber,j); x2=(S.threshold((binnumber+1),j))-1; 
        S.R(binnumber, j, 1)=mean(cg(x1:x2,1),'omitnan');
    end
end

for i=2:size(dataset,2)
    dataset1=zeros(size(S.wsALL,1), 1); dataset1(:,:)=NaN;
    dataset1S=zeros(size(S.wsALL,1), 1);dataset1S(:,:)=NaN;
    datasetR=dataset1(C);
    datasetRS=dataset1S(C); 
   
    for j=S.lapsused
        cg=datasetR(S.startlevelsR(j):S.startlevelsR(j+1)-1,1);  
        cg2=datasetRS(S.startlevelsR(j):S.startlevelsR(j+1)-1,1); 

        for binnumber=1:50
            S.R(binnumber, j, i)=mean(cg(S.threshold(binnumber,j):S.threshold((binnumber+1),j)-1,1),'omitnan');
            S.Rmax(binnumber, j, i)=max(cg2(S.threshold(binnumber,j):S.threshold((binnumber+1),j)-1,1));            
        end
    end 
end

clearvars -except S dataset*

%%
S.listofneurons = 2:size(S.R,3);

S.threestd=zeros(size(dataset,2),1); S.threestd(:,:)=NaN;
for k=transpose(S.listofneurons)
    histoY=histcounts(dataset(:,k),-20:0.1:10); 
    f = fit(transpose(-19.95:0.1:9.95),transpose(histoY),'gauss1');
    S.threestd(k,1)=3*(f.c1/sqrt(2)); 
    close all
end
 
S.MeanR=mean(S.R, 2, 'omitnan');
S.MeanRmax=mean(S.Rmax, 2, 'omitnan');

S.MeanRmax4=S.MeanRmax; 
S.MeanRshifted=S.MeanRmax; S.MeanRshifted(:,:,:)=NaN; 
S.MeanRmaxshifted=S.MeanRmax; S.MeanRmaxshifted(:,:,:)=NaN; 

S.Rshifted=zeros(50,size(S.R,2),size(S.R,3)); 
S.Rmaxshifted=zeros(50,size(S.R,2),size(S.R,3)); 
S.nshiftedbins=zeros(max(S.listofneurons),1); 
S.pfall=zeros(size(S.R,3),3); S.pfall(:,:)=nan; 
S.middlepf=zeros(size(S.R,3),1); S.middlepf(:,:)=NaN; 

for k=transpose(S.listofneurons)
    S.MeanRmax4(:,:,k)=smooth(S.MeanRmax4(:,:,k),3);
    S.MeanRmax4(:,:,k)=S.MeanRmax4(:,:,k)/max(S.MeanRmax4(:,:,k));    
    
    cgg=S.MeanRmax4(:,:,k); 
    [~, idx]=max(cgg); 
    for i=idx:-1:1
        if cgg(i)<0.2
            S.pfall(k,1)=i+1; 
        break
        end
    end
    
   if isnan(S.pfall(k,1))
    for i=50:-1:1
        if cgg(i)<0.2
            S.pfall(k,1)=i-1;
            break
        end
    end
   end
  
    for i=idx:1:50
        if cgg(i)<0.2
         S.pfall(k,2)=i; 
         break
        end
    end
    
    if isnan(S.pfall(k,2))
        for i=1:1:50
            if cgg(i)<0.2
                S.pfall(k,2)=i;
                break
            end
        end
    end
    
    if S.pfall(k,2)>S.pfall(k,1)
        S.pfall(k,3)=S.pfall(k,1)+round((S.pfall(k,2)-S.pfall(k,1))/2); 
    else
        S.pfall(k,3)=S.pfall(k,1)+(round((S.pfall(k,2)+50-S.pfall(k,1))/2));
        if S.pfall(k,3)>50
            S.pfall(k,3)=S.pfall(k,3)-50;
        end
    end
    
    S.middlepf(k,1)=S.pfall(k,3); 
            
    if S.pfall(k,3)<25
        cg=25-S.pfall(k,3); 
        S.nshiftedbins(k,1)=cg; 
        Ac=reshape(S.R(:,:,k), [], 1);
        add=NaN*ones(cg,1);
        Acnew=vertcat(add, Ac);
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.R,2))
            S.Rshifted(:,j,k)=h{j}; 
        end
        
        Ac=reshape(S.Rmax(:,:,k), [], 1);
        add=NaN*ones(cg,1);
        Acnew=vertcat(add, Ac);
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.Rmax,2))
            S.Rmaxshifted(:,j,k)=h{j}; 
        end
                     
    elseif S.pfall(k,3)>25
        cg=S.pfall(k,3)-25; 
        S.nshiftedbins(k,1)=25-S.pfall(k,3); 
        
        Ac=reshape(S.R(:,:,k), [], 1);
        Ac(1:cg)=[]; add=NaN*ones(cg,1); Acnew=vertcat(Ac, add); 
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.R,2))
            S.Rshifted(:,j,k)=h{j};
        end
        
        Ac=reshape(S.Rmax(:,:,k), [], 1);
        Ac(1:cg)=[]; add=NaN*ones(cg,1); Acnew=vertcat(Ac, add); 
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.Rmax,2))
            S.Rmaxshifted(:,j,k)=h{j};
        end
    else
        S.nshiftedbins(k,1)=0; 
        for j=1:(size(S.R,2))
            S.Rshifted(:,j,k)=S.R(:,j,k);
        end
        
        for j=1:(size(S.Rmax,2))
            S.Rmaxshifted(:,j,k)=S.Rmax(:,j,k);
        end
    end
    
    S.MeanRshifted(:,:,k)=mean(S.Rshifted(:,:,k), 2, 'omitnan');
    S.MeanRmaxshifted(:,:,k)=mean(S.Rmaxshifted(:,:,k), 2, 'omitnan');
       
    S.MeanRmaxshifted4(:,:,k)=S.MeanRmaxshifted(:,:,k);
    S.MeanRmaxshifted4(:,:,k)=smooth(S.MeanRmaxshifted4(:,:,k),3);
    S.MeanRmaxshifted4(:,:,k)=S.MeanRmaxshifted4(:,:,k)/max(S.MeanRmaxshifted4(1:50,:,k),[],1); 
end


S.firstlap=zeros(size(S.Rmax,3),3); S.firstlap(:,:)=nan;
S.pfall=zeros(size(S.Rmax,3),2); S.pfall(:,:)=nan; 
S.maxRmaxshifted=zeros(size(S.Rmax,3), size(S.Rmax,2)); 

for k=transpose(S.listofneurons)
    
    cgg=S.MeanRmaxshifted4(:,:,k);
    [~, idx]=max(cgg); 
    for i=idx:-1:1
        if cgg(i)<0.2
            S.pfall(k,1)=i; 
        break
        end
    end

    for i=idx:1:50
        if cgg(i)<0.2
         S.pfall(k,2)=i; 
         break
        end
    end
    
    if isnan(S.pfall(k,2))
        S.pfall(k,2)=50; 
    end
    
    if isnan(S.pfall(k,1))
        S.pfall(k,1)=1; 
    end
        
    start=S.pfall(k,1); stop=S.pfall(k,2);
    
    prelim=[]; 
    S.maxRmaxshifted(k,:)=max(S.Rmaxshifted(start:stop,:,k),[],1);
    
    if S.pfall(k,1)<3
        S.pfall(k,1)=3;
    end
    
    for i=1:size(S.Rmax,2)
        idxx=find(S.Rmaxshifted(:,i,k)==S.maxRmaxshifted(k,i)); 
        if isempty(idxx)==0 
            for j=idxx-1:-1:S.pfall(k,1)-2
                if S.Rmaxshifted(j,i,k)<S.threestd(k,1)
                    break
                end
            end
            
            if j<S.pfall(k,1)-1
               a=S.pfall(k,1); 
               cg=S.Rmaxshifted(:,i,k); 
               for n=S.pfall(k,1)+1:S.pfall(k,2)
                   if cg(n,1)>S.threestd(k,1)
                       a=n; 
                   else
                       break
                   end
               end
               cg(S.pfall(k,1):a)=NaN;
               S.maxRmaxshifted(k,i)=max(cg(S.pfall(k,1):S.pfall(k,2),1));
            end
        end
    end         
    M=S.maxRmaxshifted(k,:);
    
    B=M>S.threestd(k,1); 
    if (sum(B)>=1)  
    
         laps=find(M>S.threestd(k,1)); 
         llaps=laps<size(S.Rmax,2)-4;
         laps=laps(llaps);
         
         if length(laps)>0
            n=1; 
            for i=laps
                
                if isempty(laps)==1
                    break
                else
               
                    jm2=S.maxRmaxshifted(k,i+1:i+5); 
                    C=jm2>S.threestd(k,1); 
                   if i<size(S.R,2)-31
                        rc=mean(S.R(:,i+1:i+31,k), 2, 'omitnan');
                   else
                        rc=mean(S.R(:,i+1:size(S.R,2),k), 2, 'omitnan');
                   end

                   [rho, pval]=corr(rc, S.MeanR(:,:,k));
                                                           
                    if sum(C)>1 && rho>0 && pval<0.05 
                        prelim(n,1)=i;
                        field=S.maxRmaxshifted(k,i:size(S.Rmax,2));
                        idx=find(field>S.threestd(k,1)); 
                        idx2=zeros(1,length(idx)-1); 
                        for m=1:length(idx)-1
                           idx2(1,m)=idx(1,m+1)-idx(1,m);  
                        end
                        idx3=find(idx2>21);  
                        if isempty(idx3)==1
                            S.firstlap(k,1)=i;
                            break
                        else
                            b=laps>prelim(n,1)-2+idx(1,idx3(1,1)+1); 
                            laps=laps(b);
                        end
                        n=n+1; 
                    end   
                end
    
            end
            
        end
    end
    
    if isnan(S.firstlap(k,1))
         if isempty(prelim)==0
             S.firstlap(k,1)=prelim(1,1); 
         end
    end
end

for k=transpose(S.listofneurons)
    if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)
        maxLocR=find(S.Rmaxshifted(:,S.firstlap(k,1),k)==S.maxRmaxshifted(k,S.firstlap(k,1))); 
        if maxLocR-S.nshiftedbins(k,1)<1
            S.firstlap(k,2:3)=S.firstlap(k,1)-1; 
                if S.firstlap(k,1)-1<1
                    S.firstlap(k,2:3)=1; 
                end
        elseif maxLocR-S.nshiftedbins(k,1)>50
            S.firstlap(k,2:3)=S.firstlap(k,1)+1; 
        else
            S.firstlap(k,2:3)=S.firstlap(k,1);
        end
    end
end

clearvars -except S dataset*

%%
idx=size(S.wsALL,2)+1; %idx should be 10 
S.wsALL(:,idx)=S.wsALL(:,5); %%%binary frametiming;
S.wsALL(:,idx)=0; S.wsALL(S.frametiming(:,1),idx)=5; 
C=S.wsALL(:,2)>0.05; 
S.running.allframesR=S.wsALL(C,idx); 
S.frametimingR=find(S.running.allframesR>4); 

MeanRnew=zeros(50,1,size(S.R,3)); 
S.MeanRnew4=zeros(50,1,size(S.R,3)); %peak-normalized

for k=transpose(S.listofneurons)
    if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)
        MeanRnew(:,:,k)=mean(S.R(:,S.firstlap(k,2):size(S.R,2),k), 2, 'omitnan');
        S.MeanRnew4(:,:,k)=smooth(MeanRnew(:,:,k),3);
        S.MeanRnew4(:,:,k)=S.MeanRnew4(:,:,k)/max(S.MeanRnew4(:,:,k),[],1); 
     end
end

[~, I]=max(S.MeanRnew4); 

S.Rnewshifted=zeros(50,size(S.R,2),size(S.R,3)); 
S.nshiftedbinsNew=zeros(max(S.listofneurons),2); 
for k=transpose(S.listofneurons)
    if I(:,:,k)<25
        cg=25-I(:,:,k); 
        S.nshiftedbinsNew(k,1)=cg; 
        Ac=reshape(S.R(:,:,k), [], 1);
        add=NaN*ones(cg,1);
        Acnew=vertcat(add, Ac);
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.R,2))
            S.Rnewshifted(:,j,k)=h{j}; 
        end
    elseif I(:,:,k)>25
        cg=I(:,:,k)-25; 
        S.nshiftedbinsNew(k,1)=cg; 
        Ac=reshape(S.R(:,:,k), [], 1);
        Ac(1:cg)=[]; add=NaN*ones(cg,1); Acnew=vertcat(Ac, add); 
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.R,2))
            S.Rnewshifted(:,j,k)=h{j};
        end
    else
        S.nshiftedbinsNew(k,1)=0; 
        for j=1:(size(S.R,2))
            S.Rnewshifted(:,j,k)=S.R(:,j,k);
        end
    end
end

%for cells - calculate new averages of shifted data based on first lap &
%calculate SI for actual data

S.MeanRnewshifted=zeros(50,1,size(S.R,3)); 
S.SI=zeros(50,1,size(S.R,3)); 
for k=transpose(S.listofneurons)
    if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)-3
        S.MeanRnewshifted(:,:,k)=mean(S.Rnewshifted(:,S.firstlap(k,1):size(S.R,2),k), 2, 'omitnan');
        data=S.MeanRnewshifted(:,:,k);
        for i=1:50
            S.SI(i,1,k) = 1/50 * data(i,1) * log2(data(i,1)/mean(data));
        end
    end
end

S.sumSI=sum(S.SI); 

%calculate shuffle
%shuffle(number of spatial bins, cells, number of shuffles)

S.shuffle.shuffle=zeros(50,size(S.R,3),100); S.shuffle.shuffle(:,:,:)=NaN; 
S.frametimingSS=[];
MeanRshuffle=zeros(50,1,size(S.R,3)); MeanRshuffle(:,:,:)=NaN; 

for i=transpose(S.listofneurons)
    if S.firstlap(i,1)>0 && S.firstlap(i,1)<size(S.R,2)-3
        startpoint=find(S.frametimingR>S.startlevelsR(S.firstlap(i,2),1),1);    
        startFrame=S.frametimingR(1:startpoint-1,1);
       
        for shufflenumber=1:100
            shufflenumber;      
            p = randperm(6);
            S.shufflep{1,shufflenumber,i}=p;
            x=randi(1000); 
            
            Rshuffle=zeros(50, size(S.R,2), size(S.R,3)); Rshuffle(:,:,:)=NaN;
            cg2=[]; cg2=S.frametimingR(startpoint:length(S.frametimingR),1);
            frametimingS=[]; 
            frametimingS(:,1)=vertcat(cg2(x:length(cg2),1),cg2(1:(x-1),1));
            
            part=floor(length(frametimingS)/6); 

            se{1,1}=frametimingS(1:part,1);
            se{1,2}=frametimingS((part+1):part*2,1);
            se{1,3}=frametimingS((part*2)+1:part*3,1);
            se{1,4}=frametimingS((part*3)+1:part*4,1);
            se{1,5}=frametimingS((part*4)+1:part*5,1);
            se{1,6}=frametimingS((part*5)+1:length(frametimingS),1);
            
            S.frametimingSS(:, shufflenumber,i)= vertcat(startFrame, se{1,p(1)}, se{1,p(2)}, se{1,p(3)}, se{1,p(4)}, se{1,p(5)}, se{1,p(6)});  
                                                          
            datasetS=zeros(size(S.wsALL,1), 1); datasetS(:,:)=NaN;                 
            datasetS=datasetS(C); 
            datasetRS=zeros(length(S.running.allframesR),1); datasetRS(:,:)=NaN; 
            
            if length(datasetS)>length(S.frametimingR) || length(datasetS)<length(S.frametimingR)
                errorMessage = sprintf('Problem');
                uiwait(warndlg(errorMessage));
                return;  
            end
            
            datasetRS(S.frametimingSS(:,shufflenumber,i),1)=datasetS(:,1); 
            
           for j=S.firstlap(i,2):size(S.R,2)
                cg=datasetRS(S.startlevelsR(j):S.startlevelsR(j+1)-1,1);  
                for binnumber=1:50
                    Rshuffle(binnumber, j, i)=mean(cg(S.threshold(binnumber,j):S.threshold((binnumber+1),j)-1,1),'omitnan');                        
                end
            end
                                                                                                                                                        
            MeanRshuffle(:,:,i)=mean(Rshuffle(:,S.firstlap(i,2):size(S.R,2),i), 2, 'omitnan');
            S.shuffle.shuffle(:,i,shufflenumber)=MeanRshuffle(:,:,i);
       end
    end
   
end
   
S.shuffle.shuffleSI=zeros(50,size(S.R,3),100); 
S.shuffle.sumshuffleSI=zeros(1, size(S.R,3), 100);
S.shuffle.shuffleCI=zeros(1, size(S.R,3)); 

for k=transpose(S.listofneurons)
   if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)-3
    for shufflenumber=1:100
        data=S.shuffle.shuffle(:,k,shufflenumber); 
         for i=1:50
             S.shuffle.shuffleSI(i,k,shufflenumber) = 1/50 * data(i,1) * log2(data(i,1)/mean(data));
         end 
         S.shuffle.sumshuffleSI=sum(S.shuffle.shuffleSI); 
    end
    jcm=permute(S.shuffle.sumshuffleSI(1,k,:), [3 2 1]); 

    ts=tinv([0.025  0.975],length(jcm)-1); %valid for normal distribution
    S.shuffle.shuffleCI(1,k)=mean(jcm)+ts(1,2)*std(jcm)/sqrt(length(jcm)); 
   end
   
end

S.sumSI=permute(S.sumSI, [3 2 1]); 
S.shuffle.shuffleCI=permute(S.shuffle.shuffleCI, [2 1]); 
thres=S.sumSI>S.shuffle.shuffleCI; 
interim=transpose(1:size(dataset,2)); placecells=interim(thres); 

%activity in more than 1/3 of laps 
active5=zeros(size(S.R,3),1); 
MeanRnew=zeros(50,1,size(S.R,3)); MeanRnew(:,:,:)=NaN; 
MeanRmaxshifted=MeanRnew; MeanRmaxshifted(:,:,:)=NaN; 
pfall=zeros(size(S.R,3),2); pfall(:,:)=NaN; 
borders=zeros(max(S.listofneurons),2);

for k=transpose(placecells)
    
    MeanRnew(:,:,k)=mean(S.Rmax(:,S.firstlap(k,2):size(S.R,2),k), 2, 'omitnan');
    MeanRnew(:,:,k)=smooth(MeanRnew(:,:,k),3);
    MeanRnew(:,:,k)= MeanRnew(:,:,k)/max(MeanRnew(:,:,k)); 
    
    cgg=MeanRnew(:,:,k); 
    [~, idx]=max(cgg); 
    for i=idx:-1:1
        if cgg(i)<0.2
            pfall(k,1)=i; 
        break
        end
    end
    
   if isnan(pfall(k,1))
    for i=50:-1:1
        if cgg(i)<0.2
            pfall(k,1)=i;
            break
        end
    end
   end
  
    for i=idx:1:50
        if cgg(i)<0.2
         pfall(k,2)=i; 
         break
        end
    end
    
    if isnan(pfall(k,2))
        for i=1:1:50
            if cgg(i)<0.2
                pfall(k,2)=i;
                break
            end
        end
    end
    
    if pfall(k,2)>pfall(k,1)
        pfall(k,3)=pfall(k,1)+round((pfall(k,2)-pfall(k,1))/2); 
    else
        pfall(k,3)=pfall(k,1)+(round((pfall(k,2)+50-pfall(k,1))/2));
        if pfall(k,3)>50
            pfall(k,3)=pfall(k,3)-50;
        end
    end
            
    if pfall(k,3)<25
        cg=25-pfall(k,3); 
        Ac=reshape(S.Rmax(:,:,k), [], 1);
        add=NaN*ones(cg,1);
        Acnew=vertcat(add, Ac);
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.Rmax,2))
            Rmaxshifted(:,j,k)=h{j}; 
        end
    end
                     
    if pfall(k,3)>25       
        cg=pfall(k,3)-25; 
        Ac=reshape(S.Rmax(:,:,k), [], 1);
        Ac(1:cg)=[]; add=NaN*ones(cg,1); Acnew=vertcat(Ac, add); 
        n=numel(Acnew);
        h=mat2cell(Acnew,diff([0:50:n-1,n]));
        for j=1:(size(S.Rmax,2))
            Rmaxshifted(:,j,k)=h{j};
        end
    end
    
    if pfall(k,3)==25
        for j=1:(size(S.Rmax,2))
            Rmaxshifted(:,j,k)=S.Rmax(:,j,k);
        end
    end
       
    MeanRmaxshifted(:,:,k)=mean(Rmaxshifted(:,S.firstlap(k,1):size(S.R,2),k), 2, 'omitnan');
    MeanRmaxshifted(:,:,k)=smooth(MeanRmaxshifted(:,:,k),3);
    MeanRmaxshifted(:,:,k)=MeanRmaxshifted(:,:,k)/max(MeanRmaxshifted(:,:,k),[],1); 
    
    cgg=MeanRmaxshifted(:,:,k);
    [~, idx]=max(cgg); 
    
    for i=idx:-1:1 
        if cgg(i)<0.2
            borders(k,1)=i+1; 
        break
        end
    end

    for i=idx:1:50 
        if cgg(i)<0.2
            borders(k,2)=i-1; 
         break
        end
    end
    
    if borders(k,2)==0
        borders(k,2)=50; 
    end
    
    if borders(k,1)==0 
        borders(k,1)=1; 
    end
             
    subR=Rmaxshifted(borders(k,1):borders(k,2),S.firstlap(k,1):size(S.R,2),k);
    b=size(subR); 
    active=subR>S.threestd(k,1); 
    active2=sum(active);
    active3=active2>0; 
    active4=sum(active3)>floor(b(1,2)/3)-1; %reliability criterium
    active5(k,1)=active4; 
end

S.placecells=interim(active5>0);

MeanRnew=zeros(50,1,size(S.R,3)); 
S.MeanRnew4=zeros(50,1,size(S.R,3)); %peak-normalized

for k=transpose(S.listofneurons)
    if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)
        MeanRnew(:,:,k)=mean(S.R(:,S.firstlap(k,2):size(S.R,2),k), 2, 'omitnan');
        MeanRnew(:,:,k)=smooth(MeanRnew(:,:,k),3);
        S.MeanRnew4(:,:,k)=MeanRnew(:,:,k)/max(MeanRnew(:,:,k),[],1); 
     end
end

clearvars -except S dataset* 

%%

%for place field peak shift - calculate new averages of shifted data based on 
%firstlap(k,1)+1:firstlap(k,1)+31
%makes gaussian fitting easier
MeanRCOM=zeros(50,1,size(S.R,3)); 
S.MeanRCOM3=zeros(50,1,size(S.R,3)); 

for k=transpose(S.placecells)
    if S.firstlap(k,1)>0 && S.firstlap(k,1)<size(S.R,2)-31
        MeanRCOM(:,:,k)=mean(S.Rshifted(:,S.firstlap(k,1)+1:S.firstlap(k,1)+31,k), 2, 'omitnan');
    else
        MeanRCOM(:,:,k)=mean(S.Rshifted(:,S.firstlap(k,1)+1:size(S.R,2),k), 2, 'omitnan');
    end
    
    S.MeanRCOM3(:,:,k)=smooth(MeanRCOM(:,:,k),3); 
end


%place field peak shift

S.COM=zeros(size(S.Rshifted,3), 3); S.COM(:,:)=NaN; 
for k=transpose(S.placecells)
    [~, I]=max(S.MeanRCOM3(:,:,k)); 
    I=I(1); 
    if I<25
        cg=25-I; 
        cgg=vertcat(S.MeanRCOM3(50-cg:50,:,k), S.MeanRCOM3(1:49-cg,:,k));
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);                                                                                                                                                             
        S.COM(k,2)=f.b1-cg;
        if (S.COM(k,2)<1)
            S.COM(k,2)=S.COM(k,2)+50;
        end
    elseif I>25
        cg=I-25; 
        cgg=vertcat(S.MeanRCOM3(cg:50,:,k), S.MeanRCOM3(1:cg-1,:,k));
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);
        S.COM(k,2)=f.b1+cg;
        if (S.COM(k,2)>50)
            S.COM(k,2)=S.COM(k,2)-50;
        end
    else
        cgg=S.MeanRCOM3(:,:,k);
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);
        S.COM(k,2)=f.b1;
    end
   
    I=find(S.Rmaxshifted(:,S.firstlap(k,1),k)==S.maxRmaxshifted(k,S.firstlap(k,1))); 
    I=I(1); 
    if I<25
        cg=25-I; 
        cgg=vertcat(S.Rmaxshifted(50-cg:50,S.firstlap(k,1),k), S.Rmaxshifted(1:49-cg,S.firstlap(k,1),k));
        if isnan(cgg(1))
            cgg(1)=0; 
        end
        if isnan(cgg(50))
            cgg(50)=0; 
        end
        nanx = isnan(cgg);
        t = 1:numel(cgg);
        cgg(nanx) = interp1(t(~nanx), cgg(~nanx), t(nanx));
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);
        S.COM(k,1)=f.b1-cg;
        if (S.COM(k,1)<1)
            S.COM(k,1)=S.COM(k,1)+50;
        end
    elseif I>25
        cg=I-25; 
        cgg=vertcat(S.Rmaxshifted(cg:50,S.firstlap(k,1),k), S.Rmaxshifted(1:cg-1,S.firstlap(k,1),k));
        if isnan(cgg(1))
            cgg(1)=0; 
        end
        if isnan(cgg(50))
            cgg(50)=0; 
        end
        nanx = isnan(cgg);
        t = 1:numel(cgg);
        cgg(nanx) = interp1(t(~nanx), cgg(~nanx), t(nanx));
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);
        S.COM(k,1)=f.b1+cg;
        if (S.COM(k,1)>50)
            S.COM(k,1)=S.COM(k,1)-50;
        end
    else
        cgg=S.Rmaxshifted(:,S.firstlap(k,1),k);
        if isnan(cgg(1))
            cgg(1)=0; 
        end
        if isnan(cgg(50))
            cgg(50)=0; 
        end
        nanx = isnan(cgg);
        t = 1:numel(cgg);
        cgg(nanx) = interp1(t(~nanx), cgg(~nanx), t(nanx));
        f = fit(transpose(1:50),cgg,'gauss1','Lower', [max(cgg)/2 15 0],'Upper', [max(cgg) 35 25]);
        S.COM(k,1)=f.b1;
    end
end

%pf width as a function of running speed
%calculate new averages of non-shifted data based on 
%firstlap(k,3)+1:firstlap(k,3)+31

MeanRwid=zeros(50,1,size(S.R,3)); 
S.MeanRwid4=zeros(50,1,size(S.R,3));  %peak-scaled

S.pf=zeros(size(S.R,3),2); S.pf(:,:)=nan;  
for k=transpose(S.placecells)
    if S.firstlap(k,3)>0 && S.firstlap(k,3)<size(S.R,2)-31
        MeanRwid(:,:,k)=mean(S.R(:,S.firstlap(k,3)+1:S.firstlap(k,3)+31,k), 2, 'omitnan');
    else
        MeanRwid(:,:,k)=mean(S.R(:,S.firstlap(k,3)+1:size(S.R,2),k), 2, 'omitnan');
    end
    
    MeanRwid(:,:,k)=smooth(MeanRwid(:,:,k),3); 
    S.MeanRwid4(:,:,k)=MeanRwid(:,:,k)/max(MeanRwid(:,:,k)); 
    
    [~, I]=max(S.MeanRwid4(:,:,k)); 
    I=I(1); 
    if I<25
        cg=25-I; 
        cgg=vertcat(S.MeanRwid4(50-cg:50,:,k), S.MeanRwid4(1:49-cg,:,k));
        [~, idx]=max(cgg); 
        for i=idx:-1:1
            if cgg(i)<0.2
                 S.pf(k,1)=i-cg;
                 if S.pf(k,1)<1
                     S.pf(k,1)=S.pf(k,1)+50; 
                 end
                 break
            end
        end
        
        if isnan(S.pf(k,1))
            for i=50:-1:1
                if cgg(i)<0.2
                    S.pf(k,1)=i-cg;
                    if S.pf(k,1)<1
                        S.pf(k,1)=S.pf(k,2)+50; 
                    end
                    break
                end
            end
        end
                
        for i=idx:1:50
            if cgg(i)<0.2
                 S.pf(k,2)=i-2-cg;
                 if S.pf(k,2)<1
                     S.pf(k,2)=S.pf(k,2)+50; 
                 end
                 break
            end
        end
        
        if isnan(S.pf(k,2))
            for i=1:1:50
                if cgg(i)<0.2
                    S.pf(k,2)=i-2-cg;
                    if S.pf(k,2)<1
                        S.pf(k,2)=S.pf(k,2)+50; 
                    end
                    break
                end
            end
        end
        
    elseif I>25
        cg=I-25; 
        cgg=vertcat(S.MeanRwid4(cg:50,:,k), S.MeanRwid4(1:cg-1,:,k));
        [~, idx]=max(cgg); 
        for i=idx:-1:1
            if cgg(i)<0.2
                 S.pf(k,1)=i+cg; 
                 if S.pf(k,1)>50
                     S.pf(k,1)=S.pf(k,1)-50; 
                 end
                 break
            end         
        end
        
       if isnan(S.pf(k,1))
            for i=50:-1:1
                if cgg(i)<0.2
                    S.pf(k,1)=i+cg;
                    if S.pf(k,1)>50
                        S.pf(k,1)=S.pf(k,1)-50; 
                    end
                    break
                end
            end
        end
               
        for i=idx:1:50
            if cgg(i)<0.2
                 S.pf(k,2)=i-2+cg;
                 if S.pf(k,2)>50
                     S.pf(k,2)=S.pf(k,2)-50; 
                 end
                 break
            end
        end 
        
        if isnan(S.pf(k,2))
            for i=1:1:50
                if cgg(i)<0.2
                    S.pf(k,2)=i-2+cg;
                    if S.pf(k,2)>50
                        S.pf(k,2)=S.pf(k,2)-50; 
                    end
                    break
                end
            end
        end
                 
    else
        cgg=S.MeanRwid4(:,:,k);
        [~, idx]=max(cgg); 
        for i=idx:-1:1
            if cgg(i)<0.2
                 S.pf(k,1)=i+1; 
                 break
            end
        end
        
        for i=idx:1:50
            if cgg(i)<0.2
                 S.pf(k,2)=i-1; 
                 break
            end
        end 
                
    end  
end

%check whether pf size and pf peak shift were correctly detected


%% calculate pf width
 
binsize=zeros(length(S.startlevels),1); 
for i=1:(length(S.startlevels)-1)
    binsize(i,1)=max(S.wsL{1,8}{1,i})/50; 
end

S.pf(:,3)=NaN; S.pf(:,4)=NaN; %pf(:,3)=width in bins; pf(:,4)=runningspeed in V
for k=transpose(S.placecells)
    if S.pf(k,2)>S.pf(k,1)
        S.pf(k,3)=S.pf(k,2)-S.pf(k,1)+1; 
    else
        S.pf(k,3)=51-S.pf(k,1)+S.pf(k,2);
    end
    
    if S.firstlap(k,3)==1
        if S.pf(k,2)>S.pf(k,1)
            if S.pf(k,1)==1
                start=1;
            else
                startV=((S.pf(k,1)-1)*binsize(S.firstlap(k,3)));
                idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>startV;
                start = find(idxl,1); 
            end
            
            if S.pf(k,2)==50
                term=length(S.wsL{1,8}{1,S.firstlap(k,3)});
            else
                termV=S.pf(k,2)*binsize(S.firstlap(k,3));
                idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>=termV; 
                term = (find(idxl,1))-1; 
            end
            
            running=S.wsL{1,2}{1,S.firstlap(k,3)}; 
            C=running<-0.005; running(C)=NaN; 
            S.pf(k,4)=mean(running(start:term, 1),1,'omitnan');
%         else %only for earlier experiments where i do not add startmarker at the beginning;
%             %for later experiments (pt1==1) this data does not exist
%             %uncomment only if S.startmarkerChange==0
%             cd (S.parentfolder)
%             data=ws.loadDataFile(S.fileNames{1});
%             names=fieldnames(data);
%             s=names{2};
%             dataanalog=data.(s);
%             ws1=dataanalog.analogScans;
%             idxl = ws1(:,1)>=2;
%             idxl(1) = 0;
%             idx = find(idxl);
%             yest = ws1(idx-1,1)<2; 
%             startlevels1=idx(yest)-1; 
%             running=ws1(1:startlevels1(1),2);
%             pos0=cumtrapz(running);
%             pos1=S.wsL{1,8}{1,1};
%             postest=pos0+(max(pos1)-max(pos0));
%             startV=(S.pf(k,1)-1)*binsize(1); 
%             termV=S.pf(k,2)*binsize(S.firstlap(k,3)); 
%             idxl = postest>startV;
%             start = find(idxl,1); 
%             if start==1
%                S.pf(k,4)=NaN;
%             else
%                 idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>=termV; 
%                 term = (find(idxl,1))-1; 
%                 C=running<-0.005; running(C)=NaN;
%                 pt1=mean(running(start:length(running), 1),1,'omitnan');
%                 running=S.wsL{1,2}{1,S.firstlap(k,3)}; 
%                 C=running<-0.005; running(C)=NaN;
%                 pt2=mean(running(1:term, 1),1,'omitnan');
%                 S.pf(k,4)=(pt1+pt2)/2;
%             end
         end
    else
        if S.pf(k,2)>S.pf(k,1)
            if S.pf(k,1)==1
                start=1;
            else
                startV=((S.pf(k,1)-1)*binsize(S.firstlap(k,3),1));
                idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>startV;
                start = find(idxl,1); 
            end
            
            if S.pf(k,2)==50
                term=length(S.wsL{1,8}{1,S.firstlap(k,3)});
            else
                termV=S.pf(k,2)*binsize(S.firstlap(k,3));
                idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>=termV; 
                term = (find(idxl,1))-1; 
            end
            
            running=S.wsL{1,2}{1,S.firstlap(k,3)}; 
            C=running<-0.005; running(C)=NaN; 
            S.pf(k,4)=mean(running(start:term, 1),1,'omitnan');
        else
            startV=(S.pf(k,1)-1)*binsize((S.firstlap(k,3)-1)); 
            termV=S.pf(k,2)*binsize(S.firstlap(k,3)); 
            idxl = S.wsL{1,8}{1,S.firstlap(k,3)-1}>startV;
            start = find(idxl,1); 
            idxl = S.wsL{1,8}{1,S.firstlap(k,3)}>=termV; 
            term = (find(idxl,1))-1; 
            running=S.wsL{1,2}{1,S.firstlap(k,3)-1}; 
            C=running<-0.005; running(C)=NaN;
            pt1=mean(running(start:length(running), 1),1,'omitnan');
            running=S.wsL{1,2}{1,S.firstlap(k,3)}; 
            C=running<-0.005; running(C)=NaN;
            pt2=mean(running(1:term, 1),1,'omitnan');
            S.pf(k,4)=(pt1+pt2)/2; 
        end
    end
end

S.pf(:,5)=S.pf(:,3).*3.6; %pf(:,5)=width in cm; 
S.pf(:,6)=S.pf(:,4).*40;% running speed in cm/s; 

S.COM(:,3)=S.COM(:,2)-S.COM(:,1); %negative values mean left shift in the colorplot(=advancement);

S.nlaps=size(S.R,2); 
S.reclengthSec=length(S.wsALL(:,1))/S.acq; 

clearvars -except S dataset*

