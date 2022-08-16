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

%% combine reg tiffs from suite2p%%
cd (S.parentfolder)
cd suite2p/plane0/reg_tif/

datasetfield=zeros(length(dir('*.tif'))*50,2); datasetfield(:,:)=NaN; 
j=1; m=1;
for i=0:(length(dir('*.tif'))-1)
    vol=ScanImageTiffReader(['file' num2str(i,'%03i') '_chan0.tif']).data();
    vol2=rot90(fliplr(vol),1);
    for k=1:size(vol2,3)        
        datasetfield(j,1)=mean(vol2(100:140,200:240,k),'all','omitnan');        
        datasetfield(j,2)=mean(vol2(1+75:512-75,1+75:512-75,k),'all','omitnan');
        j=j+1;
    end
    S.meanImg(:,:,m)=mean(vol2,3,'omitnan'); 
    m=m+1; 
end

S.meanImage=mean(S.meanImg,3,'omitnan'); 
cd (S.parentfolder)
saveastiff(int16(S.meanImage), 'all.tif');

clearvars -except S datasetfield control;

%% noisecorrelation

%load Fall.mat

S.lapsused=1:(length(S.startlevels)-2);
S.medfiltC=500;                                                                                                                                    C=500;
S.polynom=1; 
S.noiseCorrThres=0.5; %usually 0.5 to start
S.lenfile=400; %usualy 400

S.bar=S.bar2; S.bar(:,:)=-2;
S.F=F;
S.Fneu=Fneu;
S.iscell=iscell;
S.ops=ops;
S.spks=spks;
S.stat=stat;
S.datasetfieldorig=datasetfield;
datasetALL=transpose(F); 

datasetALL=horzcat(datasetfield, datasetALL);
S.datasetALLorig=datasetALL; 

S.allROIs=vertcat(1,2,find(iscell(:,1)==1)+2); %ROI1+2 for data from datasetfield 

S.FzeroALL=zeros(size(datasetALL,2),1); S.FzeroALL(:,:)=NaN;

startframe=1:S.lenfile:size(datasetALL,1);

 for k=transpose(S.allROIs)
             
     test=datasetALL(:,k); 
    
     if k==1
         thres=test<S.bar2(1); 
     end
    
    y=test; y(thres)=NaN; 
    
    test2=y(~thres); 
    
    edges=-150:0.1:3000;
    histoX=-149.95:0.1:2999.95;
    ha=histogram(test2,edges); hold on; 
    histoY=ha.Values; 
    f = fit(transpose(histoX),transpose(histoY),'gauss1'); 
    S.FzeroALL(k,1)=f.b1; 
   
    x1=[S.FzeroALL(k,1) S.FzeroALL(k,1)]; y1=[0 15];  
    hold on; plot(x1,y1, 'r-', 'LineWidth',2); title(k);
    y=(y-S.FzeroALL(k,1))/(S.FzeroALL(k,1)); 
    datasetALL(:,k)=y; 
    
    x=1:size(datasetALL,1); x=transpose(x); 
    y=datasetALL(:,k); 
    
    y(thres)=NaN; 
    filtdata=medfilt1(y,S.medfiltC, 'omitnan'); %used to be 500
    
    c = polyfit(x,filtdata,S.polynom);
    y_est = polyval(c,x);
    datasetALL(:,k)=datasetALL(:,k)-y_est; 
    
    edges=-1000:0.02:1000; 
    histoX=-999.99:0.02:999.99;
    histoY=histcounts(datasetALL(:,k),edges); 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    datasetALL(:,k)=datasetALL(:,k)-f.b1;
    close all; 
    
    y=datasetALL(:,k); 
    
    y(thres)=NaN; 
    for i=1:length(startframe)
        y(startframe(i):startframe(i)+1,1)=NaN;
    end
    
    datasetALL(:,k)=y;
 end

S.datasetALL=datasetALL;

S.thresholdS=zeros(50,length(S.startlevels)-2); % must have 51 rows at the end
 
for i=S.lapsused 
    cg=S.wsL{1,8}{1,i};
    for binnumber=1:50
        thres=binnumber * (S.maxpos(i,1)/50); 
        idxl = cg>=thres;
        idxl(1) = 0;
        idx = find(idxl, 1);
        yest = cg(idx-1,1)<thres; 
        if idx(yest)>0
            S.thresholdS(binnumber,i)=idx(yest)-1;
        else
            S.thresholdS(binnumber,i)=numel(cg);
        end
    end
end

interim=zeros(1,(length(S.startlevels))-2); interim(:,:)=1;
S.thresholdS=vertcat(interim, S.thresholdS);

S.Rall=zeros(50, (length(S.startlevels)-2), size(datasetALL,2));
S.Rall(:,:,:)=NaN; 

datasetE=zeros(size(S.wsALL,1), 1);
B=datasetE<1;
datasetE(B==1)=nan;
for k=1:size(datasetALL,1)
    datasetE(S.frametiming(k),1)=datasetALL(k,1); 
end
    
%all data
datasetER=datasetE; 

for j=S.lapsused
    y1=S.startlevels(j); y2=(S.startlevels(j+1))-1; 
    cg=datasetER(y1:y2,1); 
    for binnumber=1:50
        x1=S.thresholdS(binnumber,j); x2=(S.thresholdS((binnumber+1),j))-1; 
        cgg=cg(x1:x2,1); 
        S.Rall(binnumber, j, 1)=mean(cgg,'omitnan');
    end
end

for i=transpose(S.allROIs)%%cellnumber
    dataset1=zeros(size(S.wsALL,1), 1);
    B=dataset1<1;
    dataset1(B==1)=nan;
    
    for k=1:size(datasetALL,1)
        dataset1(S.frametiming(k),1)=datasetALL(k,i); 
    end
    
    %all data
    datasetR=dataset1;

    for j=S.lapsused
        cg=datasetR(S.startlevels(j):S.startlevels(j+1)-1,1);  
        for binnumber=1:50
            S.Rall(binnumber, j, i)=mean(cg(S.thresholdS(binnumber,j):(S.thresholdS((binnumber+1),j))-1,1),'omitnan');           
        end
        
     end 
end

S.Ralllaps=S.Rall; 
i=2; k=1;
goodlaps=[];
for j=S.lapsused
    b=find(isnan(S.Rall(:,j,i)));
    if isempty(b)==1
        goodlaps(k,1)=j;
        k=k+1;
    end
end

S.Rall=S.Rall(:,goodlaps,:); 
S.MeanRall=mean(S.Rall, 2, 'omitnan');
S.RallSubtract=S.Rall; 
subtracteddata=zeros(binnumber*size(S.Rall,2),size(S.Rall,3)); subtracteddata(:,:)=NaN;  
for i=transpose(S.allROIs)
    for j=1:size(S.RallSubtract,2)
        S.RallSubtract(:,j,i)=S.Rall(:,j,i)-S.MeanRall(:,:,i); 
    end
    subtracteddata(:,i)=reshape(S.RallSubtract(:,:,i), [], 1);
end

tbin=1;
NN=size(subtracteddata,2); 
NT=size(subtracteddata,1); 
nt=(floor(NT)/tbin)*tbin;
nt2=reshape(subtracteddata, [tbin,nt/tbin,NN]);
nt2=mean(nt2,1,'omitnan');
nt2=permute(nt2,[2,3,1]); 
nt3=zscore(nt2,1,'omitnan'); 
S.noisecorr=(transpose(nt3)*nt3)/(nt/tbin);

% make histogram of noise correlations
S.noiseCorrValues=zeros(1,1); 
i=1;
cg=cell(1,1); 
for k=1:size(S.noisecorr,2)-1
    cg{1,i}=S.noisecorr(k,k+1:size(S.noisecorr,2));     
    i=i+1;
end

S.noiseCorrValues=[];
for k=1:length(cg)
    S.noiseCorrValues=vertcat(S.noiseCorrValues, transpose((cg{1,k}))); 
end

S.ROIlist=cell(1,1);
noisecorrdel=S.noisecorr; 
n=1;
for i=transpose(S.allROIs(3:length(S.allROIs)))
    m=1; 
    test=[]; 
    for j=3:length(S.noisecorr)
        if noisecorrdel(i,j)>S.noiseCorrThres
            test(:,m)=j;
            m=m+1;
        end
    end
    
    if isempty(test)==0
        S.ROIlist{1,n}=test;
        n=n+1; 
        for j=1:length(test)
            noisecorrdel(test(j),:)=NaN;
            noisecorrdel(:,test(j))=NaN;
        end
    end
end

clearvars -except S control datasetfield

%%

%use only after having computed datasetfield and S.ROIlist

S.medfiltSm=5;

S.ROIlist{1,(length(S.ROIlist)+1)}=horzcat(S.ROIlist{1,3:length(S.ROIlist)});

S.listofneurons=2:length(S.ROIlist); 
S.listofneurons=transpose(S.listofneurons);

S.datasetfieldorig=datasetfield;
S.F=F;
S.Fneu=Fneu;
S.iscell=iscell;
S.ops=ops;
S.spks=spks;
S.stat=stat;

datasetALL=transpose(F); 
datasetALL=horzcat(datasetfield, datasetALL);

S.datasetALLorig=datasetALL; 
S.npixels=S.ROIlist;
S.sumnpixels=zeros(length(S.ROIlist),2); S.sumnpixels(:,:)=NaN;
dataset=zeros(length(datasetALL), length(S.ROIlist)); dataset(:,:)=NaN;
for k=1:length(S.ROIlist)
        if length(S.ROIlist{1,k})>1
            S.npixels{1,k}=single(S.ROIlist{1,k}); 
            cg=zeros(length(datasetALL), length(S.ROIlist{1,k})); 
            for i=1:length(S.ROIlist{1,k})
                S.npixels{1,k}(i)=S.stat{1,S.ROIlist{1,k}(i)-2}.npix; 
                cg(:,i)=datasetALL(:, S.ROIlist{1,k}(i)); 
                cg(:,i)=S.npixels{1,k}(i).*cg(:,i);
            end
            dataset(:,k)=sum(cg,2)/sum(S.npixels{1,k});
        else
            dataset(:,k)=datasetALL(:, S.ROIlist{1,k});
            if k>2
                i=1;
                S.npixels{1,k}(i)=S.stat{1,S.ROIlist{1,k}(i)-2}.npix; 
            end
        end
    
    S.sumnpixels(k,1)=sum(S.npixels{1,k});
end

S.Fzero=zeros(size(dataset,2),1); S.Fzero(:,:)=NaN;

for k=1:size(dataset,2)
    
    test=dataset(:,k); 
    
    if k==1
        thres=test<S.bar2(1); 
    end
    
    y=test; y(thres)=NaN; 

    test2=y(~thres); 
    
    edges=-150:0.1:3000;
    histoX=-149.95:0.1:2999.95;
    ha=histogram(test2,edges); hold on; 
    histoY=ha.Values; 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    S.Fzero(k,1)=f.b1; 

    x1=[S.Fzero(k,1) S.Fzero(k,1)]; y1=[0 15];  
    hold on; plot(x1,y1, 'r-', 'LineWidth',2); title(k);
    y=(y-S.Fzero(k,1))/(S.Fzero(k,1)); 
    dataset(:,k)=y; 
    
end

for k=1:size(dataset,2)
    x=1:size(dataset,1); x=transpose(x); 
    y=dataset(:,k); 
    
    y(thres)=NaN; 
    filtdata=medfilt1(y,S.medfiltC, 'omitnan'); 
    
    c = polyfit(x,filtdata,S.polynom);
    y_est = polyval(c,x);
    dataset(:,k)=dataset(:,k)-y_est; 
    
    edges=-1000:0.02:1000; 
    histoX=-999.99:0.02:999.99;
    histoY=histcounts(dataset(:,k),edges); 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    dataset(:,k)=dataset(:,k)-f.b1;
    close all; 
    S.sumnpixels(k,2)=max(dataset(:,k));
end

startframe=1:S.lenfile:size(dataset,1);
datasetSm=dataset;
for i=1:size(dataset,2)
    datasetSm(:,i)=medfilt1(dataset(:,i),S.medfiltSm);

    for k=1:length(startframe)
        dataset(startframe(k):startframe(k)+1,i)=NaN; 
        datasetSm(startframe(k):startframe(k)+1,i)=NaN; 
    end
    
    dataset(:,i)=y; 
    datasetSm(:,i)=x; 
end

S.dataset=dataset;
S.datasetSm=datasetSm; 

S.FzeroALL=zeros(size(datasetALL,2),1); S.FzeroALL(:,:)=NaN;

for k=transpose(S.allROIs)

    test=datasetALL(:,k); 
    
    if k==1
        thres=test<S.bar2(1); 
    end
    
    y=test; y(thres)=NaN; 
    test2=y(~thres); 
    edges=-150:0.1:3000;
    histoX=-149.95:0.1:2999.95;
    ha=histogram(test2,edges); hold on; 
    histoY=ha.Values; 
    f = fit(transpose(histoX),transpose(histoY),'gauss1'); 
    S.FzeroALL(k,1)=f.b1; 
   
    x1=[S.FzeroALL(k,1) S.FzeroALL(k,1)]; y1=[0 15];  
    hold on; plot(x1,y1, 'r-', 'LineWidth',2); title(k);
    y=(y-S.FzeroALL(k,1))/(S.FzeroALL(k,1)); 
    datasetALL(:,k)=y; 
    
    close all; 
end

for k=transpose(S.allROIs)
    x=1:size(datasetALL,1); x=transpose(x); 
    y=datasetALL(:,k); 
    
    y(thres)=NaN; 
    filtdata=medfilt1(y,S.medfiltC, 'omitnan'); 
    
    c = polyfit(x,filtdata,S.polynom);
    y_est = polyval(c,x);
    datasetALL(:,k)=datasetALL(:,k)-y_est; 
    
    edges=-100:0.02:100; 
    histoX=-99.99:0.02:99.99;
    histoY=histcounts(datasetALL(:,k),edges); 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    datasetALL(:,k)=datasetALL(:,k)-f.b1;
    close all; 
end

startframe=1:S.lenfile:size(datasetALL,1);
for i=transpose(S.allROIs)
    y=datasetALL(:,i); 
    
    y(thres)=NaN; 
    for k=1:length(startframe)
        y(startframe(k):startframe(k)+1,1)=NaN;
    end
    
    datasetALL(:,i)=y;
end

S.datasetALL=datasetALL;

S.frametimesperlap=cell(1,(length(S.startlevels)-2));
S.dataperlap=cell(1,size(dataset,2));
S.dataperlapCorr=cell(1,size(dataset,2));

for i=2:size(dataset,2)
    test=dataset(:,i); 
    for j=S.lapsused
        y1=S.startlevels(j); y2=(S.startlevels(j+1))-1; 
        C=S.frametiming>y1&S.frametiming<y2; S.frametimesperlap{1,j}=S.frametiming(C); 
        S.frametimesperlap{1,j}=S.frametimesperlap{1,j}-y1; 
        D=find(S.frametiming>y1&S.frametiming<y2); 
        S.dataperlap{1,i}{1,j}=test(D); 
     end
end

clearvars -except S control dataset*

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
