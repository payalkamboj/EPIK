%% Extracting voxel information
%Epik preprocessing

%clear all
%clc

global minPoints epsilon imContour factorM
minPoints = 5;
epsilon = 2;
y = [];
num_patients = 1;
folder = dir('data/ASUAI_0*');
fileSummary = dir('data/*ASUAI_0*'); 
TrainingData = [];
TestData = []; 
trueLabel = [];
trainLabel = [];
foldA = [1:num_patients];
Rejects = [6, 15, 16, 14, 17];
factorX = ones(1,num_patients);

factorX(15) = 2/5;
factorX(22) = 1;
factorX(31) = 1;
PatID = [1:num_patients];
X = 0;
Y = 0;
size_X = 0;
size_Y = 0;
num_Row = 0;
num_Col =0;

NewScalings = [];
NewnumRowInd = [];
NewnumColInd = [];
for foldIMK = 1:num_patients
    foldI = PatID(foldIMK);
    factorM = factorX(foldI);
    fold = foldA(foldI); 
    files = dir(strcat(['data/' folder(fold).name '/MO/report/*_thresh*']));
    firstThresh = files(1);
    imgName = 'IC_1_thresh.png';
    im = imread(strcat(['data/' folder(fold).name '/MO/report/' imgName]));
    %imshow(im)
    templateImage = imread('R.png');
    [X, Y, size_X, size_Y,num_Row,num_Col] = automateSlicing(im,templateImage);
    NewScalings  = [NewScalings;X, Y, size_X, size_Y];
    save('NewScalings.mat', 'NewScalings');
    NewnumRowInd = [NewnumRowInd; num_Row];
    NewnumColInd = [NewnumColInd; num_Col];
end

NewnumRowInd = NewnumRowInd';
NewnumColInd = NewnumColInd';

Scalings = NewScalings;
need = 1;
Offset = 1;
endOff = 6;

numRowInd=NewnumRowInd;
numColInd=NewnumColInd;
problem = [8,9];
factorX = ones(1,76);
factorX(15) = 2/5;
factorX(22) = 1;
factorX(31) = 1;
PatID = [1:76];
%for foldIMK = 54:76
for foldIMK = 1:num_patients
    foldI = PatID(foldIMK);
    factorM = factorX(foldI);
    numCol = numColInd(foldI);
    numRow = numRowInd(foldI);
    labelingVector{foldI} = [];
    fold = foldA(foldI); 
    files = dir(strcat(['data/' folder(fold).name '/MO/report/*_thresh*']));
    n = length(files) ;
    Interest = [];
    listing = dir(strcat(['data/' folder(fold).name '/MO/report/t*.txt']));
    sigLength = 295;
    freq = (1:sigLength).*0.08475;
    tThresh = 70;
    fThresh = 15;
    sThresh = ceil(fThresh/0.08475);
    maxScore = [];
    sigEnergy = [];
    energyPercent = [];
    classSign = [];
    icaData = [];
    filenameString = [];
    label = [];
    label2 = [];
    label3 = [];
    label4 = [];
    label5 = [];
    closenessCluster = 20;
    minClusterIndex = zeros(1,numCol*numRow);
    maxClusterSizeI = 1000*ones(1,numCol*numRow);
    
    if(need)
        for fileNum = 1:n
            tic
            icaData{fileNum} = dlmread(strcat(['data/' folder(fold).name '/MO/report/' listing(fileNum).name]));
            filename = files(fileNum).name;
            fileNameTemp = sscanf(filename,'IC_%d_thresh');
             if(fileNameTemp == 18)
                 disp('here')
             end
            im = imread(strcat(['data/' folder(fold).name '/MO/report/' filename]));
            startX = Scalings(foldI,1); 
            s_value = Scalings(foldI,1);
            startY = Scalings(foldI,2); 
            sizeX = Scalings(foldI,3);
            sizeY = Scalings(foldI,4);
            BrainSize = 1200000^(1/3);
            VoxelSize = 3; 
            imCropped = [];
            k = 1;
            fileInfo = [];
            for ii = 1:numRow
                for jj = 1:numCol
                    imContour{k} = [];
                    imCropped{k} = imcrop(im,[startX,startY,sizeX,sizeY]);
                    imshow(imCropped{k})
                    fileInfo = [fileInfo; [ii,jj,fileNum]];
                    startX = startX + sizeX;
                    k = k + 1;
                end
                startY = startY + sizeY;
                startX = s_value;
            end
            Offset = 1;
            endOff = 6;
            sizeLimit = 5*VoxelSize;
            pixelLimit = VoxelSize * VoxelSize * 15;
            mismatchFactor = 1.3;
            numCluster = [];
            sizeMatch = [];
            bigClusterID = [];
            clusterInfo = [];

            for i = 1:size(imCropped,2)
                if(~isempty(imCropped{i}))
                    if(~Allzeros(imCropped{i},Offset,endOff))

                     
                           [Clusters{i}, newI{i}, newIB{i}] = clusterDetectV1(imCropped{i},Offset,endOff);
                           maxClusterSize = 0;
                           for p = 1:size(Clusters{i},2)
                                  if(maxClusterSize <= size(Clusters{i}{p},2))
                                      maxClusterSize = size(Clusters{i}{p},2);

                                  end

                           end

                           if(maxClusterSizeI(i) > maxClusterSize)
                                maxClusterSizeI(i) = maxClusterSize;
                                minClusterIndex(i) = fileNum;
                           end



                      

                    end
                end
            end



        end

        for fileNum = 1:n
            tic
            icaData{fileNum} = dlmread(strcat(['data/' folder(fold).name '/MO/report/' listing(fileNum).name]));
            filename = files(fileNum).name;
            fileNameTemp = sscanf(filename,'IC_%d_thresh');
             if(fileNameTemp == 18)
                 disp('here')
             end
            im = imread(strcat(['data/' folder(fold).name '/MO/report/' filename]));
            imshow(im);
            startX = Scalings(foldI,1); 
            s_value = Scalings(foldI,1);
            startY = Scalings(foldI,2); 
            sizeX = Scalings(foldI,3);
            sizeY = Scalings(foldI,4);
            BrainSize = 1200000^(1/3);
           
            VoxelSize = 3; 
            imCropped = [];
            k = 1;
            fileInfo = [];
            for ii = 1:numRow
                for jj = 1:numCol
                    
                    imCropped{k} = imcrop(im,[startX,startY,sizeX,sizeY]);
                    %imshow(imCropped{k})
                    fileInfo = [fileInfo; [ii,jj,fileNum]];
                    startX = startX + sizeX;
                    k = k + 1;
                end
                startY = startY + sizeY;
                startX = s_value;
            end
            Offset = 6;
            endOff = 5;
            sizeLimit = 5*VoxelSize;
            pixelLimit = VoxelSize * VoxelSize * 15;
            mismatchFactor = 1.3;
            numCluster = [];
            sizeMatch = [];
            bigClusterID = [];
            clusterInfo = [];

            for i = 1:size(imCropped,2)
                if(~isempty(imCropped{i}))
                    if(~Allzeros(imCropped{i},Offset,endOff))

                       
                          if(minClusterIndex(i) == fileNum)
                             imContour{i} = imCropped{i}; 

                          end



                      

                    end
                end
            end
        end
        save(strcat(['ContourImagesP' folder(fold).name 'Init.mat']),'imContour') 
    else
        
       load(strcat(['ContourImagesP' folder(fold).name 'Init.mat'])) 
       
    end
    for i = 1:numRow*numCol
       if(~isempty(imContour{i}))
            imContour{i} = removeR(imContour{i},Offset,endOff);
       end
   end
%% After contour images are produced    
    %load('imagePixels.mat');
    imagePixels = [];
    for j = 1:size(imContour{20},1)
        for i = 1:size(imContour{20},2)
            imagePixels = [imagePixels; [i j]];
        end
    end
   
    for fileNum = 1:n
        tic
        icaData{fileNum} = dlmread(strcat(['data/' folder(fold).name '/MO/report/' listing(fileNum).name]));
        filename = files(fileNum).name;
        fileNameTemp = sscanf(filename,'IC_%d_thresh');
        [G H] = find(problem == fileNameTemp);
         if(~isempty(G))
             disp('here');
         end
        im = imread(strcat(['data/' folder(fold).name '/MO/report/' filename]));
        %imshow(im);
        startX = Scalings(foldI,1); 
        s_value = Scalings(foldI,1);
        startY = Scalings(foldI,2); 
        sizeX = Scalings(foldI,3);
        sizeY = Scalings(foldI,4);
        BrainSize = 1200000^(1/3);
       
        VoxelSize = 3; %floor(BrainSize/VoxelLength);
        imCropped = [];
        k = 1;
        fileInfo = [];
        for ii = 1:numRow
            for jj = 1:numCol
                
                imCropped{k} = imcrop(im,[startX,startY,sizeX,sizeY]);
                %imshow(imCropped{k})
                [G,H] = find(rgb2gray(imCropped{k}) > 0);
                sizeIM(k) = max(size(G)); 
                fileInfo = [fileInfo; [ii,jj,fileNum]];
                startX = startX + sizeX;
                k = k + 1;
            end
            startY = startY + sizeY;
            startX = s_value;
        end
        
        %% Image Processing approach
        % Call edge detection
        
        sizeLimit = 5*VoxelSize;
        pixelLimit = factorM*VoxelSize * VoxelSize * 15;
        mismatchFactor = 1.3;
        numCluster = [];
        sizeMatch = [];
        bigClusterID = [];
        clusterInfo = [];
        Clusters = [];
        clusOverStat = [];
        peri = [];
        periNum = [];
        periMax = [];
        ClusterWidth = [];
        art = [];
        artMax = [];
        whiteM = [];
        maxID = [];
        percWhite = [];
        ClusterCentroid = [];
        majorClusID = [];
        
        for i = 1:size(imCropped,2)
            clusterInfo2{i} = []; 
            clusterInfo(i,:) = [0 0 0 0];
            if(i == 18)
               disp('here'); 
            end
            numCluster(i) = 0;
            sizeMatch(i) = 0;
            
            if(~isempty(imCropped{i}))
                
                if(~Allzeros(imCropped{i},Offset,endOff))
                    imCropped{i} = removeR(imCropped{i},Offset,endOff);

                    imshow(imCropped{i});
                    [Clusters{i},ClusterWidth{i}, peri{i}, periNum{i}, art{i}, whiteM{i}, maxID(i), clusOverStat{i}, percWhite(i)] = clusterDetect(imCropped{i},Offset,endOff,imContour{i},imagePixels, numRow, numCol, i, sizeIM);
                    imshow(imCropped{i});
                    ClusterCentroid(i,:) = [0 0];
                    periMax(i) = 0;
                    artMax(i) = 0;
                    
                    if(maxID(i)~=0)
                        if(size(Clusters{i}{maxID(i)},1) >= pixelLimit)
                            ClusterCentroid(i,:) = [mean(Clusters{i}{maxID(i)}(:,1)) mean(Clusters{i}{maxID(i)}(:,2))]; 
                            periMax(i) = peri{i}(maxID(i));
                            artMax(i) = art{i}(maxID(i));
                        end
                    end
                    maxSizeClus = 0;
                    maxClusID = 0;
                    for kM = 1:size(Clusters{i},2)
                   
                       if(size(Clusters{i}{kM},1) >= pixelLimit && i >= 8 && i <= 36 )
                           numCluster(i) = numCluster(i) + 1;
                           clusterInfo2{i} = [clusterInfo2{i}; [size(Clusters{i}{kM},1), peri{i}(kM), art{i}(kM), whiteM{i}(kM)]];
                           if(maxSizeClus < size(Clusters{i}{kM},2))
                               maxSizeClus = size(Clusters{i}{kM},2);
                               maxClusID = kM;
                           end
                       end
                    end
                    clusterInfo(i,:) = [size(Clusters{i}{kM},1), peri{i}(kM), art{i}(kM), whiteM{i}(kM)];
                    for kM = 1:size(Clusters{i},2)
                   
                       if(size(Clusters{i}{kM},1) >= pixelLimit && i >= 8 && i <= 36 )
                      
                           majorClusID = [majorClusID i];
                           if(ClusterWidth{i}(kM,1) > sizeLimit || ClusterWidth{i}(kM,2) > sizeLimit )
                               
                                   if(peri{i}(kM) == -1)
                                        sizeMatch(i) = -1;
                                        bigClusterID = [bigClusterID 0];
                                        break;
                                   else
                                       if(art{i}(kM) == -1)
                                           sizeMatch(i) = -1;
                                         
                                          break;
                                       else
                                           if(whiteM{i}(kM) == 1)
                                               sizeMatch(i) = -1;
                                           
                                               break;
                                           else
                                                sizeMatch(i) = 1;
                                                bigClusterID = [bigClusterID 1];
                                           end
                                       end
                                   end
                                   
                               %end
                           end
                           
                       end

                   
                    end 
                    
                    
                    
                end
            end
        end
        
        [G,H] = find(numCluster > 0);
        [G,H] = find(numCluster > 0 & numCluster <= 3);
        [G4,H4] = find(sizeMatch > 0);
        
        [G1,H1] = find(sizeMatch == -1);


        fileNum
        
        tThresh = 90;
        label(fileNum) = 0;
        label2(fileNum) = 0;
        label3(fileNum) = 0;
        label4(fileNum) = 0;
        label5(fileNum) = 0;
        cvFile(fileNum) = std(icaData{fileNum})/mean(icaData{fileNum});
        arteryCVThresh = 2.5;
        [G,H] = max(clusterInfo(:,1));
        [G3,H3] = find(clusterInfo(:,2) == -1);
        if(clusterInfo(H(1),2) == -1)
            label3(fileNum) = 0;
        elseif(clusterInfo(H(1),3) == -1)
            label3(fileNum) = 0;
        elseif(clusterInfo(H(1),4) == 1)
            label3(fileNum) = 0;
        else
            label3(fileNum) = 1;
        end
        if(~isempty(G4))
            label2(fileNum) = 1;
            label4(fileNum) = 1;
            if(size(G,2) >= size(G1,2))
                %if(isempty(G1))
                    label(fileNum) = 1;
                
   
            else
                label(fileNum) = 0;
            end
            if(~isempty(G3) )
               label4(fileNum) = 0; 
            end
        end
        k = 1;
        OverallClusters = [];
        ClusterSize = [];
        OverallPeri = [];
        OverallArt = [];
        OverallWhite = [];
        OverallWhiteOverlap = [];
        OverallPeriNum = [];
        percBrainArray = [];
        percWhiteArray = [];
        for i = 1:size(Clusters,2)
            OverallClusters = [OverallClusters Clusters{i}];
            for lk = 1:size(Clusters{i},2)
               ClusterSize = [ClusterSize size(Clusters{i}{lk},1)]; 
            end
            
            OverallPeri = [OverallPeri peri{i}];
            OverallPeriNum = [OverallPeriNum periNum{i}];
            OverallArt = [OverallArt art{i}];
            OverallWhite = [OverallWhite whiteM{i}];
            OverallWhiteOverlap = [OverallWhiteOverlap clusOverStat{i}];
            percWhiteArray = [percWhiteArray repmat(percWhite(i),1,size(Clusters{i},2))];
            percBrainArray = [percBrainArray repmat(100*(sizeIM(i)/max(sizeIM)),1,size(Clusters{i},2))];
        end

        
        

            

        [B, I] = sort(ClusterSize,'descend');
        %Select top P
        ClusterSizePat{fileNum} = ClusterSize;
        OverallPeriFile{fileNum} = OverallPeri;
        OverallArtFile{fileNum} = OverallArt;
        OverallWhiteFile{fileNum} = OverallWhite;
        percWhiteArrayFile{fileNum} = percWhiteArray;
        percBrainArrayFile{fileNum} = percBrainArray;
        majorityDec = [];
        lk = 1;
        while lk <= size(ClusterSize,2)
            if(ClusterSize(I(lk)) >= pixelLimit)
               if(OverallPeri(I(lk)) ~= -1 && OverallArt(I(lk)) ~= -1 && OverallWhite(I(lk)) ~= 1 && percBrainArray(I(lk)) > 10)
                  majorityDec(lk) = 1; 
               elseif(OverallPeri(I(lk)) ~= -1 && OverallArt(I(lk)) ~= -1 && OverallWhite(I(lk)) == 1 && percWhiteArray(I(lk)) <= 10)
                  majorityDec(lk) = 0;     
               else
                  majorityDec(lk) = -1;

               end
                   
            else
                break;
               %majorityDec(lk) = -1;
            end
            lk = lk + 1;
        end
        
        
        
        [G,H] = find(majorityDec == 1);
        [G2,H2] = find(majorityDec ~= 0);
        
        scoreM = sum(ClusterSize(I(1:size(majorityDec,2))).*majorityDec)./(mean(ClusterSize));
        ClusterNumLimit = 3;
        [GSize, HSize] = find(ClusterSize > pixelLimit);
        if(max(size(G))/max(size(G2)) > 0.5)
        %if(scoreM >= 0 && max(size(GSize)) >= ClusterNumLimit)
           label5(fileNum) = 1; 
        else
           label5(fileNum) = 0;
        end
        
        
            
        
        
        %% Look for temporal coherence of clusters
        needed = 0;
        if(needed)
        AGH = (ClusterCentroid(:,1)+ClusterCentroid(:,2));
        [G,H] = find(AGH > 0);
        Window = 2;
        TGTrue = 0;
        if(min(G) >= 15)
            TG = diff(G) - 1;
            startT = 1;
            endT = startT+Window-1;
            while(endT<=max(size(TG)))
               if(sum(TG(startT:endT)) == 0)
                   TGTrue = 1;
                   break;
               else
                  startT = startT + 1;
                  endT = startT+Window-1;
               end
                
            end
            
            coherenceLimit = 3;
            perOver = [];
            [GSize, HSize] = find(ClusterSize(I) > pixelLimit);
            if(max(ClusterSize(I(HSize))) > 1.6*pixelLimit)
                majorWithin = [];
                if(TGTrue == 1 && size(G,1) > coherenceLimit && (sum(OverallPeri(I(HSize))) == 0) && (sum(OverallArt(I(HSize))) == 0))
                    for kop = 1:size(G,1)-1
                       if(sum((ClusterCentroid(G(kop),:) - ClusterCentroid(G(kop+1),:)).^2 )^0.5 > closenessCluster)
                           majorWithin = [majorWithin kop];
                       end
                    end
                    if((kop == (size(G,1)-1)) && (max(size(majorWithin)) < size(G,1)/2))
                        label5(fileNum) = 1;
       

                    end
                end
            else
                
                label5(fileNum) = 0;
            end
        end
        
        [GW,HW] = find(percWhiteArray(I(1:lk)) > 20);
        
        if(max(size(GW)) <= 2)
            label5(fileNum) = 0;
        end
        
        
        [sortedClusID] = sort(majorClusID,'ascend');
        
        for mn = 1:max(size(sortedClusID))
            maxSizeClus = 0;
            potMax = 0;
            for jk = 1:size(Clusters{sortedClusID(mn)},2)
                tempSize = (size(Clusters{sortedClusID(mn)}{jk},1));
                if(maxSizeClus <= tempSize)
                   maxSizeClus = tempSize; 
                   potMax = jk;
                end
            end
            if(peri{sortedClusID(mn)}(potMax) == -1 || whiteM{sortedClusID(mn)}(potMax) == 1 || art{sortedClusID(mn)}(potMax) == -1 || maxSizeClus < pixelLimit)
                
            else
               if(sortedClusID(mn) > 25)
                  label5(fileNum) = 0; 
               end
               break;
            end
        end
        end
        
        
        
        
        
        
        %%
        
        VR = [];
        score = zeros(1,k-1);
        scoreCluster = zeros(1,k-1);
        scoreAssym = zeros(1,k-1);
        assymScore = [];
        for i = 1:k-1
            if(size(imCropped{i},1) < 160 || size(imCropped{i},1) < 160)
                VR{i} = [];
            else
                VR{i} = Colordetection(imCropped{i},VoxelSize);
            end
        
            dataOrganization = [];
            for k = 1:size(VR{i},1)
                for j = 1:size(VR{i},2)
                    if(VR{i}(k,j) ~= 0)
                        dataOrganization = [dataOrganization; [k j VR{i}(k,j)]]; 
                    end
                end
            end
            if(~isempty(dataOrganization))
                [IDX, COREPTS] = dbscan(dataOrganization, 2, 5);
                 A = unique(IDX);
                 mu = 2; % number of clusters
                 sigma = 2; % std of no of clusters
                %if(max(size(A)) == 2)
                x = max(size(A));
                scoreCluster(i) = exp(-(x-mu).^2/(sigma^2));
                 % scoreCluster(i) = x; %exp(-(x-mu).^2/(sigma^2));
                %else
                    
                    
                %end
            else
                scoreCluster(i) = 0;

            end
            
            assymScore(i) = computeAssymetry(VR{i});
            stdAssymScore = 40;
            
            scoreAssym(i) = exp(-assymScore(i)^2/(stdAssymScore^2));
             
            score(i) = 0.5*scoreCluster(i) + 0.5*scoreAssym(i);
            
            
        end

        
        
        for j = 1:sigLength
            sigEnergy{fileNum}(j) = 100*sum(icaData{fileNum}(1:j))/sum(icaData{fileNum});
        end
        stdThresh = 17.7;
       
            energyPercent(fileNum) = sigEnergy{fileNum}(sThresh);
            classSign(fileNum) = 2-1/(exp((sigEnergy{fileNum}(sThresh)-tThresh)/(stdThresh)));
        
        maxScore(fileNum,:) = [maxk(scoreCluster,3) mink(assymScore,3) sigEnergy{fileNum}(sThresh)];
        [G,H] = find(score > 0.8);
        filenameString(fileNum) = sscanf(filename,'IC_%d_thresh'); 

        toc
    end
    
    
    
    

    [B,I] = sort(filenameString');

    label = label';
    T = table(filenameString',label);
    label2 = label2';
    T2 = table(filenameString',label2);
    label3= label3';
    T3 = table(filenameString',label3);
    label4= label4';
    T4 = table(filenameString',label4);
    label5= label5';
    T5 = table(filenameString',label5);
    
    

     
    %% Computation of precision and recall
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    labelV = label(I);
    labelV2 = label2(I);
    labelV3 = label3(I);
    labelV4 = label4(I);
    labelV5 = label5(I);
   
    
    %% Large number of small clusters
    
    for i = 1:max(size(I))
        ClusterSizePatV2{i} = ClusterSizePat{I(i)};
    end
    
    for i = 1:max(size(I))
        [G,H] = find(ClusterSizePatV2{i} > pixelLimit);
        numClustersOverall(i) = size(H,2);
        if(labelV5(i) == 0)
            numClustersOverall0(i) = size(H,2);
            if(~isempty(H))
                maxSizeClusters0(i) = max(ClusterSizePatV2{i}(H));
                meanSizeClusters0(i) = mean(ClusterSizePatV2{i}(H));
            else
                maxSizeClusters0(i) = 0;
                meanSizeClusters0(i) = 0;
            end
        else
            numClustersOverall1(i) = size(H,2);
            if(~isempty(H))
                maxSizeClusters1(i) = max(ClusterSizePatV2{i}(H));
                meanSizeClusters1(i) = mean(ClusterSizePatV2{i}(H));
            else
                maxSizeClusters1(i) = 0;
                meanSizeClusters1(i) = 0;
            end
        end
    end
    
    [G,H] = find(maxSizeClusters0 > 0);
    maxSizeClusters0NZ = maxSizeClusters0(H);
    
    [G,H] = find(maxSizeClusters1 > 0);
    maxSizeClusters1NZ = maxSizeClusters1(H);
    
    [G,H] = find(meanSizeClusters0 > 0);
    meanSizeClusters0NZ = meanSizeClusters0(H);
    
    [G,H] = find(meanSizeClusters1 > 0);
    meanSizeClusters1NZ = meanSizeClusters1(H);
    for i = 1:max(size(I))
        if(labelV5(i) == 1)
            sizeNorm1(i) = normDist(meanSizeClusters1(i),mean(meanSizeClusters1NZ),std(meanSizeClusters1NZ));
            numNorm(i) = normDist(numClustersOverall(i),mean(numClustersOverall),std(numClustersOverall));
        
            sizeNorm0(i) = normDist(meanSizeClusters1(i),mean(meanSizeClusters0NZ),std(meanSizeClusters0NZ));
            
        end
    end
    
    
    %%
    sigEnergySorted = [];
    for iP = 1:max(size(I))
       sigEnergySorted{iP} = sigEnergy{I(iP)}; 
    end
    
    labelingVector{foldI} = labelV5;  
    save(strcat(['Workspace-' folder(fold).name 'V1.mat']));
end

function y = normDist(x,mu,sigma)
y = exp(-0.5*((x-mu)/sigma)^2);
end
