%EPIK 
clc
clear all
close all
numSubjects = 1
PatID = [1:numSubjects];
folder = dir('data/ASUAI_0*');
foldA = [1:numSubjects];
AgePat = [9.91666666700000	16.4166666700000	14.8333333300000	11.9166666700000	13.0833333300000	9	3.33333333300000	3.50000000000000	12.6666666700000	12.5833333300000	2.75000000000000	12.5000000000000	14.5833333300000	4.50000000000000	0.583333333000000	5.25000000000000	5.58333333300000	1.41666666700000	1.33333333300000	14.5833333300000	14.8333333300000	12.8333333300000	18	2.58333333300000	3.75000000000000	3.75000000000000	14	14.6666666700000	7.91666666700000	8.33330000000000	2	4.33333333300000	6.66666666700000	9.25000000000000	16.2500000000000	10.0833333333333	8.16666666666667	15.5000000000000	10.4166666666667	3.16666666666667	17.7500000000000	11.6666666666667	4.75000000000000	18	10.4166666666667	13.6666666666667	2.58333333333333	2.58333333333333	0.250000000000000	2.83333333333333	2.91666666666667	2.08333333333333];
OverallDataNumClusters = [];
OverallDataAssymetryScore = [];
OverallDataGiniI = [];
OverallDataGiniS = [];
OverallDataLabels = [];
OverallUSALabels = [];
OverallAge = [];
OverallGender = [];

TrainDataNumClusters = [];
TrainDataAssymetryScore = [];
TrainDataGiniI = [];
TrainDataGiniS = [];
TrainDataLabels = [];
TrainUSALabels = [];
TrainAge = [];
TrainGender = [];

perc = 60;

TestDataNumClusters = [];
TestDataAssymetryScore = [];
TestDataGiniI = [];
TestDataGiniS = [];
TestDataLabels = [];
TestUSALabels = [];
TestAge = [];
TestGender = [];
MF = [1	1	1	2	1	2	2	1	1	2	1	1	2	1	2	2	1	2	2	2	1	1	2	2	1	2	1	1	1	2	1	2	2	2	2	2	1	2	2	2	1	1	2	2	2	1	2	1	2	2	1	2];

Range1 = [0 5];
Range2 = [5 13];
Range3 = [13 20];
A1Size = [];
A2Size = [];
A3Size = [];
MSize = [];
FSize = [];
load('svmModel.mat')
for foldIMKV = 1:numSubjects
    PatID = [1:numSubjects];
    foldA = [1:numSubjects];
    foldI = PatID(foldIMKV);
    fold = foldA(foldI); 
    if(fold == 31)
        disp('problematic');
    end
    folder = dir('data/ASUAI_0*');
    files = dir(strcat(['data/' folder(fold).name '/MO/report/*_thresh*']));
    n = length(files) ;
   

    
    load(strcat(['Workspace-' folder(fold).name 'V1.mat']));
    
    PatID = [1:numSubjects];
    foldA = [1:numSubjects];
    foldI = PatID(foldIMKV);
    fold = foldA(foldI); 
    folder = dir('data/ASUAI_0*');
    files = dir(strcat(['data/' folder(fold).name '/MO/report/*_thresh*']));
    n = length(files) ;
    listing = dir(strcat(['data/' folder(fold).name '/MO/report/t*.txt']));
    tic
    numberOfClusters = [];
    numberOfClustersP = [];
    assymetryScore = [];
    GiniI = [];
    GiniS = [];
    for fileNum = 1:n
        filename = files(fileNum).name;
        icaData{fileNum} = dlmread(strcat(['data/' folder(fold).name '/MO/report/' listing(fileNum).name]));

        %% Obtaining the images

        im = imread(strcat(['data/' folder(fold).name '/MO/report/' filename]));
        %imshow(im);
        NewScalings=load('NewScalings.mat'); 
        Scalings = NewScalings.NewScalings;
        
        startX = Scalings(foldI,1); 
        startY = Scalings(foldI,2); 
        sizeX = Scalings(foldI,3);
        sizeY = Scalings(foldI,4);
          
        BrainSize = 1200000^(1/3);
        %VoxelLength = 1;
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
            startX = 1;
        end

        %% Computing assymetry

        VR = [];
        score = zeros(1,k-1);
        scoreCluster = zeros(1,k-1);
        scoreAssym = zeros(1,k-1);
         assymScore = [];


        %% BOLD SIGNAL SPARSITY IN ACTIVELET DOMAIN
        lengthSWT = 256;
        GEx{fileNum} = [];
        GiniI(fileNum) = 0;
        for g = 1:lengthSWT:size(icaData{fileNum},1)-lengthSWT
            signalB = icaData{fileNum}(g:g+lengthSWT-1,1);
            SWC = swt(signalB,3,'bior2.2');
            GExLine = lasso(SWC',signalB,'Lambda',0.5);
            BN = GiniIndex(GExLine);
            if(GiniI(fileNum) < BN)
                GiniI(fileNum) = BN;
                GEx{fileNum} = GExLine;
            end
            
            

        end
        

        %% Sparsity in sine domain

        dictionary2 = {'sin'};
        [mpdict2,nbvect2] = wmpdictionary(length(icaData{fileNum}),'lstcpt',dictionary2);
        y2 = wmpalg('BMP',icaData{fileNum},mpdict2,'itermax',35);
        GiniS(fileNum) = GiniIndex(y2);

         
            [G,H] = find(ClusterSizePat{fileNum} > pixelLimit);
            numberOfClusters(fileNum) = max(size(G));

            
            [nL,scores] = predict(cl,[numberOfClusters(fileNum), GiniS(fileNum),GiniI(fileNum)]);
            labelSVM(fileNum) = nL;
%uncomment below to get soz as 1s (currenty l we have networks as 1) and
%comment next if condition if you uncomment this
             if(label5(fileNum) == 1 && nL == 0 && scores(1,1) > 0.4)
                 label5(fileNum) = 0; %making RSN as label 0 here so that remaining 1s become SOZ automatically?
                 
 
             end



         filenameString(fileNum) = sscanf(filename,'IC_%d_thresh');
    end
    toc


        
        %%
        [B,I] = sort(filenameString');
        labelV6 = label5(I);
        GiniISorted = GiniI(I);
        sum(labelV6-labelV5)
        OverallDataGiniI = [OverallDataGiniI GiniISorted];
        TrainDataGiniI = [TrainDataGiniI GiniISorted(1:floor((perc/100)*max(size(GiniISorted))))];
        TestDataSVMGiniI{foldI} = GiniISorted;
        TestDataSVMGiniS{foldI} = GiniISorted;
        TestDataSVMNum{foldI} = numberOfClusters;
        TestDataGiniI{foldI} = GiniISorted(floor((perc/100)*max(size(GiniISorted)))+1:end);
        GiniSSorted = GiniS(I);
        OverallDataGiniS = [OverallDataGiniS GiniSSorted];
        TrainDataGiniS = [TrainDataGiniS GiniSSorted(1:floor((perc/100)*max(size(GiniSSorted))))];
        TestDataGiniS{foldI} = GiniSSorted(floor((perc/100)*max(size(GiniSSorted)))+1:end);
        numberOfClusters = numberOfClusters(I);
        OverallDataNumClusters = [OverallDataNumClusters numberOfClusters];
        TrainDataNumClusters = [TrainDataNumClusters numberOfClusters(1:floor((perc/100)*max(size(GiniSSorted))))];
        TestDataNumClusters{foldI} = numberOfClusters(floor((perc/100)*max(size(GiniSSorted)))+1:end);
        
       
        PatID = [1:76];
        foldA = [1:76];
        foldI = PatID(foldIMKV);
        fold = foldA(foldI); 
        folder = dir('data/ASUAI_0*');
        files = dir(strcat(['data/' folder(fold).name '/MO/report/*_thresh*']));
        n = length(files) ;
        listing = dir(strcat(['data/' folder(fold).name '/MO/report/t*.txt']));
        OverallUSALabels = [OverallUSALabels labelV6'];
        
        indices_of_SOZ = find(labelV6 == 1)
       
        filenameExcel = strcat([folder(fold).name '_PotentialSOZICs.csv']); 
       
        % Write data to a CSV file
        data = [indices_of_SOZ.']; 
        csvwrite(filenameExcel, data);


 

       
    
end




