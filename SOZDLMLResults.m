function [Precision, Recall, F1Score, Accuracy, Specificity, TP, TN, FP, FN] = SOZDLMLResults(trueLabels,PredictedLabels)
TP = 0;
FN = 0;
FP = 0;
TN = 0;

for label = 1:size(trueLabels,1)
    if(trueLabels(label)~= 2)
       if(trueLabels(label) == 3 && PredictedLabels(label) == 3 )
           TP = TP + 1;
       elseif(trueLabels(label) == 3 && PredictedLabels(label) == 0 )
           FN = FN + 1;
       elseif(trueLabels(label) == 3 && PredictedLabels(label) == 1 )
           FN = FN + 1;
       elseif(trueLabels(label) == 0 && PredictedLabels(label) == 3 )
           FP = FP + 1;
       elseif(trueLabels(label) == 1 && PredictedLabels(label) == 3 )
           FP = FP + 1;
       elseif(trueLabels(label) == 0 && PredictedLabels(label) == 0)
           TN = TN + 1;
       elseif(trueLabels(label) == 0 && PredictedLabels(label) == 1)
           TN = TN + 1;
       elseif(trueLabels(label) == 1 && PredictedLabels(label) == 0)
           TN = TN + 1;
       elseif(trueLabels(label) == 1 && PredictedLabels(label) == 1)
           TN = TN + 1;
           
       end
       
    end
           
end
    Precision = TP/(TP + FP);
    Recall = TP/(TP + FN);
    F1Score = 2*(1/(1/Precision + 1/Recall));
    Accuracy = (TP+TN)/(TP+TN+FP+FN);
    Specificity = TN/(TN+FP);