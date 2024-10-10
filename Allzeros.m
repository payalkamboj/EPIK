function y = Allzeros(image,Offset,endOff)
newI = [];
for i = Offset:size(image,1)-endOff
    for j = Offset:size(image,2)-endOff
        newI(i,j) = sum(image(i,j,:));
        
    end
end
if(sum(sum(newI)) == 0)
    y = 1;
else
    y = 0;
end
