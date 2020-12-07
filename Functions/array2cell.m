function str=array2cell(A)
    str=cell(size(A));
    for j=1:length(A(1,:))
        for i=1:length(A(:,1))
            str{i,j}=num2str(A(i,j));
        end
    end
end