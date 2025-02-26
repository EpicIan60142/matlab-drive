function latexFormatting = arrToBMat(A)
    strMat = string(A);
    rows = size(A,1);
    cols = size(A,2);
    latexFormatting = "\begin{bmatrix}"+newline;
    for i = 1:rows
        latexFormatting = latexFormatting+"    ";
        for j = 1:cols-1
            latexFormatting = latexFormatting + A(i,j)+" & ";
        end
        latexFormatting = latexFormatting+A(i,cols)+"\\"+newline;
    end
    latexFormatting = latexFormatting + "\end{bmatrix}";
end