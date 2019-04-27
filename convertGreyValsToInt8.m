function convertedGrayVals= convertGreyValsToInt8(greyValuesDicom)
%CONVERTGREYVALS scale dicom image values to 0-255
greyValuesDicom = double(greyValuesDicom);
greyValuesDicom = greyValuesDicom / max(greyValuesDicom(:));
convertedGrayVals= uint8(255* greyValuesDicom);
end

