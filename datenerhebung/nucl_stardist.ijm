folder = "folder"
input = "img\\"
output = "outputfolder"
list = getFileList(folder + input);

function overlayROImultimeasure(input, output, filename)
{
	run("Set Measurements...", "area mean min centroid redirect=None decimal=3");
	run("Bio-Formats Importer", "open=" + folder + input + filename + ".czi color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
	selectImage(filename + ".czi - C=1");
	close();
	selectImage(filename + ".czi - C=0");
	close();
	selectImage(filename + ".czi - C=2");
	run("Set Scale...", "distance=1 known=1 unit=pixe");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'" + filename + ".czi - C=2', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.2000000000000002', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Both', 'nTiles':'2', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectImage("Label Image");
	close();
	selectImage(filename + ".czi - C=2");
	roiManager("Show All");
	roiManager("multi-measure measure_all");
	saveAs("Results", output + "nucleus\\" + "nuclei_stardist_" + filename + ".csv");
	run("Close");
	close();
}

for (i = 0; i < list.length; i++){
	name = list[i];
	nameinput = substring(name,0,lastIndexOf(name,"."));
	overlayROImultimeasure(input, output, nameinput);
}