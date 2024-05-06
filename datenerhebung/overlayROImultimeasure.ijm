folder = "inputfolder"
input = "img\\"
inputoverlay = "overlay\\"
output = "outputfolder"
list = getFileList(folder + input);

function overlayROImultimeasure(input, output, filename)
{
	run("Set Measurements...", "area mean min centroid redirect=None decimal=3");
	run("Bio-Formats Importer", "open=" + folder + input + filename + ".czi color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	selectImage(filename + ".czi");
	open(folder + inputoverlay + filename + "_cp_masks.png");
	selectImage(filename + "_cp_masks.png");
	run("Label image to ROIs", "rm=[RoiManager[size=147, visible=true]]");
	close();
	selectImage(filename + ".czi");
	roiManager("Show All");
	run("Set Scale...", "distance=1 known=1 unit=pixe");
	roiManager("multi-measure measure_all");
	for (row=0; row<nResults/3; row++) {
		roiManager("Select", row);
		Roi.getContainedPoints(x,y);
    	Table.create("ROIpixels");
		Table.setColumn("ROIpixelsX", x);
		Table.setColumn("ROIpixelsY", y);
		rows = newArray(x.length);
		for (i=0; i<x.length; i++){
			rows[i] = row+1;
		}
		Table.setColumn("iROI", rows);
		Table.save(output + "roi\\" + "ROI_" + filename + row + ".csv", "ROIpixels");
		NumberofRows=Table.size("ROIpixels");
		Table.deleteRows(0, NumberofRows-1, "ROIpixels")
	}
	updateResults();
	saveAs("Results", output + "protein\\" + filename + ".csv");
	run("Close");
	close();
}

for (i = 0; i < list.length; i++){
	name = list[i];
	nameinput = substring(name,0,lastIndexOf(name,"."));
	overlayROImultimeasure(input, output, nameinput);
	
//    	Table.create("ROIpixels");
  //  	Table.save("ROIpixels", output + "roi\\" + "ROI_" + filename + row + ".csv");
}