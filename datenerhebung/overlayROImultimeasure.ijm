folder = "basepath"
input = "img\\"
inputoverlay = "overlay\\"
output = "result\\"
list = getFileList(folder + input);

function overlayROImultimeasure(input, output, filename)
{
	run("Bio-Formats Importer", "open=" + folder + input + filename + ".czi color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	selectImage(filename + ".czi");
	open(folder + inputoverlay + filename + "_cp_masks.png");
	selectImage(filename + "_cp_masks.png");
	run("Label image to ROIs", "rm=[RoiManager[size=147, visible=true]]");
	close();
	selectImage(filename + ".czi");
	roiManager("Show All");
	roiManager("multi-measure measure_all");
	saveAs("Results", folder + output + filename + ".csv");
	run("Close");
	close();
}

for (i = 0; i < list.length; i++){
	name = list[i];
	nameinput = substring(name,0,lastIndexOf(name,"."));
	overlayROImultimeasure(input, output, nameinput);
}