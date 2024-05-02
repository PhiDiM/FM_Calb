folder = "basepath"
input = "img\\"
output = "tif\\"
list = getFileList(folder + input);


function cziToTif(input, output, filename)
{
	run("Bio-Formats Importer", "open=" + folder + input + filename + ".czi color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	selectImage(filename + ".czi");
	saveAs("Tiff", folder + output + filename + ".tif");
	close();
}

for (i = 0; i < list.length; i++){
	name = list[i];
	nameinput = substring(name,0,lastIndexOf(name,"."));
	cziToTif(input, output, nameinput);
}