//Super bare bones function to split channels and subtract background from images as required
//Change channels and naming as needed

function bgsub (input,output,filename){
	open(input + filename);
	run("Split Channels");
	selectWindow("C2-" + filename);
	run("Subtract Background...", "rolling=5");
	saveAs("Tiff", output + filename + "_R5BGs");
	close();
	selectWindow("C1-" + filename);
	close();
	selectWindow("C3-" + filename);
	saveAs("Tiff", output + filename + "_EGFP");
	close();
}

//Add input and output directory paths
input = ""
output= ""

//Get list of files in input directory
list=getFileList(input)

//run on input directory, save output
for(i=0;i<list.length;i++){
	bgsub(input,output,list[i]);
}
