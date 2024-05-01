close("*");
close("ROI Manager");

// get base file name and directory
imFilePath=File.openDialog("Select an image");
open(imFilePath);
fn=getTitle(); 
print(fn);
filename=getInfo("image.filename");
fileNoExtension = File.nameWithoutExtension;
workDir=getInfo("image.directory");
outDir=workDir+"/rotated_cropped_mask";
denoisedDir=workDir+"/n2v_denoise/denoised";
File.makeDirectory(outDir);
greenHeight=getHeight();
close();



// get denoised file to make masks etc
open(denoisedDir+"/"+fileNoExtension+"_green_n2v.tif");
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast", "saturated=0.35");

Dialog.create("Is the image of a head (1=yes)?");
Dialog.addChoice("Head (1=yes)?:", newArray(true,false));
Dialog.show()

isHead=Dialog.getChoice();

//only process if it is a head region
if(isHead){

// get rotation angle
run("Rotate... ");
r=getValue("rotation.angle");
print(r);
rotation=newArray(toString(r));
Array.show("rotation",rotation);
Table.save(outDir+"/rotation_MAX_"+fileNoExtension+".txt");
close("rotation");
close();
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast", "saturated=0.35");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");

// crop worm region
//makeRectangle(0, 335, 976, 274);
makeRectangle(0, 335, 1952, 548);
waitForUser("adjust crop box");
roiManager("add");
run("Crop");


// select pharynx
setTool("polygon");
waitForUser("select head up to pharynx");
roiManager("add");
roiManager("Select", 1);
RoiManager.setGroup(1);
//makeRectangle(72, 241, 863, 26);
makeRectangle(72, 241, 1726, 52);
waitForUser("select background region");
roiManager("add");
roiManager("Select", 2);
RoiManager.setGroup(2);
roiManager("Save", outDir+"/roi_"+fileNoExtension+".zip");

// make head mask
run("Multi-class mask(s) from Roi(s)", "show_mask(s) save_in=[] suffix=[] save_mask_as=tif rm=[RoiManager[size=3, visible=true]]");
saveAs("Tiff",outDir+"/phaMask_"+fileNoExtension+".tif");
close("*");


// show selection on brightfield
open(workDir+"/"+fileNoExtension+"_bf.nd2");
run("Size...", "width="+greenHeight+" height="+greenHeight+" depth=1 constrain average interpolation=Bilinear");
run("Enhance Contrast", "saturated=0.35");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");
roiManager("Select", 0);
run("Crop");
roiManager("Select", 1);
run("Flatten");
saveAs("Tiff",outDir+"/phaBF_"+fileNoExtension+"_bf.tif");
waitForUser;
close("*");

// show selection on green
open(denoisedDir+"/"+fileNoExtension+"_green_n2v.tif");
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast", "saturated=0.35");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");
roiManager("Select", 0);
run("Crop");
roiManager("Select", 1);
run("Flatten");
saveAs("Tiff",outDir+"/phaG_"+fileNoExtension+".tif");
close("*");

// process denoised image for nuclei detection
open(denoisedDir+"/"+fileNoExtension+"_green_n2v.tif");
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast", "saturated=0.35");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");
roiManager("Select", 0);
run("Crop");
saveAs("Tiff",outDir+"/n2v_"+fileNoExtension+".tif");
close("*");


// process denoised image for display
open(denoisedDir+"/"+fileNoExtension+"_green_n2v.tif");
run("Z Project...", "projection=[Max Intensity]");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");
roiManager("Select", 0);
run("Crop");
run("Invert LUT");
saveAs("Tiff",outDir+"/invLUT_n2v_"+fileNoExtension+".tif");
close("*");

// process original raw image
open(workDir+"/"+fileNoExtension+".nd2");
run("Z Project...", "projection=[Max Intensity]");
run("Rotate... ", "angle="+r+" grid=1 interpolation=Bilinear fill enlarge");
roiManager("Select", 0);
run("Crop");
saveAs("Tiff",outDir+"/raw_"+fileNoExtension+".tif");
close("*");
//isHead=false;

}

close("*");
close("ROI Manager");
