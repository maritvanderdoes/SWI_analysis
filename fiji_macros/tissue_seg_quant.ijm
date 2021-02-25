directory="/Users/Marit/Desktop/untitled folder/";
Table.create("result_table");
start=0;

choices = newArray("Yes", "No");
var choice ="";



filelist = getFileList(directory) 
for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], ".nd")) { 
        open(directory + File.separator + filelist[i]);

        title=getTitle();


		setSlice(1);
		run("Enhance Contrast", "saturated=0.1");
		waitForUser("action required", "inspect the worms age");
		
    	Dialog.create("Select choice");
    	Dialog.addMessage("Does the worm look the good age?");
    	Dialog.addChoice("-->", choices);
    	Dialog.show();
    	choice = Dialog.getChoice();
    	if  (choice == "Yes") {
        	setSlice(2);
			run("Enhance Contrast", "saturated=0");

			waitForUser("action required", "select good z slide");

			run("Duplicate...", "title=formask");
			run("Gaussian Blur...", "sigma=2");
			run("8-bit");
			run("Threshold...");
		
			waitForUser("action required", "set theshold and press ok");
			run("Close");

			run("Convert to Mask");
			run("Fill Holes");
			run("Options...", "iterations=3 count=1 pad do=Erode");
			run("Options...", "iterations=3 count=1 pad do=Dilate");
			run("Analyze Particles...", "size=300-1000 pixel exclude include add");
			selectWindow("formask");
			close();
			run("Z Project...", "projection=[Max Intensity]");
			roiManager("Show All with labels");

			setSlice(3);
			run("Enhance Contrast", "saturated=0.1");
			setSlice(2);
			run("Enhance Contrast", "saturated=0.1");
			waitForUser("action required", "delete incorrect cells, press afterwards ok");
			run("Set Measurements...", "area mean display redirect=None decimal=3");
			roiManager("Measure");

			rows = nResults;
			for (row = 0; row < rows; row++) {
				selectWindow("result_table");
				Table.set("label", row+start, filelist[i]);
				Table.set("channelGFP", row+start, getResult("Mean", row));
				Table.update;
						}	
			selectWindow("Results");
			run("Close");

			setSlice(3);
			run("Enhance Contrast", "saturated=0");
			roiManager("Measure");
	
			for (row = 0; row < rows; row++) {
				selectWindow("result_table");
				Table.set("label", row+start, getResult("Label", row));
				Table.set("channel_mcherry", row+start, getResult("Mean", row));
				Table.update;
						}

			start= start+nResults;
			selectWindow("Results");
			run("Close");
			selectWindow(title);
			run("Close");
			selectWindow("MAX_"+title);
			run("Close");
			roiManager("reset")
    	}
    	else{
        	selectWindow(title);
			run("Close");
			roiManager("reset");
    	}
		
        
    } 
}

Table.save(directory+"results.csv");
print('done')
