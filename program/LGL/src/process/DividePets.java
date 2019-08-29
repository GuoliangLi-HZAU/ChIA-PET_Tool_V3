package process;

import java.io.File;
import java.io.IOException;

import LGL.chiapet.PetClassification;

public class DividePets {
	
	private Path p;
	private String outPrefix;
	private int sPetNum = 0;
	private int iPetNum = 0;
	private int oPetNum = 0;
	
	public DividePets(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
    public void dividePets() {// divide the PETs into different categories
    	try {
    		File file = new File(outPrefix+".bedpe.selected.pet.txt");
    		if (file.exists()) {
				PetClassification.main(new String[]{outPrefix+".bedpe.selected.pet.txt", outPrefix+".ipet", outPrefix+".spet", outPrefix+".opet", 
						p.SELF_LIGATION_CUFOFF});
    		} else {
    			System.out.println("Error: "+file+" doesn't exist");
    		}
		} catch (IOException e) {
			e.printStackTrace();
		}
    	LinkerFiltering lf = new LinkerFiltering(p);
    	sPetNum = lf.lineNum(outPrefix+".spet");
    	lf.writeFile(outPrefix+".basic_statistics.txt", "Self-ligation PETs\t"+String.valueOf(sPetNum), true);
    	iPetNum = lf.lineNum(outPrefix+".ipet");
    	lf.writeFile(outPrefix+".basic_statistics.txt", "Inter-ligation PETs\t"+String.valueOf(iPetNum), true);
    	oPetNum = lf.lineNum(outPrefix+".opet");
    	lf.writeFile(outPrefix+".basic_statistics.txt", "Other PETs with short distance\t"+String.valueOf(oPetNum), true);
    }
}
