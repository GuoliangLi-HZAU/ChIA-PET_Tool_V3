package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import LGL.chiapet.LinkerFiltering_FastQ_PET;
import LGL.chiapet.LinkerFiltering_FastQ_PET_longread;

public class LinkerFiltering {
	
	private Path p;
	private String outPrefix;
	
	public LinkerFiltering(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
	public void run() {
		try {
			if (p.MODE.equals("0")) {
				new LinkerFiltering_FastQ_PET(p);
				int lineNum = lineNum(outPrefix+".1_1.R1.fastq") + lineNum(outPrefix+".2_2.R1.fastq");
		        writeFile(outPrefix+".basic_statistics.txt", "Same-linker PETs after linker filtering\t"+String.valueOf(lineNum/4), true);
			} else {
				new LinkerFiltering_FastQ_PET_longread(p);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		parameters();
	    if (!p.MODE.equals("0")) {
	    	PeakCalling peakCalling = new PeakCalling(p);
	    	String line = peakCalling.fileLine(outPrefix+".linker_composition_distribution.txt", 2);
	    	if (line != null) {
		    	line = line.trim();
		    	String[] strs = line.split("[ \t]+");
		    	if (strs.length >= 8) {
		    		int total = Integer.valueOf(strs[1])+Integer.valueOf(strs[2])+Integer.valueOf(strs[3])+Integer.valueOf(strs[5])+Integer.valueOf(strs[6])+
		    				Integer.valueOf(strs[7]);
		    		writeFile(outPrefix+".basic_statistics.txt", "Same-linker PETs after linker filtering\t"+String.valueOf(total), true);
		    	}
	    	}
	    }
	}
	
	/**
	 * @author sun
	 * @function write into file
	 * @param f the file to write
	 * @param s the content to write
	 */
	public void writeFile(String f, String s, boolean append) {
		File file = new File(f);
		try {
	        if(!file.exists()){
			    file.createNewFile();
	        }
		    FileWriter fw = new FileWriter(file, append);
		    fw.write(s+"\n");
		    fw.close();
		} catch (IOException e) {
		    e.printStackTrace();  
		}
	}
	
	/**
     * @author sun
     * @function return the number of rows in the file
     * @param file dest file
     * @return
     */
    public int lineNum(String file) {
    	File f = new File(file);
    	int lineNum = 0;
    	if (f.exists()) {
	        try {
	        	BufferedReader reader = new BufferedReader(new FileReader(f));
	        	String line = reader.readLine();
			    while (line != null) {
			    	lineNum++;
			        line = reader.readLine();
			    }
			    reader.close();
		    } catch (IOException e) {
			    e.printStackTrace();
		    }
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
	    return lineNum;
    }
    
    public void parameters() {// runningInformation
		String mapping = "Number of threads\t"+p.NTHREADS+
				"\nGENOME_INDEX directory and prefix\t"+p.GENOME_INDEX;
    	String removingRedundancy = "Mapping quality score cutoff\t"+p.MAPPING_CUTOFF+
    			"\nMerge distance\t"+p.MERGE_DISTANCE;
    	String categories = "Self-ligation cutoff\t"+p.SELF_LIGATION_CUFOFF;
    	String clustering = "Tag extension length\t"+p.EXTENSION_LENGTH+
    			"\nP-value cutoff\t"+p.PVALUE_CUTOFF_INTERACTION;
    	String peakCalling = "Self-ligation cutoff\t"+p.SELF_LIGATION_CUFOFF+
    			"\nMin distance between peak\t"+p.MIN_DISTANCE_BETWEEN_PEAK+
    			"\nPeak mode (1: region, 2: summit)\t"+p.PEAK_MODE+
    			"\nMin coverage for peak\t"+p.MIN_COVERAGE_FOR_PEAK+
    			"\nGenome length\t"+p.GENOME_LENGTH+
    			"\nGenome coverage ratio\t"+p.GENOME_COVERAGE_RATIO+
    			"\nP-value cutoff\t"+p.PVALUE_CUTOFF_PEAK;
    	String report = "Cytoband data\t"+p.CYTOBAND_DATA+ "\nSpecies (1: human, 2: mouse, 3: other)\t"+p.SPECIES;
    	writeFile(outPrefix+".runningInformation.Mapping.txt", mapping, false);
		writeFile(outPrefix+".runningInformation.RemovingRedundancy.txt", removingRedundancy, false);
		writeFile(outPrefix+".runningInformation.Categories.txt", categories, false);
		writeFile(outPrefix+".runningInformation.clustering.txt", clustering, false);
		writeFile(outPrefix+".runningInformation.Peakcalling.txt", peakCalling, false);
		writeFile(outPrefix+".runningInformation.Report.txt", report, false);
	}
    
    public void reset(int step) {
    	if (step > 1 && step < 7) {
    		File f = new File(outPrefix+".basic_statistics.txt");
    		Pvalues pValues = new Pvalues(p);
    		if (f.exists()) {
    			pValues.copyFile(f, new File(outPrefix+".basic_statistics_copy.txt"));
    		} else {
    			f = new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/"+p.OUTPUT_PREFIX+".basic_statistics.txt");
    			if (f.exists()) {
    				pValues.copyFile(f, new File(outPrefix+".basic_statistics_copy.txt"));
    			} else {
        			return;
        		}
    		}
		    try {
		    	BufferedReader reader = new BufferedReader(new FileReader(outPrefix+".basic_statistics_copy.txt"));
		    	String line = reader.readLine();
		    	int n = 0;
		    	new File(outPrefix+".basic_statistics.txt").delete();
			    while (line != null) {
			    	n++;
			    	writeFile(outPrefix+".basic_statistics.txt", line, true);
			    	if (step < 4 && n == 2) {
			    		break;
			    	} else if (step < 5 && n == 5) {
			    		break;
			    	} else if (step < 7 && n == 8) {
			    		break;
			    	}
			        line = reader.readLine();
			    }
			    reader.close();
		    } catch (IOException e) {
			    e.printStackTrace();
		    }
		    new File(outPrefix+".basic_statistics_copy.txt").delete();
    	}
    }
}
