package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;

import LGL.chiapet.BindingSitesFromPETs;
import LGL.chiapet.TagCountForPeaks;
import LGL.chiapet.spetCountForPeaks;

public class PeakCalling {
	
	private Path p;
	private String outPrefix;
	private LinkerFiltering lf;
	
	public PeakCalling(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
		lf = new LinkerFiltering(p);
	}
	
	public void run() {
		try {
			File file = new File(outPrefix+".spet");
			if (file.exists()) {
				BindingSitesFromPETs.main(new String[]{outPrefix+".spet", outPrefix+".peak", p.EXTENSION_LENGTH, p.SELF_LIGATION_CUFOFF, p.MIN_COVERAGE_FOR_PEAK,
						p.PEAK_MODE, p.MIN_DISTANCE_BETWEEN_PEAK});
			} else {
				System.out.println("Error: "+file+" doesn't exist");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
    public void localTag() throws IOException {// generate local tag counts for local density
        File file = new File(outPrefix+".peak");
    	if (file.exists()) {
	        BufferedReader reader = new BufferedReader(new FileReader(file));
	        String line = reader.readLine();
		    while (line != null) {
		    	line = line.trim();
			    String[] strs = line.split("[ \t]+");
			    if (strs.length >= 6) {
		    	    int center = (Integer.valueOf(strs[1])+Integer.valueOf(strs[2]))/2;
		    	    String str1 = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+(center-5000)+"\t"+(center+5000)+"\t"+strs[5];
		    	    String str2 = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+(center-10000)+"\t"+(center+10000)+"\t"+strs[5];
		    	    String[] strs1 = str1.split("\t");
		        	String[] strs2 = str2.split("\t");
	    	        int a1 = 0;
	    	        int b1 = 0;
	    	        int a2 = 0;
	    	        int b2 = 0;
	    	        if (Integer.valueOf(strs1[3]) <= 0) {
	    		        a1 = 1;
	    	        } else {
	    		        a1 = Integer.valueOf(strs1[3]);
	    	        }
	    	        if (Integer.valueOf(strs1[1]) < 0) {
	    		        b1 = 1;
	    	        } else {
	    		        b1 = Integer.valueOf(strs1[1]);
	    	        }
	    	        if (Integer.valueOf(strs2[3]) <= 0) {
	    		        a2 = 1;
	    	        } else {
	    		        a2 = Integer.valueOf(strs2[3]);
	    	        }
	    	        if (Integer.valueOf(strs2[1]) < 0) {
	    		        b2 = 1;
	    	        } else {
	    		        b2 = Integer.valueOf(strs2[1]);
	    	        }
	    	        str1 = strs1[0]+"\t"+b1+"\t"+strs1[2]+"\t"+a1+"\t"+strs1[4]+"\t"+strs1[5];
	    	        str2 = strs2[0]+"\t"+b2+"\t"+strs2[2]+"\t"+a2+"\t"+strs2[4]+"\t"+strs2[5];
	    	        lf.writeFile(outPrefix+".peak.5K_5K.temp", str1, true);
	    	        lf.writeFile(outPrefix+".peak.10K_10K.temp", str2, true);
		        }
			    line = reader.readLine();
		    }
		    reader.close();
    	} else {
    		System.out.println("Error: "+outPrefix+".peak"+" doesn't exist");
    	}
    	new File(outPrefix+".peak.5K_5K").delete();
    	new File(outPrefix+".peak.10K_10K").delete();
    	awkArgind(p.CHROM_SIZE_INFO, outPrefix+".peak.5K_5K.temp", outPrefix+".peak.5K_5K");
    	awkArgind(p.CHROM_SIZE_INFO, outPrefix+".peak.10K_10K.temp", outPrefix+".peak.10K_10K");
		new File(outPrefix+".peak.5K_5K.temp").delete();
		new File(outPrefix+".peak.10K_10K.temp").delete();
		String fileName[] = {"5K_5K", "10K_10K"};
		File file1 = new File(outPrefix+".spet");
		if (file1.exists()) {
			for (int i = 0; i < 2; i++) {
				File file2 = new File(outPrefix+".peak."+fileName[i]);
	    		if (file2.exists()) {
	    			spetCountForPeaks.main(new String[]{outPrefix+".peak."+fileName[i], outPrefix+".spet", outPrefix+".peak.spetCounts."+fileName[i]});
	    		} else {
	    			System.out.println("Error: "+file2+" doesn't exist");
				}
			}
		} else {
    		System.out.println("Error: "+file1+" doesn't exist");
    	}
	}
	
    /**
     * @author sun
     * @function awk
     * @param file1
     * @param file2
     * @param destfile
     */
    public void awkArgind(String file1, String file2, String destFile) {
		File f1 = new File(file1);
		BufferedReader reader = null;
        String line = null;
		ArrayList<String> chromName = new ArrayList<String>();
		ArrayList<String> chromSize = new ArrayList<String>();
		if (f1.exists()) {
	        try {
	        	reader = new BufferedReader(new FileReader(f1));
	        	line = reader.readLine();
			    while (line != null) {
			    	line = line.trim();
				    String[] strs = line.split("[ \t]+");
				    chromName.add(strs[0]);
				    chromSize.add(strs[1]);
				    line = reader.readLine();
			    }
			    reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
    		System.out.println("Error: "+file1+" doesn't exist");
    	}
		File f2 = new File(file2);
		if (f2.exists()) {
	        try {
	        	reader = new BufferedReader(new FileReader(f2));
	        	line = reader.readLine();
			    while (line != null) {
			    	line = line.trim();
				    String[] strs = line.split("[ \t]+");
				    if (strs.length >= 6) {
			        	String b = null;
			        	String c = null;
			        	String d = null;
			        	for (int i = 0; i < chromName.size(); i++) {// find
			        		if (strs[0].equals(chromName.get(i))) {
			        			b = chromSize.get(i);
			        			break;
			        		}
			        	}
			        	if (Integer.valueOf(strs[4]) > Integer.valueOf(b)) {
			        		c = b;
			        	} else {
			        		c = strs[4];
			        	}
                        if (Integer.valueOf(strs[2]) > Integer.valueOf(b)) {
			        		d = b;
			        	} else {
			        		d = strs[2];
			        	}
			        	lf.writeFile(destFile, strs[0]+"\t"+strs[3]+"\t"+c+"\t"+strs[1]+"\t"+d+"\t"+strs[5], true);
			        }
				    line = reader.readLine();
			    }
			    reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
    		System.out.println("Error: "+file2+" doesn't exist");
    	}
	}
    
    public void globalSpet() {// generate the global spet count for global density
    	int sPetNum = lf.lineNum(outPrefix+".spet");
	    lf.writeFile(outPrefix+".nSpets.txt", String.valueOf(sPetNum), false);
	}
    
    public void calculatePvalue() {// calculate p-value
    	String file1 = outPrefix+".peak.spetCounts.5K_5K";
    	String file2 = outPrefix+".peak.spetCounts.10K_10K";
    	pasteSpetCounts(file1, file2);
    	String line = "R --vanilla --slave --args genomeLengthStr="+p.GENOME_LENGTH+"  genomeCoverageRatioStr="+p.GENOME_COVERAGE_RATIO+
    			" extensionLengthStr="+p.EXTENSION_LENGTH+" outPrefix="+outPrefix+" < "+p.PROGRAM_DIRECTORY+"/RScript/pois.r";
    	//lf.writeFile(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh", line, false);
    	lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh", line, false);
    	Shell shell = new Shell();
        //shell.runShell(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh");
    	shell.runShell(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh");
    }
    
    /**
     * @author sun
     * @function pasteSpetCounts
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     */
    public void pasteSpetCounts(String file1, String file2) {
    	File f1 = new File(file1);
    	File f2 = new File(file2);
    	if (f1.exists()) {
    		if (f2.exists()) {
    			try {
	    			BufferedReader reader1 = new BufferedReader(new FileReader(f1));
	    	    	BufferedReader reader2 = new BufferedReader(new FileReader(f2));
	    	    	String line1 = reader1.readLine();
	    	    	String line2 = reader2.readLine();
	    	    	new File(outPrefix+".peak.spetCounts.txt").delete();
	    	    	while (line1 != null && line2 != null) {
	    	    		line1 = line1 + "\t" + line2;
	    	    		line1 = line1.trim();
	    		        String[] strs = line1.split("[ \t]+");
	    		        if (strs.length >= 13) {
	    		        	line1 =  strs[4]+"\t"+strs[5]+"\t"+strs[12];
	    		        	lf.writeFile(outPrefix+".peak.spetCounts.txt", line1, true);
	    		        }
	    	    		line1 = reader1.readLine();
	    	    		line2 = reader2.readLine();
	    	    	}
	    	    	//assume two files have same number of line
	    	    	reader1.close();
	    	    	reader2.close();
    			} catch (IOException e) {
    			    e.printStackTrace();  
    			}
    		} else {
    			System.out.println("Error: "+file1+" doesn't exist");
    		}
    	} else {
			System.out.println("Error: "+file2+" doesn't exist");
		}
    }
    
    public void chromosomalPet() throws IOException {// intra- and inter-chromosomal pet counts for peaks
    	File f1 = new File(outPrefix+".peak");
    	File f2 = new File(outPrefix+".ipet");
    	if (f1.exists()) {
    		if (f2.exists()) {
		    	TagCountForPeaks.main(new String[]{outPrefix+".peak", outPrefix+".ipet", outPrefix+".peak.withTagCounts"});
    		} else {
    			System.out.println("Error: "+f2+" doesn't exist");
    		}
    	} else {
    		System.out.println("Error: "+f1+" doesn't exist");
    	}
    	String file = outPrefix+".peak.withTagCounts";
    	File f = new File(file);
    	BufferedReader reader = null;
        String line = null;
        String str = null;
    	if (f.exists()) {
	        reader = new BufferedReader(new FileReader(f));
		    line = reader.readLine();
		    while (line != null) {
		    	line = line.trim();
		    	String[] strs = line.split("[ \t]+");
		    	if (strs.length >= 7) {
		        	str = strs[5]+"\t"+strs[6]+"\t"+Integer.valueOf(strs[5])+Integer.valueOf(strs[6]);
		        	lf.writeFile(outPrefix+".temp.txt", str, true);
		        }
		    	line = reader.readLine();
	    	}
	    	reader.close();
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    	String file1 = outPrefix+".peak";
    	String file2 = outPrefix+".pvalue.pois.txt";
    	pastePeak1(file1, file2);
    	file = outPrefix+".peak";
    	f = new File(file);
    	if (f.exists()) {
    		reader = new BufferedReader(new FileReader(f));
		    line = reader.readLine();
		    new File(outPrefix+".peak.compact").delete();
		    while (line != null) {
		    	line = line.trim();
		    	String[] strs = line.split("[ \t]+");
		    	if (strs.length >= 6) {
    		    	str = strs[0]+":"+strs[1]+"-"+strs[2]+"\t"+strs[5];
    		    	lf.writeFile(outPrefix+".peak.compact", str, true);
    		    }
		    	line = reader.readLine();
		    }
		    reader.close();
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    	file1 = outPrefix+".peak.compact";
        file2 = outPrefix+".pvalue.pois.txt";
     	file = outPrefix+".temp.txt";
    	pastePeak2(file1, file2, file);
    	file1 = p.PROGRAM_DIRECTORY+"/"+"peakHeader.txt";
        file2 = outPrefix+".peak.long";
        Pvalues pvalues = new Pvalues(p);
        pvalues.mergeFiles(new String[]{file1, file2}, outPrefix+".peak.xls");
        file1 = outPrefix+".peak.5K_5K";
    	file2 = outPrefix+".pvalue.pois.txt";
    	int peakNum = pastePeak3(file1, file2);
    	str = "Peaks from self-ligation\t"+String.valueOf(peakNum);
    	lf.writeFile(outPrefix+".basic_statistics.txt", str, true);
    	int clusterNum = lf.lineNum(outPrefix+".cluster.FDRfiltered.txt");
    	str = "Interacting pairs\t"+String.valueOf(clusterNum);
    	lf.writeFile(outPrefix+".basic_statistics.txt", str, true);
	}
    
    /**
     * @author sun
     * @function pastePeak
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     */
    public void pastePeak1(String file1, String file2) {
    	File f1 = new File(file1);
    	File f2 = new File(file2);
    	if (f1.exists()) {
    		if (f2.exists()) {
    			try {
	    			BufferedReader reader1 = new BufferedReader(new FileReader(f1));
	    	    	BufferedReader reader2 = new BufferedReader(new FileReader(f2));
	    	    	String line1 = reader1.readLine();
	    	    	String line2 = reader2.readLine();
	    	    	new File(outPrefix+".peak.tsv").delete();
	    	    	new File(outPrefix+".peak.bed").delete();
	    	    	new File(outPrefix+".peak.aln").delete();
	    	    	while (line1 != null && line2 != null) {
	    	    		line1 = line1 + "\t" + line2;
	    	    		line1 = line1.trim();
	    		        String[] strs = line1.split("[ \t]+");
	    		        if (strs.length >= 8) {
	    		        	if (Double.valueOf(strs[7]) < Double.valueOf(p.PVALUE_CUTOFF_PEAK)) {
	    		        		lf.writeFile(outPrefix+".peak.tsv", line1, true);
	    		        		double tmp = 0-100*Math.log(Double.valueOf(strs[7]))/Math.log(10);
		        		    	line1 = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t.\t"+formatDouble(tmp)+"\t.\t"+strs[3]+"\t"+strs[4];
		        		    	lf.writeFile(outPrefix+".peak.bed", line1, true);
		        		    	int sum = Integer.valueOf(strs[1])+Integer.valueOf(strs[2]);
		        		    	line1 = strs[0]+"\t"+(sum/2)+"\t"+strs[5];
		        		    	lf.writeFile(outPrefix+".peak.aln", line1, true);
	    		        	}
	    		        }
	    	    		line1 = reader1.readLine();
	    	    		line2 = reader2.readLine();
	    	    	}
	    	    	//assume two files have same number of line
	    	    	reader1.close();
	    	    	reader2.close();
    			} catch (IOException e) {
    			    e.printStackTrace();  
    			}
    		} else {
    			System.out.println("Error: "+file2+" doesn't exist");
    		}
    	} else {
    		System.out.println("Error: "+file1+" doesn't exist");
    	}
    }
    
    /**
     * @author sun
     * @function process decimal digits
     * @param a decimal
     * @return
     */
    public String formatDouble(double a) {
    	String num = String.valueOf(a);
    	int pos = num.indexOf(".");
		pos = 6 - pos;// result pos is the dest digits
		if (pos == 0) {
			num = String.format("%.0f", a);
		} else if (pos == 1) {
			num = String.format("%.1f", a);
		} else if (pos == 2) {
			num = String.format("%.2f", a);	
		} else if (pos == 3) {
			num = String.format("%.3f", a);
		} else if (pos == 4) {
			num = String.format("%.4f", a);
		} else if (pos == 5) {
			num = String.format("%.5f", a);
		}
		double b = Double.valueOf(num);// delete xxx.x0
	    num = String.valueOf(b);
	    pos = num.indexOf(".");// from 999.9 to 1000.0
	    if (num.charAt(pos + 1) == '0' && pos + 1 == num.length() - 1) {// xxx.0 or xxx.0x
	    	num = num.substring(0, pos);
	    }
    	return num;
    }
    
    /**
     * @author sun
     * @function pastePeak
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     * @param file2 the source file to paste
     */
    public void pastePeak2(String file1, String file2, String file3) {
    	File f1 = new File(file1);
    	File f2 = new File(file2);
    	File f3 = new File(file3);
    	if (f1.exists()) {
    		if (f2.exists()) {
    			if (f3.exists()) {
    				try {
	    				BufferedReader reader1 = new BufferedReader(new FileReader(f1));
		    	    	BufferedReader reader2 = new BufferedReader(new FileReader(f2));
		    	    	BufferedReader reader3 = new BufferedReader(new FileReader(f3));
		    	    	String line1 = reader1.readLine();
		    	    	String line2 = reader2.readLine();
		    	    	String line3 = reader2.readLine();
		    	    	new File(outPrefix+".peak.long").delete();
		    	    	while (line1 != null && line2 != null && line3 != null) {
		    	    		line1 = line1 + "\t" + line2 + "\t" + line3;
		    	    		line1 = line1.trim();
		    		        String[] strs = line1.split("[ \t]+");
		    		        if (strs.length >= 4) {
	            		    	if (Double.valueOf(strs[3]) < Double.valueOf(p.PVALUE_CUTOFF_PEAK)) {
	            		    		lf.writeFile(outPrefix+".peak.long", line1, true);
	            		    	}
	            		    }
		    		        line1 = reader1.readLine();
		    	    		line2 = reader2.readLine();
		    	    		line3 = reader3.readLine();
		    	    	}
		    	    	//assume three files have same number of line
		    	    	reader1.close();
		    	    	reader2.close();
		    	    	reader3.close();
    				} catch (IOException e) {
	    			    e.printStackTrace();  
	    			}
    			} else {
    				System.out.println("Error: "+file3+" doesn't exist");
    			}
    		} else {
    			System.out.println("Error: "+file2+" doesn't exist");
    		}
    	} else {
    		System.out.println("Error: "+file1+" doesn't exist");
    	}
    }
    
    /**
     * @author sun
     * @function pastePeak
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     * @return the number of line in *.peak.FDRfiltered.txt
     */
    public int pastePeak3(String file1, String file2) {
    	int num = 0;
    	File f1 = new File(file1);
    	File f2 = new File(file2);
    	if (f1.exists()) {
    		if (f2.exists()) {
    			try {
	    			BufferedReader reader1 = new BufferedReader(new FileReader(f1));
	    	    	BufferedReader reader2 = new BufferedReader(new FileReader(f2));
	    	    	String line1 = reader1.readLine();
	    	    	String line2 = reader2.readLine();
	    	    	new File(outPrefix+".peak.FDRfiltered.txt").delete();
	    	    	while (line1 != null && line2 != null) {
	    	    		line1 = line1 + "\t" + line2;
	    	    		line1 = line1.trim();
	    		        String[] strs = line1.split("[ \t]+");
	    		        if (strs.length >= 8) {
	    		        	if (Double.valueOf(strs[7]) < Double.valueOf(p.PVALUE_CUTOFF_PEAK)) {
	    		        		line1 = strs[0]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t"+strs[6]+"\t"+strs[7];
	    		        		lf.writeFile(outPrefix+".peak.FDRfiltered.txt", line1, true);
	    		        		num++;
	    		        	}
	    		        }
	    	    		line1 = reader1.readLine();
	    	    		line2 = reader2.readLine();
	    	    	}
	    	    	//assume two files have same number of line
	    	    	reader1.close();
	    	    	reader2.close();
    			} catch (IOException e) {
    			    e.printStackTrace();  
    			}
    		} else {
    			System.out.println("Error: "+file2+" doesn't exist");
    		}
    	} else {
    		System.out.println("Error: "+file1+" doesn't exist");
    	}
    	return num;
    }
    
    public void runningTime() {// calculate running time
    	File tfile = new File(outPrefix+".time.txt");
    	if (tfile.exists()) {
	        try {
	        	BufferedReader reader = new BufferedReader(new FileReader(tfile));
	        	String line = reader.readLine();
	        	String content = "";
			    while (line != null) {
				    content = content + line + "\t";
			    	line = reader.readLine();
		    	}
		    	reader.close();
		    	String[] strs = content.split("\t");
		    	if (strs.length >= 7) {
		        	double tmp = (Double.valueOf(strs[1])-Double.valueOf(strs[0]))/60;
		        	String str1 = "Linker filtering\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[2])-Double.valueOf(strs[1]))/60;
		    	    String str2 = "Mapping to genome\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[3])-Double.valueOf(strs[2]))/60;
		    	    String str3 = "Removing redundancy\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[4])-Double.valueOf(strs[3]))/60;
		    	    String str4 = "Categorization of PETs\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[5])-Double.valueOf(strs[4]))/60;
		    	    String str5 = "Interaction calling\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[6])-Double.valueOf(strs[5]))/60;
		    	    String str6 = "Peak calling\t"+String.format("%.2f", tmp);
		    	    tmp = (Double.valueOf(strs[6])-Double.valueOf(strs[0]))/60;
		    	    String str7 = "Total\t"+String.format("%.2f", tmp);
	    	    	lf.writeFile(outPrefix+".running_time.txt", str1, false);
		    	    lf.writeFile(outPrefix+".running_time.txt", str2, true);
		    	    lf.writeFile(outPrefix+".running_time.txt", str3, true);
		    	    lf.writeFile(outPrefix+".running_time.txt", str4, true);
		    	    lf.writeFile(outPrefix+".running_time.txt", str5, true);
		    	    lf.writeFile(outPrefix+".running_time.txt", str6, true);
		    	    lf.writeFile(outPrefix+".running_time.txt", str7, true);
		        }
	        } catch (IOException e) {
			    e.printStackTrace();
	        }
	        tfile.delete();
    	} else {
    		System.out.println("Error: "+tfile+" doesn't exist");
    	}  	
    }
    
    public void summary() throws IOException {// write to *.summary.txt
    	String line = fileLine(outPrefix+".linker_composition_distribution.txt", 2);
    	if (line != null) {
    		line = line.trim();
    	    String[] strs = line.split("[ \t]+");
    	    if (p.MODE.equals("0")) {
    	    	if (strs.length >= 6) {
        	        double a = (Double.valueOf(strs[0])+Double.valueOf(strs[3]))/Double.valueOf(strs[5]);
        	        String str = "(same-linker PETs)/(Total PETs)\t"+String.valueOf(a);
        	        lf.writeFile(outPrefix+".summary.txt", str, false);
    		    }
    	    } else {
    	    	if (strs.length >= 10) {
    	    		double a = (Double.valueOf(strs[1])+Double.valueOf(strs[2])+Double.valueOf(strs[3])+Double.valueOf(strs[5])+Double.valueOf(strs[6])+
    	    				Double.valueOf(strs[7]))/(Double.valueOf(strs[9]));
    	    		String str = "(same-linker PETs)/(Total PETs)\t"+String.valueOf(a);
        	        lf.writeFile(outPrefix+".summary.txt", str, false);
    	    	}
    	    }
    	}
    	File file = new File(outPrefix+".basic_statistics.txt");
    	BufferedReader reader = null;
    	if (file.exists()) {
    	    reader = new BufferedReader(new FileReader(file));
    		line = reader.readLine();
    		String content = "";
    		while (line != null) {
    			String[] strs1 = line.split("\t");
    			content = content + strs1[1] + "\t";
    			line = reader.readLine();
    		}
    		reader.close();
    		String[] strs2 = content.split("[ \t]+");
    		if (strs2.length >= 7) {
		    	String str1 = "(Unique mapped same-linker PETs)/(Total PETs)\t";
		    	String str2 = "(PETs after removing redundancy)/(Total PETs)\t";
		    	String str3 = "(inter-ligation PETs)/(PETs after removing redundancy)\t";
	    		double a = Double.valueOf(strs2[2])/Double.valueOf(strs2[0]);
			    str1 = str1 + String.valueOf(a);
			    double b = Double.valueOf(strs2[4])/Double.valueOf(strs2[0]);
			    str2 = str2 + String.valueOf(b);
			    double c = Double.valueOf(strs2[6])/Double.valueOf(strs2[4]);
			    str3 = str3 + String.valueOf(c);
			    lf.writeFile(outPrefix+".summary.txt", str1, true);
			    lf.writeFile(outPrefix+".summary.txt", str2, true);
			    lf.writeFile(outPrefix+".summary.txt", str3, true);
		    }
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    	file = new File(outPrefix+".cluster.FDRfiltered.txt");
    	if (file.exists()) {
    	    reader = new BufferedReader(new FileReader(file));
    		line = reader.readLine();
    		int n = 0;
    		int linenum1 = 0;
    		int linenum2 = 0;
    		while (line != null) {
    			n++;
    			line = line.trim();
    			String[] strs = line.split("[ \t]+");
    			if (strs.length >= 9) {
    			    if (Integer.valueOf(strs[7]) == 1) {
    			    	linenum1++;
    			    	if (Integer.valueOf(strs[8]) < 1000000) {
    			    		linenum2++;
    			    	}
    			    }
    		    }
    			line = reader.readLine();
    		}
    		reader.close();
	        double a = (double)linenum1/n;
	        String str = "(intra-chromosomal inter-ligation PETs)/(inter-ligation PETs)\t"+String.valueOf(a);
	        lf.writeFile(outPrefix+".summary.txt", str, true);
	        a = (double)linenum2/n;
	        str = "(intra-chromosomal inter-ligation PETs within 1Mb)/" +"(intra-chromosomal inter-ligation PETs)\t"+String.valueOf(a);
	        lf.writeFile(outPrefix+".summary.txt", str, true);
    	} else {
    		String str = "(intra-chromosomal inter-ligation PETs)/(inter-ligation PETs)\t"+"-nan";
    		lf.writeFile(outPrefix+".summary.txt", str, true);
    		str = "(intra-chromosomal inter-ligation PETs within 1Mb)/" +"(intra-chromosomal inter-ligation PETs)\t"+"-nan";
    		lf.writeFile(outPrefix+".summary.txt", str, true);
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    	Deletion deletion = new Deletion(p);
    	deletion.move();
    }
    
    /**
     * @author sun
     * @function get one line of file
     * @param file the path of file
     * @param n the number of line
     * @return line
     */
    public String fileLine(String file, int n) {
    	File f = new File(file);
	 	String line = null;
    	if (f.exists()) {
    		try {
    			BufferedReader reader = new BufferedReader(new FileReader(f));
        	    line = reader.readLine();
	            int linenum = 0;
	            while (line != null) {
	        	    linenum++;
	        	    if (linenum == n) {
	        	    	break;
	        	    }
	        	    line = reader.readLine();
	            }
	            reader.close();
    		} catch (IOException e) {
    		    e.printStackTrace();
    	    }
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    	return line;
    }

    public void statisticsReport() {
    	copy(new File(p.PROGRAM_DIRECTORY+"/ChIA-PET_Report"), new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report"));
    	String line = "R --vanilla --slave --args "+p.PROGRAM_DIRECTORY+" "+p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report "+p.OUTPUT_PREFIX+" "+
    	p.CYTOBAND_DATA+" "+p.SPECIES+" "+p.MODE+" < "+p.PROGRAM_DIRECTORY+"/RScript/ChIA-PET_Report.r";
    	//lf.writeFile(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh", line, false);
    	lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh", line, false);
    	Shell shell = new Shell();
        //shell.runShell(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh");
    	shell.runShell(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh");
        File f = new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".ChIA-PET_Report");
        if (f.exists()) {
        	deleteFile(f);
        }
        copy(new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/ChIA-PET_Report/"), 
        		new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".ChIA-PET_Report"));
        deleteFile(new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/ChIA-PET_Report/"));
    }
    
    /**
	 * @author sun
	 * @function copy file or directory
	 * @param source source file or directory
	 * @param dest dest file or directory
	 */
	public void copy(File source, File dest) {
    	if (source.exists()) {
    		if (source.isDirectory()) {// copy directory
                if (dest.exists()) {
                	 if (dest.isDirectory()) {
                		 File f = new File(dest.getAbsolutePath()+"/"+source.getName());
                		 f.mkdir();
                		 String files[] = source.list();
                         for (String file : files) {
                             File childFile = new File(source, file);
                             copy(childFile, f);
                         }
                	 } else {
                		 System.out.println("Error: "+"can't copy this file");
                	 }
                } else {
                	if (dest.getParentFile().exists()) {
                		dest.mkdir();
                		String files[] = source.list();
                        for (String file : files) {
                            File childFile = new File(source, file);
                            copy(childFile, dest);
                        }
                	} else {// parent directory doesn't exist
            			System.out.println("Error: "+dest+" path doesn't exist");
            		}
                }
            } else {// copy file
            	if (dest.exists()) {
            		if (dest.isDirectory()) {
            			File f = new File(dest.getAbsolutePath()+"/"+source.getName());// use getAbsolutePath to delete /
            			copyFile(source,f);
            		} else {
            			copyFile(source,dest);
            		}
            	} else {
            		if (dest.getParentFile().exists()) {// dest file or directory doesn't exist, default copy file
    					copyFile(source,dest);
            		} else {// parent directory doesn't exist
            			System.out.println("Error: "+dest+" path doesn't exist");
            		}
            	}
            }  
    	} else {
    		System.out.println("Error: "+source+" doesn't exist");
    	}
    }
    
	/**
     * @author sun
     * @function copy file
     * @param source source file
     * @param dest dest file
     */
	public void copyFile(File source, File dest) {
    	if (source.exists()) {
    		try {
    			FileInputStream fis = new FileInputStream(source);
    			FileChannel inputChannel = fis.getChannel();
    			FileOutputStream fos = new FileOutputStream(dest);
    			FileChannel outputChannel = fos.getChannel();
    			outputChannel.transferFrom(inputChannel, 0, inputChannel.size());
    			inputChannel.close();
    			outputChannel.close();
    			fis.close();
    			fos.close();
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
    	} else {
    		System.out.println("Error: "+source+" doesn't exist");
    	}
    }
	
	/**
	 * @author sun
	 * @function delete the previous files
	 * @param f the path of file
	 */
	public void deleteFile(File f) {
		if (f.exists()) {
		    if (f.isFile()) {
		      f.delete(); 
		    } else if (f.isDirectory()) {
		        File[] files = f.listFiles();// put all files of directory into files[]
		        for (int i = 0;i < files.length; i++) {
		            deleteFile(files[i]);
		        }  
		        f.delete();// delete folder
		    }  
		}
	}
	
    public void genomicBrowser() throws IOException {// UCSC
    	String str = "browser position chr1:9997500-10922500\n"+
    "browser hide all\n"+
    "track name=ChIAPET description=\"ChIA-PET "+p.OUTPUT_PREFIX+" \"";
    	lf.writeFile(outPrefix+".cluster.toUCSC.gff", str, false); 	
    	File file = new File(outPrefix+".cluster.FDRfiltered.txt");
    	BufferedReader reader = null;
    	String line = null;
    	if (file.exists()) {
    	    reader = new BufferedReader(new FileReader(file));
		    line = reader.readLine();
		    int linenum = 0;
		    while (line != null) {
			    linenum++;
			    line = line.trim();
		        String[] strs = line.split("[ \t]+");
		        if (strs.length >= 6) {
		    	    if (strs[0].equals(strs[3])) {
		    		    str = strs[0]+"\tChIAPET\t"+p.OUTPUT_PREFIX+"\t"+strs[1]+"\t"+strs[2]+"\t.\t.\t.\ttouch"+String.valueOf(linenum)+"\n"+strs[3]+
		    		    		"\tChIAPET\t"+p.OUTPUT_PREFIX+"\t"+strs[4]+"\t"+strs[5]+"\t.\t.\t.\ttouch"+String.valueOf(linenum);
		    		    lf.writeFile(outPrefix+".cluster.toUCSC.gff", str, true);
		    	    }
		        }
		        line = reader.readLine();
		    }
		    reader.close();
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
		str = "browser position chr17:73626165-73912870\n"+
    	"browser hide all\n"+
    	"browser pack refGene encodeRegions\n"+
    	"browser full altGraph\n"+
		"#300 base wide bar graph, autoScale is on by default == graphing\n"+
		"#limits will dynamically change to always show full range of data\n"+
		"#in viewing window, priority = 20 positions this as the second graph\n"+
		"#Note, zero-relative, half-open coordinate system in use for bedGraph format\n"+
		"track type=bedGraph name=\""+p.OUTPUT_PREFIX+" ChIAPET peak\" description=\""+p.OUTPUT_PREFIX+" ChIA-PET peaks\""+
		" visibility=full color=32,178,170 altColor=0,255,127 priority=20";
		lf.writeFile(outPrefix+".peak.toUCSC.bedGraph", str, false);
		file = new File(outPrefix+".peak.FDRfiltered.txt");
		if (file.exists()) {
		    reader = new BufferedReader(new FileReader(file));
		    line = reader.readLine();
		    while (line != null) {
		    	line = line.trim();
		        String[] strs = line.split("[ \t]+");
		        if (strs.length >= 4) {
		    	    str = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3];
		    	    lf.writeFile(outPrefix+".peak.toUCSC.bedGraph", str, true);
		        }
		        line = reader.readLine();
		    }
		   reader.close();
		} else {
			System.out.println("Error: "+file+" doesn't exist");
		}
    }
}
