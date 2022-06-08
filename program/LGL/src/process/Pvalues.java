package process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import LGL.shortReads.TagCountInGivenRegions;

public class Pvalues {
	
    private Path p;
    private String outPrefix;
    private LinkerFiltering lf;
    
	public Pvalues(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
		lf = new LinkerFiltering(p);
	}
	
	public void calculation() {// calculation of p-values
		p.checkPath(outPrefix+".ipet");
		p.checkPath(outPrefix+".cluster.filtered");
		
		cutIPet(outPrefix+".ipet");
		cutCluster(outPrefix+".cluster.filtered");
		String fileName[] = {"anchor1", "anchor2"};
		for (int i = 0; i < 2; i++) {
		    try {
		    	File file1 = new File(outPrefix+".aln");
		    	File file2 = new File(outPrefix+".cluster.filtered."+fileName[i]);
		    	if (file1.exists()) {
		    		if (file2.exists()) {
		    			if (p.MODE.equals("0")) {
		    				System.out.println("Short reads mode process clusters filter!!!!\n");
		    				TagCountInGivenRegions.main(new String[]{outPrefix+".aln", outPrefix+".cluster.filtered."+fileName[i], outPrefix+
		    						".cluster.filtered."+fileName[i]+".tagCount.txt", p.EXTENSION_LENGTH, p.EXTENSION_MODE}); //"500"
		    			} else {
		    				System.out.println("Long reads mode process clusters filter!!!!\n");
		    				//TagCountInGivenRegions.main(new String[]{outPrefix+".aln", outPrefix+".cluster.filtered."+fileName[i], outPrefix+
		    				//		".cluster.filtered."+fileName[i]+".tagCount.txt", "501", "1"});
		    				TagCountInGivenRegions.main(new String[]{outPrefix+".aln", outPrefix+".cluster.filtered."+fileName[i], outPrefix+
				    				".cluster.filtered."+fileName[i]+".tagCount.txt", p.EXTENSION_LENGTH, p.EXTENSION_MODE}); //"500"
		    			}
		    			
		    		} else {
		    			System.out.println("Error: "+file2+" doesn't exist");
					}
		    	} else {
		    		System.out.println("Error: "+file1+" doesn't exist");
		    	}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * @author sun
	 * @function cut ipet
	 * @param file the path of ipet file
	 */
	public void cutIPet(String source) {
		File f = new File(source);
		if (f.exists()) {
		    try {
		    	BufferedReader reader = new BufferedReader(new FileReader(f));
		    	String line = reader.readLine();
		    	new File(outPrefix+".aln1").delete();
		    	new File(outPrefix+".aln2").delete();
		    	BufferedWriter outaln1 = new BufferedWriter(new FileWriter(outPrefix+".aln1", false));
		    	BufferedWriter outaln2 = new BufferedWriter(new FileWriter(outPrefix+".aln2", false));
			    while (line != null) {
			    	line = line.trim();
			    	String[] strs = line.split("[ \t]+");
			    	if (strs.length >= 6) {
			    		String str = strs[0]+"\t"+strs[1]+"\t"+strs[2];
			    		//lf.writeFile(outPrefix+".aln1", str, true);
			    		outaln1.write(str);
			    		outaln1.newLine();
			    		str = strs[3]+"\t"+strs[4]+"\t"+strs[5];
			    		//lf.writeFile(outPrefix+".aln2", str, true);
			    		outaln2.write(str);
			    		outaln2.newLine();
			    	}
			    	line = reader.readLine();
			    }
			    outaln1.close();
			    outaln2.close();
			    mergeFiles(new String[]{outPrefix+".aln1", outPrefix+".aln2"}, outPrefix+".aln");
			    new File(outPrefix+".aln1").delete();
			    new File(outPrefix+".aln2").delete();
			    reader.close();
		    } catch (IOException e) {
			    e.printStackTrace();  
			}
		} else {
			System.out.println("Error: "+source+" doesn't exist");
		}
	}
	
	/**
	 * @author sun
	 * @function merge files, one by one
	 * @param files source files
	 * @param file dest file
	 */
	public void mergeFiles(String[] files, String file) {
        try {
        	FileOutputStream fos = new FileOutputStream(file);
        	FileChannel destfc = fos.getChannel();
            for (String f : files) {
            	File srcFile = new File(f);
            	if (srcFile.exists()) {
	            	FileInputStream fis = new FileInputStream(f);
	            	FileChannel srcfc = fis.getChannel();
	                ByteBuffer buf = ByteBuffer.allocate(1024);  
	                while (srcfc.read(buf) != -1) {  
	                    buf.flip();// reset index of buf
	                    destfc.write(buf);  
	                    buf.clear();  
	                }  
	                srcfc.close();
	                fis.close();
            	} else {
            		System.out.println("Error: "+f+" doesn't exist");
            	}
            }
            destfc.close();
            fos.close();
        } catch (IOException e) {
		    e.printStackTrace();  
		}
    }
	
	/**
	 * @author sun
	 * @function cut cluster
	 * @param file the path of cluster file
	 */
	public void cutCluster(String file) {
		File f = new File(file);
		if (f.exists()) {
		    try {
		    	BufferedReader reader = new BufferedReader(new FileReader(f));
		    	String line = reader.readLine();
		    	new File(outPrefix+".cluster.filtered.anchor1").delete();
		    	new File(outPrefix+".cluster.filtered.anchor2").delete();
		    	BufferedWriter outanchor1 = new BufferedWriter(new FileWriter(outPrefix+".cluster.filtered.anchor1", false));
		    	BufferedWriter outanchor2 = new BufferedWriter(new FileWriter(outPrefix+".cluster.filtered.anchor2", false));
			    while (line != null) {
			    	line = line.trim();
			    	String[] strs = line.split("[ \t]+");
			    	if (strs.length >= 13) {
			    		String str = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[12];
			    		//lf.writeFile(outPrefix+".cluster.filtered.anchor1", str, true);
			    		outanchor1.write(str);
			    		outanchor1.newLine();
			    		str = strs[6]+"\t"+strs[7]+"\t"+strs[8]+"\t"+strs[12];;
			    		//lf.writeFile(outPrefix+".cluster.filtered.anchor2", str, true);
			    		outanchor2.write(str);
			    		outanchor2.newLine();
			    	}else {
			    		System.out.println("Unexpected cluster: " + line);
			    	}
			    	line = reader.readLine();
			    }
			    reader.close();
			    outanchor1.close();
			    outanchor2.close();
		    } catch (IOException e) {
			    e.printStackTrace();  
			}
		} else {
			System.out.println("Error: "+file+" doesn't exist");
		}
	}
	
	public void globalTag() {// generate the global tag count for global density
		int iPetNum = lf.lineNum(outPrefix+".ipet");
		lf.writeFile(outPrefix+".nTags.txt", String.valueOf(iPetNum * 2), false);
	}
	
    public void calculate() {// calculate p-value
    	//p.checkPath(outPrefix+".cluster"); ??? why qw
    	
    	String file1 = outPrefix+".cluster.filtered.anchor1.tagCount.txt";
    	String file2 = outPrefix+".cluster.filtered.anchor2.tagCount.txt";
    	pasteTagCount(file1, file2);
    	String line = "R --vanilla --slave --args "+outPrefix+" < "+p.PROGRAM_DIRECTORY+"/RScript/hypergeometric.r";
    	//lf.writeFile(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".pvalues.sh", line, false);
    	lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".pvalues.sh", line, false);
    	Shell shell = new Shell();
        //shell.runShell(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".pvalues.sh");
    	shell.runShell(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".pvalues.sh");
    	file1 = outPrefix+".cluster.filtered";
    	file2 = outPrefix+".petCount.tagCount.txt";
    	String file3 = outPrefix+".pvalue.hypergeo.txt";
    	p.checkPath(file1);
    	p.checkPath(file2);
    	p.checkPath(file3);
    	
    	pasteCluster(file1, file2, file3);
		File file = new File(outPrefix+".cluster.FDRfiltered.txt");
		if (file.exists()) {
			try {
		        BufferedReader reader = new BufferedReader(new FileReader(file));
		        line = reader.readLine();
		        Map<String, Integer> map = new HashMap<String, Integer>();
		        String destFile1 = outPrefix+".calculate.temp1.txt";// delete at last
		        String destFile2 = outPrefix+".calculate.temp2.txt";
			    File df1 = new File(destFile1);
			    File df2 = new File(destFile2);
			    while (line != null) {
			    	line = line.trim();
			        String[] strs = line.split("[ \t]+");
			        if (strs.length >= 8) {
			        	if (map.containsKey(strs[6]+"\t"+strs[7])) {
		    	        	map.put(strs[6]+"\t"+strs[7], map.get(strs[6]+"\t"+strs[7]) + 1);
		    	        } else {
		    	        	map.put(strs[6]+"\t"+strs[7], 1);
		    	        }
			        }
			        line = reader.readLine();
			    }
			    reader.close();
			    Iterator<Map.Entry<String, Integer>> iterator = map.entrySet().iterator();
			    Map.Entry<String, Integer> entry;
			    while (iterator.hasNext()) {
			    	entry = iterator.next();
			    	String k = entry.getKey();
			    	Integer v = entry.getValue();
			    	lf.writeFile(destFile1, v + " " + k, true);
			    }
			    reader = new BufferedReader(new FileReader(df1));
			    line = reader.readLine();
			    int a = 0, b = 0, f = 0, g = 0;
			    int[] c = new int [10];
			    int[] d = new int [10];
			    boolean[] change = new boolean [10];
			    for (int i = 0; i < c.length; i++) {
			    	c[i] = 0;
			    	d[i] = 0;
			    	change[i] = false;
			    }
			    while (line != null) {// process the file from uniq -c
			    	line = line.trim();
			        String[] strs = line.split("[ \t]+");
			        if (strs.length >= 3) {
			        	if (Integer.valueOf(strs[1]) >= 10) {
			        		a = a + Integer.valueOf(strs[0]);
			        		b = b + Integer.valueOf(strs[0]) * Integer.valueOf(strs[2]);
			        		f = f + Integer.valueOf(strs[0]);
			        		g = g + Integer.valueOf(strs[0]) * Integer.valueOf(strs[2]);
			        	} else {
			        		change[Integer.valueOf(strs[1])] = true;
			        		c[Integer.valueOf(strs[1])] = c[Integer.valueOf(strs[1])] + Integer.valueOf(strs[0]);
			        		d[Integer.valueOf(strs[1])] = d[Integer.valueOf(strs[1])] + Integer.valueOf(strs[0]) * Integer.valueOf(strs[2]);		        		
			        		f = f + Integer.valueOf(strs[0]);
			        		g = g + Integer.valueOf(strs[0]) * Integer.valueOf(strs[2]);
			        	}
			        }
			        line = reader.readLine();
			    }
			    reader.close();
			    for (int i = 0; i < change.length; i++) {
			    	if (change[i] == true) {
			    		lf.writeFile(destFile2, String.valueOf(c[i])+"\t"+i+"\t"+String.valueOf(d[i]), true);
			    	}
			    }
			    lf.writeFile(destFile2, a+"\t10\t"+b, true);
			    lf.writeFile(destFile2, f+"\t11\t"+g, true);
			    df1.delete();
			    Purifying purifying = new Purifying(p);
			    purifying.sortN(destFile2, 2);
			    String str = "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters";
			    lf.writeFile(outPrefix+".PET_count_distribution.txt", str, false);
			    reader = new BufferedReader(new FileReader(df2));
			    line = reader.readLine();
			    while (line != null) {
			    	line = line.trim();
			        String[] strs = line.split("[ \t]+");
			        int strs0 = 0, strs1 = 0, strs2 = 0;
			        if (strs[0].equals("0") && strs[2].equals("0")) {
			        	strs[0] = "10";
			        	strs[1] = "";
			        	strs[2] = "";
			        	strs0 = 10;
			        } else {
			        	strs0 = Integer.valueOf(strs[0]);
			        	strs1 = Integer.valueOf(strs[1]);
			        	strs2 = Integer.valueOf(strs[2]);
			        }
			        String lastColumn = String.valueOf((double)strs2/strs0);
		        	if (strs1 == 10) {
		        		str = ">="+strs[1]+"\t"+strs[0]+"\t"+strs[2]+"\t"+(strs0-strs2)+"\t"+lastColumn;
		        	} else if (strs1 == 11) {
		        		str = "Total\t"+strs[0]+"\t"+strs[2]+"\t"+(strs0-strs2)+"\t"+lastColumn;
		        	} else {
		        		str = strs[1]+"\t"+strs[0]+"\t"+strs[2]+"\t"+(strs0-strs2)+"\t"+lastColumn;
		        	}
		        	lf.writeFile(outPrefix+".PET_count_distribution.txt", str, true);
			        line = reader.readLine();
			    }
			    reader.close();
			    df1.delete();
			    df2.delete();
			} catch (IOException e) {
			    e.printStackTrace();
		    }
		} else {
			System.out.println("Error: "+file+" doesn't exist");
			String str = "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters";
			lf.writeFile(outPrefix+".PET_count_distribution.txt", str, false);
			str = "\t"+"10"+"\t"+"\t"+"10"+"\t"+"0";
			lf.writeFile(outPrefix+".PET_count_distribution.txt", str, true);
			str = "\t"+"11"+"\t"+"\t"+"11"+"\t"+"0";
			lf.writeFile(outPrefix+".PET_count_distribution.txt", str, true);
		}
    }
    
    /**
     * @author sun
     * @function pasteTagCount
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     */
    public void pasteTagCount(String file1, String file2) {
    	File f1 = new File(file1);
    	File f2 = new File(file2);
    	if (f1.exists()) {
    		if (f2.exists()) {
    			try {
	    			BufferedReader reader1 = new BufferedReader(new FileReader(f1));
	    	    	BufferedReader reader2 = new BufferedReader(new FileReader(f2));
	    	    	String line1 = reader1.readLine();
	    	    	String line2 = reader2.readLine();
	    	    	new File(outPrefix+".petCount.tagCount.txt").delete();
	    	    	while (line1 != null && line2 != null) {
	    	    		line1 = line1 + "\t" + line2;
	    	    		line1 = line1.trim();
	    		        String[] strs = line1.split("[ \t]+");
	    		        if (strs.length >= 10) {
	    		        	line1 =  strs[3]+"\t"+strs[4]+"\t"+strs[9];
	    		        	lf.writeFile(outPrefix+".petCount.tagCount.txt", line1, true);
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
     * @function pasteCluster
     * @param file1 the source file to paste
     * @param file2 the source file to paste
     * @param file3 the source file to paste
     * @return the number of line in *.cluster.FDRfiltered.txt
     */
    public void pasteCluster(String file1, String file2, String file3) {
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
		    	    	String line3 = reader3.readLine();
		    	    	new File(outPrefix+".cluster.withpvalue.txt").delete();
		    	    	new File(outPrefix+".cluster.FDRfiltered.txt").delete();
		    	    	while (line1 != null && line2 != null && line3 != null) {
		    	    		line1 = line1 + "\t" + line2 + "\t" + line3;
		    	    		line1 = line1.trim();
		    		        String[] strs = line1.split("[ \t]+");
		    		        if (strs.length >= 24) {
		    		        	line1 = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[6]+"\t"+strs[7]+"\t"+strs[8]+"\t"+strs[12]+"\t"+strs[13]+"\t"+strs[14]+"\t"+
		    		        strs[18]+"\t"+strs[19]+"\t"+strs[20]+"\t"+strs[21]+"\t"+strs[22]+"\t"+strs[23];
		    		        	lf.writeFile(outPrefix+".cluster.withpvalue.txt", line1, true);
		    		        	if (Double.valueOf(strs[21]) < Double.valueOf(p.PVALUE_CUTOFF_INTERACTION)) {
		    		        		lf.writeFile(outPrefix+".cluster.FDRfiltered.txt", line1, true);
		    		        	}
		    		        }else if (strs.length == 20) {
		    		        	int samechrom = 0; int distance = Integer.MAX_VALUE;
		    		        	if(strs[0].equals(strs[6])) {
		    		        		samechrom = 1;
		    		        		distance = Math.abs(Integer.parseInt(strs[8])-Integer.parseInt(strs[1]));
		    		        	}
		    		        	line1 = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[6]+"\t"+strs[7]+"\t"+strs[8]+"\t"+strs[12]+"\t"+samechrom+"\t"+distance+"\t"+
		    		        strs[14]+"\t"+strs[15]+"\t"+strs[16]+"\t"+strs[17]+"\t"+strs[18]+"\t"+strs[19];
		    		        	lf.writeFile(outPrefix+".cluster.withpvalue.txt", line1, true);
		    		        	if (Double.valueOf(strs[17]) < Double.valueOf(p.PVALUE_CUTOFF_INTERACTION)) {
		    		        		lf.writeFile(outPrefix+".cluster.FDRfiltered.txt", line1, true);
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
}
