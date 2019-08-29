package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import LGL.chiapet.MergeSimilarPETs2;

public class Purifying {
	
	private Path p;
	private String outPrefix;
	private LinkerFiltering lf;
	
	public Purifying(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
		lf = new LinkerFiltering(p);
	}
	
    public void Purify() {// purifying the mapped reads
    	File file = new File(outPrefix+".bedpe");
    	if (file.exists()) {
    		try {
    			BufferedReader reader = new BufferedReader(new FileReader(file));
    			String line = reader.readLine();
    			String content = null;
    			int selectNum = 0;
    			ArrayList<String> tmpPrefixList = new ArrayList<String>();
    			new File(outPrefix+".bedpe.selected.txt").delete();
			    while (line != null) {
			    	line = line.trim();
    		        String[] strs = line.split("[ \t]+");
    		        if (strs.length >= 10) {
    		        	content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t.\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
    		        	lf.writeFile(outPrefix+".bedpe.selected.txt", content, true);
    		        	selectNum++;
    		        	if(!(new File(outPrefix+"."+strs[0]+".bedpe.txt").exists())) {
    		        		tmpPrefixList.add(outPrefix+"."+strs[0]);
    		        	}
    		        	lf.writeFile(outPrefix+"."+strs[0]+".bedpe.txt", content, true);
    		        }
	    	        line = reader.readLine();
			    }
			    reader.close();

    			int uniqueNum = 0;
    			Map<Integer, Integer> map2 = new HashMap<Integer, Integer>();
    			Map<String, Integer> map3 = new HashMap<String, Integer>();
    			new File(outPrefix+".bedpe.selected.unique.txt").delete();
    			for (int i = 0; i < tmpPrefixList.size(); i++) {
    				String tmpPrefix = tmpPrefixList.get(i);
    				reader = new BufferedReader(new FileReader(tmpPrefix+".bedpe.txt"));
    				line = reader.readLine();
    				Map<String, Integer> map1 = new HashMap<String, Integer>();
    				while (line != null) {
    					String[] strs = line.split("[ \t]+");
		    	        if (map1.containsKey(line)) {
		    	        	map1.put(line, map1.get(line) + 1);
		    	        } else {
		    	        	map1.put(line, 1);
		    	        	lf.writeFile(outPrefix+".bedpe.selected.unique.txt", line, true);
		    	        	uniqueNum++;
		    	        	if (strs[8].equals("+")) {
		    	        		content = strs[0]+"\t"+strs[2]+"\t"+strs[8]+"\t";
		    	        	} else {
		    	        		content = strs[0]+"\t"+strs[1]+"\t"+strs[8]+"\t";
		    	        	}
		    	        	if (strs[9].equals("+")) {
		    	        		content = content + strs[3]+"\t"+strs[5]+"\t"+strs[9];
		    	        	} else {
		    	        		content = content + strs[3]+"\t"+strs[4]+"\t"+strs[9];
		    	        	}
		    	        	lf.writeFile(tmpPrefix+".purify.temp.txt", content, true);// sort -k
		    	        	if (strs[0].equals(strs[3])) {
		    	        		if (map3.containsKey(strs[8]+"\t"+strs[9])) {
		    	        			map3.put(strs[8]+"\t"+strs[9], map3.get(strs[8]+"\t"+strs[9]) + 1);
		    	        		} else {
		    	        			map3.put(strs[8]+"\t"+strs[9], 1);
		    	        		}
		    	        	}
		    	        }
    					line = reader.readLine();
    				}
    				reader.close();
    				new File(tmpPrefix+".bedpe.txt").delete();
    				Iterator<Map.Entry<String, Integer>> iterator = map1.entrySet().iterator();
    				Map.Entry<String, Integer> entry;
    			    while (iterator.hasNext()) {
    			    	entry = iterator.next();
    			    	Integer count = entry.getValue();
    			    	if (map2.containsKey(count)) {
    			    		map2.put(count, map2.get(count) + 1);
    			    	} else {
    			    		map2.put(count, 1);
    			    	}
    			    }
    			    
    				sortK(tmpPrefix+".purify.temp.txt", tmpPrefix+".bedpe.selected.unique.pet.txt", new int[]{1, 4, 3, 6, 2});
    			    new File(tmpPrefix+".purify.temp.txt").delete();
    			    MergeSimilarPETs2.main(new String[]{tmpPrefix+".bedpe.selected.unique.pet.txt", tmpPrefix+".bedpe.selected.merged.pet.txt", 
    			    		p.MERGE_DISTANCE});
    			    new File(tmpPrefix+".bedpe.selected.unique.pet.txt").delete();
    			}
    			
    			Iterator<Map.Entry<Integer, Integer>> iterator2 = map2.entrySet().iterator();
			    Map.Entry<Integer, Integer> entry2;
			    new File(outPrefix+".bedpe.selected.dist.txt").delete();
			    while (iterator2.hasNext()) {
			    	entry2 = iterator2.next();
			    	Integer k = entry2.getKey();
			    	Integer v = entry2.getValue();
			    	lf.writeFile(outPrefix+".bedpe.selected.dist.txt", v + " " + k, true);
			    }
			    sortN(outPrefix+".bedpe.selected.dist.txt", 1);
			    
			    mergeFiles(tmpPrefixList, ".bedpe.selected.merged.pet.txt", outPrefix+".bedpe.selected.pet.txt");
			    
			    Iterator<Map.Entry<String, Integer>> iterator = map3.entrySet().iterator();
			    Map.Entry<String, Integer> entry;
			    new File(outPrefix+".bedpe.selected.unique.intra-chrom.strand.dist.txt").delete();
			    while (iterator.hasNext()) {
			    	entry = iterator.next();
			    	String strand = entry.getKey();
			    	Integer value = entry.getValue();
			    	lf.writeFile(outPrefix+".bedpe.selected.unique.intra-chrom.strand.dist.txt", value + " " + strand, true);
			    }
			    lf.writeFile(outPrefix+".basic_statistics.txt", "Uniquely Mapped same-linker PETs\t"+selectNum, true);
			    lf.writeFile(outPrefix+".basic_statistics.txt", "Merging same same-linker PETs\t"+uniqueNum, true);
    		} catch (IOException e) {
			    e.printStackTrace();
		    }
    	} else {
    		System.out.println("Error: "+file+" doesn't exist");
    	}
    }

    /**
     * @author sun
     * @function sort by column
     * @param source source file
     * @param dest dest file
     * @param l column num array
     */
    @SuppressWarnings("unchecked")
	public void sortK(String source, String dest, int []l) {
		File f = new File(source);
    	if (f.exists()) {
	        try {
	        	BufferedReader reader = new BufferedReader(new FileReader(f));
	        	String line = reader.readLine();
	        	ArrayList<SortK> list = new ArrayList<SortK>();
		        while (line != null) {
		        	SortK sortk = new SortK();
		        	sortk.setContent(line);
		        	line = line.trim();
		        	line = line.replace(" ", "\t");
    		        String[] strs = line.split("\t");
    		        for (int i = 0; i < l.length; i++) {
    		        	sortk.setColumns(strs[l[i] - 1]);
    		        }
		        	list.add(sortk);
		        	line = reader.readLine();
		        }
		        reader.close();
		        Collections.sort(list);
		        for (int i = 0; i < list.size(); i++) {
		        	lf.writeFile(dest, list.get(i).getContent(), true);
		        }
	        } catch (IOException e) {
			    e.printStackTrace();
		    }
    	} else {
    		System.out.println("Error: "+source+" doesn't exist");
    	}
	}
    
    /**
     * @author sun
     * @function sort by num of which column
     * @param file source file
     * @param n which column
     */
    @SuppressWarnings("unchecked")
	public void sortN(String file, int n) {
		File f = new File(file);
		if (f.exists()) {
	        try {
	        	BufferedReader reader = new BufferedReader(new FileReader(f));
	        	String line = reader.readLine();
	        	ArrayList<SortN> list = new ArrayList<SortN>();
		        while (line != null) {
		        	SortN sortn = new SortN();
		        	sortn.setContent(line);
		        	line = line.trim();
		        	line = line.replace(" ", "\t");
		        	String[] strs = line.split("\t");
			        sortn.setColumn(Integer.valueOf(strs[n - 1]));
			        list.add(sortn);
		        	line = reader.readLine();
		        }
		        reader.close();
		        Collections.sort(list);
		        new File(file).delete();
		        for (int i = 0; i < list.size(); i++) {
		        	lf.writeFile(file, list.get(i).getContent(), true);
		        }
	        } catch (IOException e) {
			    e.printStackTrace();
		    }
		} else {
			System.out.println("Error: "+file+" doesn't exist");
		}
	}
    
    public void mergeFiles(ArrayList<String> fileList, String tailName, String file) {
        try {
        	FileOutputStream fos = new FileOutputStream(file);
        	FileChannel destfc = fos.getChannel();
        	for (int i = 0; i < fileList.size(); i++) {
        		String f = fileList.get(i) + tailName;
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
	                srcFile.delete();
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
    
    public void combiningData() {
    	String file = outPrefix+".bedpe.selected.pet.txt";
    	int lineNum = spanDistribution(file);
    	String str = "Merging similar same-linker PETs\t";
    	lf.writeFile(outPrefix+".basic_statistics.txt", str+lineNum, true);
    }
    
    /**
	 * @author sun
	 * @function spanDistribution
	 * @param file source file
	 * @param lineNum count
	 * @return the number of line of source file
	 */
    public int spanDistribution(String file) {
    	int lineNum = 0;
    	File f = new File(file);
    	if (f.exists()) {
    	    try {
    	        BufferedReader reader = new BufferedReader(new FileReader(f));
    	        String line = reader.readLine();
	        	String destFile = outPrefix+".bedpe.selected.intra-chrom.distance.txt";
	        	new File(outPrefix+".bedpe.selected.intra-chrom.distance.txt").delete();
	        	new File(outPrefix+".bedpe.selected.intra-chrom.distance.plusplus.txt").delete();
	        	new File(outPrefix+".bedpe.selected.intra-chrom.distance.plusminus.txt").delete();
	        	new File(outPrefix+".bedpe.selected.intra-chrom.distance.minusminus.txt").delete();
			    while (line != null) {
			    	lineNum++;
			    	line = line.trim();
    		        String[] strs = line.split("[ \t]+");
    		        if (strs.length >= 6) {
    		        	if (strs[0].equals(strs[3])) {
    		            	line = String.valueOf(Integer.valueOf(strs[4])-Integer.valueOf(strs[1]));
    		            	String l = line+"\t"+strs[2]+"\t"+strs[5];   
    		    	        lf.writeFile(destFile, l, true);
    		    	        if (strs[2].equals("+") && strs[5].equals("+")) {
 	    		        	   lf.writeFile(outPrefix+".bedpe.selected.intra-chrom.distance.plusplus.txt", line, true);
 	    		            } else if (strs[2].equals("+") && strs[5].equals("-")) {
 	    		        	    lf.writeFile(outPrefix+".bedpe.selected.intra-chrom.distance.plusminus.txt", line, true);
 	    		            } else if (strs[2].equals("-") && strs[5].equals("+")) {
 	    		        	    lf.writeFile(outPrefix+".bedpe.selected.intra-chrom.distance.minusplus.txt", line, true);
 	    		            } else if (strs[2].equals("-") && strs[5].equals("-")) {
 	    		        	    lf.writeFile(outPrefix+".bedpe.selected.intra-chrom.distance.minusminus.txt", line, true);
 	    		            }
    		            } 
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
    	return lineNum;
    }
}
