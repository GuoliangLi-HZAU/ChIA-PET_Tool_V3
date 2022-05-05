package process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import LGL.dataConcise.RegionConcise;

import LGL.chiapet.MergeSimilarPETs2;
import LGL.data.REGION;
import LGL.shortReads.TagCountInGivenRegions;

public class Purifying {
	
	private Path p;
	private String outPrefix;
	private String yesHiChIP;
	private LinkerFiltering lf;
	
	public Purifying(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
		lf = new LinkerFiltering(p);
		yesHiChIP = p.hichipM;
	}
	
	int[] getResSite(String newchr, String nowchr, Hashtable<String, ArrayList<RegionConcise>> hashChrom2regionConcise,
			int start, int end, Vector<REGION> regions, Hashtable<String, Integer> hashChrom2maxRegionSpans) {
		ArrayList regionConciseList = null;
		int[] finalR = new int[3];
		finalR[0] = -1;
		finalR[1] = -1;
		finalR[2] = -1;
		if(newchr.equals("") || !newchr.equals(nowchr)) {
        	newchr = nowchr;
        	regionConciseList = hashChrom2regionConcise.get(nowchr);
        }

        if (regionConciseList == null) {
        	System.out.println("Unexpected error: pet chr is not found in restriction site file!!!" + 
               nowchr + " Please check your restriction site file!!!");
        	finalR[0] = -1;
            return finalR;
        }
       // int start = Integer.parseInt(fields[1]), end = Integer.parseInt(fields[2]);
        RegionConcise regionConcise = new RegionConcise(start, end);
        //RegionConcise regionConciseR = new RegionConcise();
        
        int index = Collections.binarySearch(regionConciseList, regionConcise);
        if (index > 0) { 
        	finalR[0] = index;
        	finalR[1] = ((ArrayList<RegionConcise>) regionConciseList).get(index).getStart();
        	finalR[2] = ((ArrayList<RegionConcise>) regionConciseList).get(index).getEnd();
            return finalR;
        }
        
        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regionConciseList.size()) {
            index = regionConciseList.size() - 1;
        }

        for (int iRegion = index; iRegion < regionConciseList.size(); iRegion++) {
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) > 0) {
            	finalR[0] = iRegion;
            	finalR[1] = ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getStart();
            	finalR[2] = ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getEnd();
                return finalR;
            } else {
                break;
            }
        }
        
        // exhaustive search, in case that there is an extreme region to cover the whole chromosome
        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (start - ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getStart() > hashChrom2maxRegionSpans.get(newchr).intValue()) {
                break;
            }
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) > 0) {
            	finalR[0] = iRegion;
            	finalR[1] = ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getStart();
            	finalR[2] = ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getEnd();
                return finalR;
            }
        }
        
        return finalR;
	}
	
	public void removePETinsameblock(String bedpefile, String filterbedpefile, String sameresbedpefile) throws IOException {
		Hashtable<RegionConcise, REGION> hashRegionConcise2Region = new Hashtable<RegionConcise, REGION>();
	    Hashtable<String, Integer> hashChrom2maxRegionSpans = new Hashtable<String, Integer>();
	    Hashtable<String, ArrayList<RegionConcise>> hashChrom2regionConcise = new Hashtable<String, ArrayList<RegionConcise>>();;
		//Vector<REGION> regions = REGION.load(p.restrictionsiteFile);
	    Vector<REGION> regions = REGION.loadwithfilter(p.restrictionsiteFile, p.minfragsize, p.maxfragsize);
		TagCountInGivenRegions.generateRegionConciseHash(regions, hashChrom2regionConcise, 
				hashChrom2maxRegionSpans, hashRegionConcise2Region);
		
		//
		BufferedWriter filterbedpeBufferedWriter = new BufferedWriter(new FileWriter(filterbedpefile));
		BufferedWriter sameresbedpeBufferedWriter = new BufferedWriter(new FileWriter(sameresbedpefile));
        BufferedReader bedpefileP = new BufferedReader(
        		new InputStreamReader(new FileInputStream(new File( bedpefile ))));
        String line;
        String newchr = "";
        int ResSite1 = -1, ResSite2 = -1;
        String newline = "";
        int[] finalR1 = new int[3];
        int[] finalR2 = new int[3];
        int fragmentsize = 0, frag1 = 0, frag2 = 0;
        while ((line = bedpefileP.readLine()) != null) {
            if (line.length() <= 0) {  // skip the short lines
                continue;
            }
            String fields[] = line.split("\t");
            finalR1 = getResSite(newchr, fields[0], hashChrom2regionConcise,
            		Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), regions, hashChrom2maxRegionSpans);
            
            finalR2 = getResSite(newchr, fields[3], hashChrom2regionConcise,
            		Integer.parseInt(fields[4]), Integer.parseInt(fields[5]), regions, hashChrom2maxRegionSpans);
            
            ResSite1 = finalR1[0];
            ResSite2 = finalR2[0];
            if(ResSite1 != ResSite2 && ResSite1>=0 && ResSite2>=0) {
            	if(fields[8].equals("+")) {
            		fragmentsize = finalR1[2] - Integer.parseInt(fields[1]);
            	}else {
            		fragmentsize = Integer.parseInt(fields[2]) - finalR1[1];
            	}
            	if(fields[9].equals("+")) {
            		fragmentsize += (finalR2[2] - Integer.parseInt(fields[4]));
            	}else {
            		fragmentsize += (Integer.parseInt(fields[5]) - finalR2[1]);
            	}
            	/*
            	if(fragmentsize < p.minInsertsize || fragmentsize > p.maxInsertsize) {
            		continue;
            	}
            	*/
                newline = line + "\t" + ResSite1 + "\t" + ResSite2 + "\t" + fragmentsize;
                filterbedpeBufferedWriter.write(newline);
                filterbedpeBufferedWriter.newLine();
            }else {
            	newline = line + "\t" + ResSite1 + "\t" + ResSite2 + "\tN";
            	sameresbedpeBufferedWriter.write(newline);
            	sameresbedpeBufferedWriter.newLine();
            }
            
        }
        bedpefileP.close();
        filterbedpeBufferedWriter.close();
        sameresbedpeBufferedWriter.close();
	}
	
    public void Purify() throws IOException {// purifying the mapped reads
    	
		String bedpefile = outPrefix+".bedpe";
		if(p.hichipM.equals("Y")) {
			//remove and filter pet in same restriction block
			if(p.removeResblock.equals("Y")) {
				removePETinsameblock(outPrefix+".bedpe", outPrefix+".bedpe.filter.byres", outPrefix+".bedpe.insameres");
				bedpefile = outPrefix+".bedpe.filter.byres";
			}
		}
    	File file = new File(bedpefile);
    	if (file.exists()) {
    		try {
    			BufferedReader reader = new BufferedReader(new FileReader(file));
    			String line = reader.readLine();
    			String content = null;
    			int selectNum = 0;
    			ArrayList<String> tmpPrefixList = new ArrayList<String>();
    			new File(outPrefix+".bedpe.selected.txt").delete();
    			BufferedWriter bedpeselectBufferedWriter = new BufferedWriter(new FileWriter(outPrefix+".bedpe.selected.txt"));
    			BufferedWriter[] bedpeChrBufferedWriter = new BufferedWriter[p.Ngenome];
    			for(int k=0; k<p.Ngenome; k++) {
    				//new File(outPrefix+"." + p.chrMAP_r.get(k) + ".bedpe.txt").delete();
    				bedpeChrBufferedWriter[k] =  new BufferedWriter(new FileWriter(outPrefix+"." + p.chrMAP_r.get(k) + ".bedpe.txt"));
    			}
			    while (line != null) {
			    	line = line.trim();
    		        String[] strs = line.split("[ \t]+");
    		        if (strs.length >= 10) {
    		        	//System.out.println("strs[6] " + strs[6] + " " + p.printreadID);
    		        	if(p.printreadID.equals("Y")) { 
    		        		//System.out.println("AAA ---- ");
    		        		content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t" + strs[6] + "\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
    		        	}else
    		        	    content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t.\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
    		        	//System.out.println("AAA "+ content);
    		        	//lf.writeFile(outPrefix+".bedpe.selected.txt", content, true);
    		        	bedpeselectBufferedWriter.write(content);
    		        	bedpeselectBufferedWriter.newLine();
    		        	selectNum++;
    		        	Collections.sort(tmpPrefixList);
    		        	if(Collections.binarySearch(tmpPrefixList, outPrefix+"."+strs[0])<0) {
    		        		tmpPrefixList.add(outPrefix+"."+strs[0]);
    		        	}
    		        	//lf.writeFile(outPrefix+"."+strs[0]+".bedpe.txt", content, true);
    		        	if(p.chrMAP.get(strs[0]) == null) {
    		                System.out.println("Error, not found " + strs[0] + " in --CHROM_SIZE_INFO file, please make sure your CHROM SIZE file chromosome name "
    		                		+ "\nis same as chromosome name in your genome file !!!");
    		                System.exit(0);
    		            }
    		        	int k = p.chrMAP.get(strs[0]);
    		        	bedpeChrBufferedWriter[k].write(content);
    		        	bedpeChrBufferedWriter[k].newLine();
    		        }
	    	        line = reader.readLine();
			    }
			    reader.close();
			    bedpeselectBufferedWriter.close();
			    for(int k=0; k<p.Ngenome; k++) {
    				bedpeChrBufferedWriter[k].close();
    			}

    			int uniqueNum = 0;
    			Map<Integer, Integer> map2 = new HashMap<Integer, Integer>();
    			Map<String, Integer> map3 = new HashMap<String, Integer>();
    			new File(outPrefix+".bedpe.selected.unique.txt").delete();
    			BufferedWriter purifysort = new BufferedWriter(new FileWriter(outPrefix + 
    					".purify.sort.sh", false));
    			
    			for (int i = 0; i < tmpPrefixList.size(); i++) {
    				String tmpPrefix = tmpPrefixList.get(i);
    				
    				//process bedpe file and get uniq pet
    				reader = new BufferedReader(new FileReader(tmpPrefix+".bedpe.txt"));
    				line = reader.readLine();
    				Map<String, Integer> map1 = new HashMap<String, Integer>();
    				Map<String, String> mapID = new HashMap<String, String>();
    				while (line != null) {
    					String[] strs = line.split("[ \t]+");
    					if(p.printreadID.equals("Y")) {
    		        		//System.out.println("AAA ---- ");
    		        		content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t.\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
    		        	}else
    		        		content = line;
		    	        if (map1.containsKey(content)) {
		    	        	map1.put(content, map1.get(content) + 1);
		    	        	if(p.printreadID.equals("Y")) {
		    	        	    mapID.put(content, mapID.get(content) + "," + strs[6]);
		    	        	}
		    	        } else {
		    	        	map1.put(content, 1);
		    	        	if(p.printreadID.equals("Y")) mapID.put(content, strs[6]);
		    	        	if(!p.printreadID.equals("Y")) {
		    	        		//content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t" + strs[6] + "\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
		    	        		lf.writeFile(outPrefix+".bedpe.selected.unique.txt", line, true);
		    	        	}
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
    				
    				// write uniq pet with ID
    				Iterator<Map.Entry<String, String>> iteratorID = mapID.entrySet().iterator();
    				Map.Entry<String, String> entryID;
    			    while (iteratorID.hasNext()) {
    			    	entryID = iteratorID.next();
    			    	String readsID = entryID.getValue();
    			    	String[] strs = entryID.getKey().split("[ \t]+");
    			    	content = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t" + readsID + "\t"+p.MAPPING_CUTOFF+"\t"+strs[8]+"\t"+strs[9];
    	        		lf.writeFile(outPrefix+".bedpe.selected.unique.txt", content, true);
    			    }
    			    //end
    				
    				//new File(tmpPrefix+".bedpe.txt").delete();
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
    			    //System.out.println("\nI think here is normal case!!!!\n");
    			    
    			    //using system sort
    			    String runsortcmd="sort -k1,1 -k4,4 -k3,3r -k6,6r -k2,2n -k5,5n "+tmpPrefix+".purify.temp.txt > " + tmpPrefix+".bedpe.selected.unique.pet.txt";
    			    purifysort.write(runsortcmd);
    			    purifysort.newLine();
    				//sortK(tmpPrefix+".purify.temp.txt", tmpPrefix+".bedpe.selected.unique.pet.txt", new int[]{1, 4, 3, 6, 2});
    			    //new File(tmpPrefix+".purify.temp.txt").delete();
    			    //MergeSimilarPETs2.main(new String[]{tmpPrefix+".bedpe.selected.unique.pet.txt", tmpPrefix+".bedpe.selected.merged.pet.txt", 
    			    //		p.MERGE_DISTANCE});
    			    //new File(tmpPrefix+".bedpe.selected.unique.pet.txt").delete();
    			}// end tmpPrefix
    			
    			// rm sort shell and delete tempfile
    			//String outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
    			purifysort.close();
    			Shell shell = new Shell();
    			shell.runShell(outPrefix + ".purify.sort.sh");
    			for (int i = 0; i < tmpPrefixList.size(); i++) {
    				String tmpPrefix = tmpPrefixList.get(i);
    				MergeSimilarPETs2.main(new String[]{tmpPrefix+".bedpe.selected.unique.pet.txt", tmpPrefix+".bedpe.selected.merged.pet.txt", 
    	    			    		p.MERGE_DISTANCE});
    				new File(tmpPrefix+".purify.temp.txt").delete();
    				new File(tmpPrefix+".bedpe.selected.unique.pet.txt").delete();
    			}
    			new File(outPrefix+".purify.sort.sh").delete();
    			
    			
    			//delete temp chrom bedpe file
    			for(int k=0; k<p.Ngenome; k++) {
    				new File(outPrefix+"." + p.chrMAP_r.get(k) + ".bedpe.txt").delete();
    				//bedpeChrBufferedWriter[k] =  new BufferedWriter(new FileWriter(outPrefix+"." + p.chrMAP_r.get(k) + ".bedpe.txt"));
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
			    if(p.hichipM.equals("Y") || p.hichipM.equals("O")) {
			    	lf.writeFile(outPrefix+".basic_statistics.txt", "Uniquely Mapped PETs\t"+selectNum, true);
			    	lf.writeFile(outPrefix+".basic_statistics.txt", "Merging same PETs\t"+uniqueNum, true);
			    }else {
			    	lf.writeFile(outPrefix+".basic_statistics.txt", "Uniquely Mapped same-linker PETs\t"+selectNum, true);
			    	lf.writeFile(outPrefix+".basic_statistics.txt", "Merging same same-linker PETs\t"+uniqueNum, true);
			    }
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
    	if(!yesHiChIP.equalsIgnoreCase("Y")) {
    		String str = "Merging similar same-linker PETs\t";
    		lf.writeFile(outPrefix+".basic_statistics.txt", str+lineNum, true);
    	}
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
