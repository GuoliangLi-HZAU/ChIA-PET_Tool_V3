package process;

import LGL.data.CLUSTER;
import LGL.data.HIT2;
import LGL.data.PET;
import process.LinkerFiltering;
import process.Path;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;

public class Pet2Cluster1 {
	
	private Path p;
	private String outPrefix;
	
	public Pet2Cluster1(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
	public void preCluster() {
		int extensionLength = Integer.parseInt(p.EXTENSION_LENGTH);
		String chromLengthInfo = p.CHROM_SIZE_INFO;
		Hashtable<String, Integer> chromLength = new Hashtable<String, Integer>();
		try {
			BufferedReader fileChromLength = new BufferedReader(new InputStreamReader(new FileInputStream(new File(chromLengthInfo))));
			String lineTemp;
		    while ((lineTemp = fileChromLength.readLine()) != null) {
		        String[] temp = lineTemp.split("[ \t][ \t]*");
		        if (temp.length > 1) {
		            chromLength.put(temp[0], Integer.valueOf(temp[1]));
		        }
		    }
		    fileChromLength.close();
		    
		    int sortingLabel = 1; // ascending order [default]
	        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(outPrefix+".ipet"))));
	        String line;
	        ArrayList<String> tmpPrefixList = new ArrayList<String>();
	        LinkerFiltering lf = new LinkerFiltering(p);
	        Purifying purifying = new Purifying(p);
	        while ((line = fileIn.readLine()) != null) {
	            if (line.length() <= 1) // skip the short lines
	            {
	                continue;
	            }
	            String[] fields = line.split("\t");
	            HIT2 hit1 = new HIT2(fields[0], Integer.parseInt(fields[1]), fields[2].charAt(0));
	            HIT2 hit2 = new HIT2(fields[3], Integer.parseInt(fields[4]), fields[5].charAt(0));
	            double weight = 1.0;
	            if (fields.length >= 7) {
	                weight = Double.parseDouble(fields[6]);
	            }
	            PET newPet = new PET(hit1, hit2, weight);
	            if (fields.length >= 8) {
	                newPet.setIndex(Integer.parseInt(fields[7]));
	            }
	            CLUSTER cluster = new CLUSTER(newPet, extensionLength, sortingLabel, chromLength);
	            if(!(new File(outPrefix+"."+fields[0]+".pre_cluster.txt").exists())) {
	        		tmpPrefixList.add(outPrefix+"."+fields[0]);
	        	}
	        	lf.writeFile(outPrefix+"."+fields[0]+".pre_cluster.txt", cluster.toString(), true);
	        }
	        fileIn.close();
	        
	        BufferedWriter purifysort = new BufferedWriter(new FileWriter(outPrefix +
                    ".pcluster.sort.sh", false));
	        String runsortcmd= "";
	        for (int i = 0; i < tmpPrefixList.size(); i++) {
	        	String tmpPrefix = tmpPrefixList.get(i);
	        	runsortcmd="sh " + p.PROGRAM_DIRECTORY + "/psort-cluster.sh " + tmpPrefix+".pre_cluster.txt " + tmpPrefix+".sort.pre_cluster.txt > /dev/null 2>&1";
	        	purifysort.write(runsortcmd);
                purifysort.newLine();
                runsortcmd="rm " + tmpPrefix+".pre_cluster.txt";
                purifysort.write(runsortcmd);
                purifysort.newLine();
	        	//purifying.sortK(tmpPrefix+".pre_cluster.txt", tmpPrefix+".sort.pre_cluster.txt", new int[]{1, 4, 2});
	        	//new File(tmpPrefix+".pre_cluster.txt").delete();
	        }
	        purifysort.close();
	        Shell shell = new Shell();
            shell.runShell(outPrefix + ".pcluster.sort.sh");
	        
	        purifying.mergeFiles(tmpPrefixList, ".sort.pre_cluster.txt", outPrefix+".pre_cluster.sorted");
		} catch (IOException e) {
		    e.printStackTrace();
	    }
	}
}
