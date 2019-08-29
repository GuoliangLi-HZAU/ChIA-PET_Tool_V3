package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import LGL.chiapet.PetCluster2;
import LGL.chiapet.PetClusterWithGivenAnchors;

public class InteractionCalling {
	
	private Path p;
	private String outPrefix;
	
	public InteractionCalling(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
    public void run() {
    	try {
	    	File file = new File(outPrefix+".ipet");
	    	if (file.exists()) {
	    		Pet2Cluster1 pet2Cluster1 = new Pet2Cluster1(p);
	    		pet2Cluster1.preCluster();
	    	} else {
	    		System.out.println("Error: "+file+" doesn't exist");
	    	}
		    file = new File(outPrefix+".pre_cluster.sorted");
		    if (file.exists()) {
		    	if (!p.INPUT_ANCHOR_FILE.equals("null")) {
		    		PetClusterWithGivenAnchors.main(new String[]{outPrefix+".pre_cluster.sorted", p.INPUT_ANCHOR_FILE, outPrefix+".cluster", "1"});
		    	} else {
		    		PetCluster2.main(new String[]{outPrefix+".pre_cluster.sorted", outPrefix+".cluster"});
		    	}
		    } else {
	    		System.out.println("Error: "+file+" doesn't exist");
	    	}
	    	file = new File(outPrefix+".cluster");
	    	if (file.exists()) {
		        BufferedReader reader = new BufferedReader(new FileReader(file));
		        String line = reader.readLine();
		        new File(outPrefix+".cluster.filtered").delete();
		        new File(outPrefix+".cluster.txt").delete();
			    while (line != null) {
			    	String l = line;
			    	line = line.trim();
			        String[] strs = line.split("[ \t]+");
			        if (strs.length > 15) {
			        	if (Integer.valueOf(strs[12]) >= 2) {
			        		LinkerFiltering lf = new LinkerFiltering(p);
			        		lf.writeFile(outPrefix+".cluster.filtered", l, true);
			        		l = strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[6]+"\t"+strs[7]+"\t"+strs[8]+"\t"+strs[12]+"\t"+strs[13]+"\t"+strs[14];
			        		lf.writeFile(outPrefix+".cluster.txt", l, true);
			        	}
			        }
			       line = reader.readLine();
			    }
			    reader.close();
	    	} else {
	    		System.out.println("Error: "+file+" doesn't exist");
	    	}
    	} catch (IOException e) {
		    e.printStackTrace();
	    }
    }
}
