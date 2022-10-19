package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import LGL.chiapet.PetCluster2;
import LGL.chiapet.PetClusterWithGivenAnchors;

public class findpeak {
	
	private Path p;
	private String outPrefix;
	
	public findpeak(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
    public void run() {
    	try {
	    	File file = new File(outPrefix+".ipet");
	    	if (file.exists()) {
	    		//ipet 2 .pre_cluster.sorted
	    		Pet2Cluster1 pet2Cluster1 = new Pet2Cluster1(p);
	    		pet2Cluster1.preCluster();
	    	} else {
	    		System.out.println("Error: "+file+" doesn't exist");
	    	}
		    file = new File(outPrefix+".pre_cluster.sorted");
		    if (file.exists()) {
		    	PetCluster2.main(new String[]{outPrefix+".pre_cluster.sorted", outPrefix+".cluster"});
		    } else {
	    		System.out.println("Error: "+file+" doesn't exist");
	    	}
	    	file = new File(outPrefix+".cluster");
	    	if (file.exists()) {
		        BufferedReader reader = new BufferedReader(new FileReader(file));
		        String line = reader.readLine();
		        //new File(outPrefix+".cluster.filtered").delete();
		        new File(outPrefix+".panchor.bed").delete();
			    while (line != null) {
			    	String l = line;
			    	line = line.trim();
			        String[] strs = line.split("[ \t]+");
			        if (strs.length > 15) {
			        	if (Integer.valueOf(strs[12]) >= 2) {
			        		LinkerFiltering lf = new LinkerFiltering(p);
			        		//lf.writeFile(outPrefix+".cluster.filtered", l, true);
			        		l = strs[0]+"\t"+strs[1]+"\t"+strs[2];
			        		lf.writeFile(outPrefix+".panchor.bed", l, true);
			        		l = strs[6]+"\t"+strs[7]+"\t"+strs[8]; //+"\t"+strs[12];
			        		lf.writeFile(outPrefix+".panchor.bed", l, true);
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
