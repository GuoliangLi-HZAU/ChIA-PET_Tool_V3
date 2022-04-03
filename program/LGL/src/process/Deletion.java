package process;

import java.io.File;

import process.Path;

public class Deletion {
	
	private Path p;
	private String outPrefix;
	
	public Deletion(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
	}
	
    public void move() {// move files
    	File f = new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report");
    	if (!f.exists()) {
    		f.mkdir();
    	}	
		String path1 = outPrefix;
		String path2 = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/"+p.OUTPUT_PREFIX;
		String[] name1 = {".runningInformation.LinkerFiltering_FastQ_PET.txt", ".tag_length_distribution.txt", ".linker_alignment_score_distribution.txt", 
				".linker_alignment_score_difference_distribution.txt", ".runningInformation.Mapping.txt", ".runningInformation.RemovingRedundancy.txt", 
				".runningInformation.Categories.txt",".runningInformation.clustering.txt", ".runningInformation.Peakcalling.txt", 
				".runningInformation.Report.txt"};
		Pvalues pvalues = new Pvalues(p);
		for (int i = 0; i < name1.length; i++) {
			f = new File(path1 + name1[i]);
			if (f.exists()) {
    			pvalues.copyFile(f, new File(path2 + name1[i]));
        		f.delete();
			}
		}
		String[] name2 = {".bedpe.qc.dist.txt", ".bedpe.strand.dist.txt",".bedpe.selected.dist.txt", ".bedpe.selected.unique.intra-chrom.strand.dist.txt"};
		if (p.MODE.equals("0") || p.hichipM.equals("Y")) {
			String map = ".mapping_statistics.txt";
			f = new File(path1+map);
			if (f.exists()) {
				pvalues.copyFile(f, new File(path2+map));
				f.delete();
			}
		} else {
			String[] prefix = {".1_2", ".2_1"};
			for (int i = 0; i < prefix.length; i++) {
				String map1 = ".1_2.mapping_statistics.txt";
				String map2 = ".3_4.mapping_statistics.txt";
				String map3 = ".5_6.mapping_statistics.txt";
				String map4 = ".7_8.mapping_statistics.txt";
				f = new File(path1+prefix[i]+map1);
				if (f.exists()) {
					pvalues.copyFile(f, new File(path2+prefix[i]+map1));
					f.delete();
				}
    			f = new File(path1+prefix[i]+map2);
    			if (f.exists()) {
    				pvalues.copyFile(f, new File(path2+prefix[i]+map2));
    				f.delete();
    			}
    			f = new File(path1+prefix[i]+map3);
    			if (f.exists()) {
    				pvalues.copyFile(f, new File(path2+prefix[i]+map3));
    				f.delete();
    			}
    			f = new File(path1+prefix[i]+map4);
    			if (f.exists()) {
    				pvalues.copyFile(f, new File(path2+prefix[i]+map4));
    				f.delete();
    			}
			}
		}
		for (int j = 0; j < name2.length; j++) {
			f = new File(path1+name2[j]);
			if (f.exists()) {
    			pvalues.copyFile(f, new File(path2+name2[j]));
        		f.delete();
			}
		}
		String[] name3 = {".bedpe.selected.intra-chrom.distance.txt", ".bedpe.selected.intra-chrom.distance.plusplus.txt", 
				".bedpe.selected.intra-chrom.distance.plusminus.txt", ".bedpe.selected.intra-chrom.distance.minusplus.txt", 
				".bedpe.selected.intra-chrom.distance.minusminus.txt", ".PET_count_distribution.txt", ".basic_statistics.txt",".running_time.txt",
				".summary.txt"};
		for (int i = 0; i < name3.length; i++) {
			f = new File(path1+name3[i]);
			if (f.exists()) {
    			pvalues.copyFile(f, new File(path2+name3[i]));
        		f.delete();
			}
		}	
		copy();
		remove();
    }
    
    public void copy () {// copy files
    	Pvalues pvalues = new Pvalues(p);
    	File f = new File(outPrefix+".linker_composition_distribution.txt");
    	if (f.exists()) {
    		pvalues.copyFile(f, new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/"+p.OUTPUT_PREFIX+".linker_composition_distribution.txt"));
    	}
    	f = new File(outPrefix+".cluster.FDRfiltered.txt");
    	if (f.exists()) {
    		pvalues.copyFile(f, new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/"+p.OUTPUT_PREFIX+".cluster.FDRfiltered.txt"));
    	}
    	f = new File(outPrefix+".peak.FDRfiltered.txt");
    	if (f.exists()) {
    		pvalues.copyFile(f, new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/files_for_report/"+p.OUTPUT_PREFIX+".peak.FDRfiltered.txt"));
    	}
    }
    
    public void remove () {// delete files
    	new File(outPrefix+".pre_cluster.sorted").delete();
    	new File(outPrefix+".cluster.filtered").delete();
    	new File(outPrefix+".cluster.txt").delete();
    	new File(outPrefix+".aln").delete();
    	new File(outPrefix+".cluster.filtered.anchor1").delete();
    	new File(outPrefix+".cluster.filtered.anchor2").delete();
    	new File(outPrefix+".cluster.filtered.anchor1.tagCount.txt").delete();
    	new File(outPrefix+".cluster.filtered.anchor2.tagCount.txt").delete();
    	new File(outPrefix+".nTags.txt").delete();
    	new File(outPrefix+".petCount.tagCount.txt").delete();
    	new File(outPrefix+".pvalue.hypergeo.txt").delete();
    	new File(outPrefix+".peak.5K_5K").delete();
    	new File(outPrefix+".peak.10K_10K").delete();
    	new File(outPrefix+".peak.spetCounts.5K_5K").delete();
    	new File(outPrefix+".peak.spetCounts.10K_10K").delete();
    	new File(outPrefix+".nSpets.txt").delete();
    	new File(outPrefix+".peak.spetCounts.txt").delete();
    	new File(outPrefix+".pvalue.pois.txt").delete();
    	new File(outPrefix+".peak.withTagCounts").delete();
    	new File(outPrefix+".temp.txt").delete();
    	new File(outPrefix+".peak.aln").delete();
    	
    	//new File(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".mapping.sh").delete();
    	//new File(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh").delete();
    	//new File(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".pvalues.sh").delete();
    	new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".mapping.sh").delete();
    	new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling1.sh").delete();
    	new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".pvalues.sh").delete();
    }
}
