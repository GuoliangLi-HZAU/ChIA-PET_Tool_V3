package process;

import java.io.File;
import java.io.IOException;
import java.util.Calendar;

public class Main {
	
	public static void main(String []args) throws IOException {
		if (args.length < 22) {
			System.out.println("Error: please set the necessary parameters");
			System.out.println("Usage: java -jar <path of ChIA_PET.jar> [options]");
			System.out.println("Necessary options:");
			System.out.println("    --mode\tmode of tool, 0: short read; 1: long read");
			System.out.println("    --fastq1\tpath of read1 fastq file");
			System.out.println("    --fastq2\tpath of read2 fastq file");
			System.out.println("    --linker\tpath of linker file");
			System.out.println("    --minimum_linker_alignment_score\tminimum alignment score");
			System.out.println("    --GENOME_INDEX\tthe path of BWA index");
			System.out.println("    --GENOME_LENGTH\tthe number of base pairs in the whole genome");
			System.out.println("    --CHROM_SIZE_INFO\tthe file that contains the length of each chromosome, example file is in ChIA-PET_Tool_V3/chrInfo");
			System.out.println("    --CYTOBAND_DATA\tthe ideogram data used to plot intra-chromosomal peaks and interactions, example file is in "
					+ "ChIA-PET_Tool_V3/chrInfo");
			System.out.println("    --SPECIES\t1: human; 2: mouse; 3: others");
			System.out.println("Other options:");
			System.out.println("    --start_step\tstart with which step, 1: linker filtering; 2: mapping to genome; 3: removing redundancy; 4: categorization of "
					+ "PETs; 5: peak calling; 6: interaction calling; 7: visualizing, default: 1");
			System.out.println("    --output\tpath of output, default: ChIA-PET_Tool_V3/output");
			System.out.println("    --prefix\tprefix of output files, default: out");
			System.out.println("    --minimum_tag_length\t minimum tag length, default: 18");
			System.out.println("    --maximum_tag_length\t maximum tag length, default: 1000");
			System.out.println("    --minSecondBestScoreDiff\tthe score difference between the best-aligned and the second-best aligned linkers, default: 3");
			System.out.println("    --output_data_with_ambiguous_linker_info\twhether to print the linker-ambiguity PETs, 0: not print; 1: print, default: 1");
			System.out.println("    --thread\tthe number of threads used in linker filtering and mapping to genome, default: 1");
			System.out.println("    --MAPPING_CUTOFF\tcutoff of mapping quality score for filtering out low-quality or multiply-mapped reads, default: 30");
			System.out.println("    --MERGE_DISTANCE\tthe distance limit to merge the PETs with similar mapping locations, default: 2");
			System.out.println("    --SELF_LIGATION_CUFOFF\tthe distance threshold between self-ligation PETs and intra- chromosomal inter-ligation PETs, "
					+ "default: 8000");
			System.out.println("    --EXTENSION_LENGTH\tthe extension length from the location of each tag, default: 5000");
			System.out.println("    --MIN_COVERAGE_FOR_PEAK\tthe minimum coverage to define peak regions, default: 5");
			System.out.println("    --PEAK_MODE\t1: peak region mode, which takes all the overlapping PET regions above the coverage threshold as peak "
					+ "regions; 2: peak summit mode, which takes the highest coverage of overlapping regions as peak regions, default: 2");
			System.out.println("    --MIN_DISTANCE_BETWEEN_PEAK\tthe minimum distance between two peaks, default: 500");
			System.out.println("    --GENOME_COVERAGE_RATIO\tthe estimated proportion of the genome covered by the reads, default: 0.8");
			System.out.println("    --PVALUE_CUTOFF_PEAK\tp-value to filter peaks that are not statistically significant, default: 0.00001");
			System.out.println("    --INPUT_ANCHOR_FILE\ta file which contains user-specified anchors for interaction calling, default: null");
			System.out.println("    --PVALUE_CUTOFF_INTERACTION\tp-value to filter false positive interactions, default: 0.5");
			System.exit(0);
		}
		
		Path p = new Path();
		p.setParameter(args);
		
		File file = new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX);
		if (!file.exists()) {
			file.mkdirs();
		}
		
		Calendar rightNow = Calendar.getInstance();
		System.out.println("[" + rightNow.getTime().toString() +"] start ChIA-PET analysis");
		LinkerFiltering lf = new LinkerFiltering(p);
		lf.reset(Integer.valueOf(p.START_STEP));
		long time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), false);
		if (Integer.valueOf(p.START_STEP) == 1) {
			System.out.println("[" + rightNow.getTime().toString() +"] Step1: Linker filtering ...");
			lf.run();
		}
		
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 2) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step2: Mapping to genome ...");
			Mapping mapping = new Mapping(p);
			mapping.Map();
		}
		
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 3) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step3: Removing redundancy ...");
			Purifying purifying = new Purifying(p);
			purifying.Purify();
			purifying.combiningData();
		}
		
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 4) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step4: Categorization of PETs ...");
			DividePets dvidePets = new DividePets(p);
			dvidePets.dividePets();
		}
		
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 5) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step5: Interaction Calling ...");
			InteractionCalling interactionCalling = new InteractionCalling(p);
			interactionCalling.run();
			
			Pvalues pValues = new Pvalues(p);
			pValues.calculation();
			pValues.globalTag();
			pValues.calculate();
		}
		
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		PeakCalling peakCalling = new PeakCalling(p);
		if (Integer.valueOf(p.START_STEP) <= 6) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step6: Peak Calling ...");
			peakCalling.run();
			peakCalling.localTag();
			peakCalling.globalSpet();
			peakCalling.calculatePvalue();
			peakCalling.chromosomalPet();
		}
		time = System.currentTimeMillis() / 1000;
		lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		
		if (Integer.valueOf(p.START_STEP) <= 7) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step7: Visualizing ...");
			peakCalling.runningTime();
			peakCalling.summary();
			peakCalling.statisticsReport();
			peakCalling.genomicBrowser();
			new File(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh").delete();
		}
		rightNow = Calendar.getInstance();
		System.out.println("[" + rightNow.getTime().toString() +"] finish ChIA-PET analysis");
	}	
}
