package process;

import java.io.File;

public class Mapping {

	private Path p;
	private String outPrefix;
	private LinkerFiltering lf;
	
	public Mapping(Path path) {
		p = path;
		outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
		lf = new LinkerFiltering(p);
	}
	
    public void Map() {// mapping reads to a reference genome
    	String file = p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".mapping.sh";
    	new File(file).delete();
    	if (p.MODE.equals("0")) {
    		mapShortRead(file);
    	} else {
    		mapLongRead(file);
    	}
        Shell shell = new Shell();
        shell.runShell(file);
    }
    
    public void mapShortRead(String file) {
    	String line = "for x in R1 R2";
    	lf.writeFile(file, line, true);
    	line = "do";
    	lf.writeFile(file, line, true);
    	line = "    cat "+outPrefix+".1_1.${x}.fastq "+outPrefix+".2_2.${x}.fastq > "+outPrefix+".${x}.fastq";
    	lf.writeFile(file, line, true);
    	line = "    bwa aln -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.fastq 1> "+outPrefix+".${x}.fastq.sai 2> "+outPrefix+
    			".${x}.sai.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa samse "+p.GENOME_INDEX+" "+outPrefix+".${x}.fastq.sai "+outPrefix+".${x}.fastq 1> "+outPrefix+".${x}.sam 2> "+outPrefix+
    			".${x}.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    rm "+outPrefix+".${x}.fastq";
    	lf.writeFile(file, line, true);
    	line = "    rm "+outPrefix+".${x}.fastq.sai";
    	lf.writeFile(file, line, true);
    	line = "done";
    	lf.writeFile(file, line, true);
    	line = "java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.MappingStatistics "+outPrefix+".R1.sam "+outPrefix+".R2.sam "+outPrefix+" "+
    	p.MAPPING_CUTOFF+" "+p.NTHREADS;
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+".R1.sam";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+".R2.sam";
    	lf.writeFile(file, line, true);
    	line = "cut -f8    < "+outPrefix+".bedpe | LANG=C sort -n | uniq -c > "+outPrefix+".bedpe.qc.dist.txt";
    	lf.writeFile(file, line, true);
    	line = "cut -f9,10 < "+outPrefix+".bedpe | LANG=C sort -n | uniq -c > "+outPrefix+".bedpe.strand.dist.txt";
    	lf.writeFile(file, line, true);
    }
    
    public void mapLongRead(String file) {
    	String line = "for x in  1_2 2_1";
    	lf.writeFile(file, line, true);
    	line = "do";
    	lf.writeFile(file, line, true);
    	line = "    bwa aln -n 2 -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R1.fastq 1> "+outPrefix+".${x}.R1.sai 2>"+outPrefix+
    			".${x}.R1.sai.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa samse "+p.GENOME_INDEX+" "+outPrefix+".${x}.R1.sai "+outPrefix+".${x}.R1.fastq 1>"+outPrefix+".${x}.R1.sam 2>"+outPrefix+
    			".${x}.R1.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa aln -n 2 -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R2.fastq 1> "+outPrefix+".${x}.R2.sai 2>"+outPrefix+
    			".${x}.R2.sai.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa samse "+p.GENOME_INDEX+" "+outPrefix+".${x}.R2.sai "+outPrefix+".${x}.R2.fastq 1>"+outPrefix+".${x}.R2.sam 2>"+outPrefix+
    			".${x}.R2.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa aln -n 2 -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R3.fastq 1> "+outPrefix+".${x}.R3.sai 2>"+outPrefix+
    			".${x}.R3.sai.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa samse "+p.GENOME_INDEX+" "+outPrefix+".${x}.R3.sai "+outPrefix+".${x}.R3.fastq 1>"+outPrefix+".${x}.R3.sam 2>"+outPrefix+
    			".${x}.R3.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa mem -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R4.fastq 1>"+outPrefix+".${x}.R4.sam 2>"+outPrefix+
    			".${x}.R4.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa mem -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R5.fastq 1>"+outPrefix+".${x}.R5.sam 2>"+outPrefix+
    			".${x}.R5.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa aln -n 2 -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R6.fastq 1> "+outPrefix+".${x}.R6.sai 2>"+outPrefix+
    			".${x}.R6.sai.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa samse "+p.GENOME_INDEX+" "+outPrefix+".${x}.R6.sai "+outPrefix+".${x}.R6.fastq 1>"+outPrefix+".${x}.R6.sam 2>"+outPrefix+
    			".${x}.R6.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa mem -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R7.fastq 1>"+outPrefix+".${x}.R7.sam 2>"+outPrefix+
    			".${x}.R7.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    bwa mem -t "+p.NTHREADS+" "+p.GENOME_INDEX+" "+outPrefix+".${x}.R8.fastq 1>"+outPrefix+".${x}.R8.sam 2>"+outPrefix+
    			".${x}.R8.sam.output.info.txt";
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.UniqueSam "+outPrefix+".${x}.R4.sam "+outPrefix+".${x}.R4.clean.sam";
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.UniqueSam "+outPrefix+".${x}.R5.sam "+outPrefix+".${x}.R5.clean.sam";
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.UniqueSam "+outPrefix+".${x}.R7.sam "+outPrefix+".${x}.R7.clean.sam";
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.UniqueSam "+outPrefix+".${x}.R8.sam "+outPrefix+".${x}.R8.clean.sam";
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.MappingStat "+outPrefix+".${x}.R1.sam "+outPrefix+".${x}.R2.sam "+outPrefix+
    			".${x}.1_2 "+p.MAPPING_CUTOFF;
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.MappingStat "+outPrefix+".${x}.R3.sam "+outPrefix+".${x}.R4.clean.sam "+outPrefix+
    			".${x}.3_4 "+p.MAPPING_CUTOFF;
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.MappingStat "+outPrefix+".${x}.R5.clean.sam "+outPrefix+".${x}.R6.sam "+outPrefix+
    			".${x}.5_6 "+p.MAPPING_CUTOFF;
    	lf.writeFile(file, line, true);
    	line = "    java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.MappingStat "+outPrefix+".${x}.R7.clean.sam "+outPrefix+".${x}.R8.clean.sam "+
    	outPrefix+".${x}.7_8 "+p.MAPPING_CUTOFF;
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.1_2.bedpe.selected.uniq.1.sam |bamToBed -i > "+outPrefix+".${x}.R1.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.1_2.bedpe.selected.uniq.2.sam |bamToBed -i > "+outPrefix+".${x}.R2.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.3_4.bedpe.selected.uniq.1.sam |bamToBed -i > "+outPrefix+".${x}.R3.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.3_4.bedpe.selected.uniq.2.sam |bamToBed -i > "+outPrefix+".${x}.R4.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.5_6.bedpe.selected.uniq.1.sam |bamToBed -i > "+outPrefix+".${x}.R5.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.5_6.bedpe.selected.uniq.2.sam |bamToBed -i > "+outPrefix+".${x}.R6.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.7_8.bedpe.selected.uniq.1.sam |bamToBed -i > "+outPrefix+".${x}.R7.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    samtools view -Sb "+outPrefix+".${x}.7_8.bedpe.selected.uniq.2.sam |bamToBed -i > "+outPrefix+".${x}.R8.uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "    paste "+outPrefix+".${x}.R1.uniq.bed "+outPrefix+".${x}.R2.uniq.bed |"+
    	"awk '{score=$5;if(score>$11){score=$11};if(($1==$7 && $2<$8) || ($1<$7))"+
    			"{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$4\"\\t\"score\"\\t\"$6\"\\t\"$12}"+
    			"else{print $7\"\\t\"$8\"\\t\"$9\"\\t\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"score\"\\t\"$12\"\\t\"$6}}' > "+
    			outPrefix+".${x}.bedpe";
    	lf.writeFile(file, line, true);
    	line = "    paste "+outPrefix+".${x}.R3.uniq.bed "+outPrefix+".${x}.R4.uniq.bed |"+
    	"awk '{score=$5;if(score>$11){score=$11};if(($1==$7 && $2<$8) || ($1<$7))"+
    			"{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$4\"\\t\"score\"\\t\"$6\"\\t\"$12}"+
    			"else{print $7\"\\t\"$8\"\\t\"$9\"\\t\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"score\"\\t\"$12\"\\t\"$6}}' >> "+
    			outPrefix+".${x}.bedpe";
    	lf.writeFile(file, line, true);
    	line = "    paste "+outPrefix+".${x}.R5.uniq.bed "+outPrefix+".${x}.R6.uniq.bed |"+
    	"awk '{score=$5;if(score>$11){score=$11};if(($1==$7 && $2<$8) || ($1<$7))"+
    			"{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$4\"\\t\"score\"\\t\"$6\"\\t\"$12}"+
    			"else{print $7\"\\t\"$8\"\\t\"$9\"\\t\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"score\"\\t\"$12\"\\t\"$6}}' >> "+
    			outPrefix+".${x}.bedpe";
    	lf.writeFile(file, line, true);
    	line = "    paste "+outPrefix+".${x}.R7.uniq.bed "+outPrefix+".${x}.R8.uniq.bed |"+
    	"awk '{score=$5;if(score>$11){score=$11};if(($1==$7 && $2<$8) || ($1<$7))"+
    			"{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$4\"\\t\"score\"\\t\"$6\"\\t\"$12}"+
    			"else{print $7\"\\t\"$8\"\\t\"$9\"\\t\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"score\"\\t\"$12\"\\t\"$6}}' >> "+
    			outPrefix+".${x}.bedpe";
    	lf.writeFile(file, line, true);
    	line = "done";
    	lf.writeFile(file, line, true);
    	line = "cat "+outPrefix+".1_2.bedpe "+outPrefix+".2_1.bedpe > "+outPrefix+".bedpe";
    	lf.writeFile(file, line, true);
    	line = "cut -f8    < "+outPrefix+".bedpe"+" | LANG=C sort -n | uniq -c > "+outPrefix+".bedpe.qc.dist.txt";
    	lf.writeFile(file, line, true);
    	line = "cut -f9,10 < "+outPrefix+".bedpe"+" | LANG=C sort -n | uniq -c > "+outPrefix+".bedpe.strand.dist.txt";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+"*sai";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+"*clean.sam";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+"*uniq*sam";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+"*uniq.bed";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+".1_2.bedpe";
    	lf.writeFile(file, line, true);
    	line = "rm "+outPrefix+".2_1.bedpe";
    	lf.writeFile(file, line, true);
    }
}
