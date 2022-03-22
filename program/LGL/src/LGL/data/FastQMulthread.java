package LGL.data;


/**
 * @author sun
 */
import multhread.BigFileProcess;

public class FastQMulthread {

    private String[] fastq = new String[BigFileProcess.SIZE];
    
    
    public void setFastq(int i, String str) {
    	//System.out.println("-_++--AAAA\n");
    	fastq[i] = str;
    }
    
    public String[] getFastq() {
    	//System.out.println("---AAAA\n");
        return fastq;
    }
}