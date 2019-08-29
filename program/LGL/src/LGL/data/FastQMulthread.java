package LGL.data;

/**
 * @author sun
 */
import multhread.BigFileProcess;

public class FastQMulthread {

    private String[] fastq = new String[BigFileProcess.SIZE];
    
    public void setFastq(int i, String str) {
    	fastq[i] = str;
    }
    
    public String[] getFastq() {
        return fastq;
    }
}