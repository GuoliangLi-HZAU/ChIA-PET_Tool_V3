/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.io.BufferedReader;
import java.io.IOException;

/**
 *
 * @author ligl
 */
public class FastQ {

    private String[] fastq = null;//fastq字符串数组

    public FastQ() {

    }

    public static FastQ load(BufferedReader fastqFile) throws IOException {//读取fastq
        String[] fastq_temp = new String[4];//暂存读出字符串
        String line;
        int nLines = 0;//读出字符串行数
        if ((line = fastqFile.readLine()) != null) {
            fastq_temp[0] = new String(line);
            nLines++;
        }
        if ((line = fastqFile.readLine()) != null) {
            fastq_temp[1] = new String(line);
            nLines++;
        }
        if ((line = fastqFile.readLine()) != null) {
            fastq_temp[2] = new String(line);
            nLines++;
        }
        if ((line = fastqFile.readLine()) != null) {
            fastq_temp[3] = new String(line);
            nLines++;
        }

        FastQ fastQ_1 = new FastQ();
        if (nLines == 4) {
            fastQ_1.setFastq(fastq_temp);
        }
        return fastQ_1;
    }

    /**
     * @return the fastq
     */
    public String[] getFastq() {//获取fastq字符串数组
        return fastq;
    }

    /**
     * @param fastq the fastq to set
     */
    public void setFastq(String[] fastq) {//设置fastq字符串数组
        this.fastq = fastq;
    }
}
