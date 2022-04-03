package multhread;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.concurrent.BlockingQueue;
//gz
import java.io.BufferedOutputStream;
//import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
//import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
//if gz
import java.io.RandomAccessFile;
//time
import java.util.Calendar;

import LGL.data.FastQMulthread;

/**

* Checks if a file is gzipped.

*

* @param f

* @return

*/

public class LinkerFilteringRead implements Runnable {
	
	private File fastq1;
	private File fastq2;
	private BlockingQueue<Object> queue;
	private Object END;
	
	public LinkerFilteringRead(File fastq1, File fastq2, BlockingQueue<Object> queue, Object END) {
		this.fastq1 = fastq1;
		this.fastq2 = fastq2;
		this.queue = queue;
		this.END = END;
	}
	
	public static boolean isGZipped(File f) {
		int magic = 0;
	
		try {
			RandomAccessFile raf = new RandomAccessFile(f, "r");
			magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
			raf.close();
		} catch (Throwable e) {
			e.printStackTrace(System.err);
		}
		
		return magic == GZIPInputStream.GZIP_MAGIC;
	}
    
	public void run() {
		try {
			BufferedReader reader1;
			BufferedReader reader2;
			Calendar rightNow = Calendar.getInstance();
			if(isGZipped(fastq1)) {
				System.out.println("[" + rightNow.getTime().toString() +"] Step1: Linker filtering with gzip fastq file! ...");
				reader1  = new BufferedReader(
		                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq1))));
				reader2  = new BufferedReader(
		                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq2))));
			}else {
				System.out.println("[" + rightNow.getTime().toString() +"] Step1: Linker filtering with fastq file! ...");
				reader1 = new BufferedReader(new FileReader(fastq1));
				reader2 = new BufferedReader(new FileReader(fastq2));
			}
			//BufferedReader reader1 = new BufferedReader(new FileReader(fastq1));
			//BufferedReader reader2 = new BufferedReader(new FileReader(fastq2));
			String line1 = "";
			String line2 = "";
			ArrayList<String> line = new ArrayList<String>();
			ArrayList<FastQMulthread> fastqList = new ArrayList<FastQMulthread>();
			int lineNum = 0;
			int mark = 0;
			while ((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null) {
				line.add(line1);
				line.add(line2);
				lineNum++;
				mark = 1;
				if (lineNum % BigFileProcess.SIZE == 0) {
					FastQMulthread fq1 = new FastQMulthread();
					FastQMulthread fq2 = new FastQMulthread();
					for (int i = 0; i < line.size(); i += 2) {
						fq1.setFastq(i / 2, line.get(i));
						fq2.setFastq(i / 2, line.get(i + 1));
					}
					line.clear();
					fastqList.add(fq1);
					fastqList.add(fq2);
				}
				if (lineNum % 2000 == 0) {
					queue.put(fastqList);
					fastqList = new ArrayList<FastQMulthread>();
					mark = 2;
				}
			}
			if (mark == 1) {// put <2000 into queue
				queue.put(fastqList);
			}
			queue.put(END);
			reader1.close();
			reader2.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}