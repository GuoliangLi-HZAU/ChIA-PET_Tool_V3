package multhread;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import LGL.data.FastQMulthread;

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
    
	public void run() {
		try {
			BufferedReader reader1 = new BufferedReader(new FileReader(fastq1));
			BufferedReader reader2 = new BufferedReader(new FileReader(fastq2));
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