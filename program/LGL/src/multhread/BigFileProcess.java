package multhread;

import java.io.BufferedWriter;
import java.io.File;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import LGL.chiapet.LinkerFiltering_FastQ_PET;
import LGL.chiapet.LinkerFiltering_FastQ_PET_longread;
import LGL.util.MappingStatistics;

public class BigFileProcess {
	
	private Thread linkerFilteingRead;
	private Thread[] linkerFilteingProccess;
	private BlockingQueue<Object> queue;
	public final Object END = new Object();
	public static final int SIZE = 4;
	private Thread mappingRead;
	private Thread[] mappingProcess;
		
	public BigFileProcess(File fastq1, File fastq2, LinkerFiltering_FastQ_PET lfp, BufferedWriter[] fileOut, int threadNum) {
		queue = new LinkedBlockingQueue<Object>(threadNum);
		linkerFilteingRead = new Thread(new LinkerFilteringRead(fastq1, fastq2, queue, END));
		linkerFilteingProccess = new Thread[threadNum];
		for (int i = 0; i < linkerFilteingProccess.length; i++) {
			linkerFilteingProccess[i] = new Thread(new LinkerFilteringProcess(queue, END, lfp, fileOut));
		}
	}
	
	public BigFileProcess(File fastq1, File fastq2, LinkerFiltering_FastQ_PET_longread lfp, BufferedWriter[] arrayOfBufferedWriter, int threadNum) {
		queue = new LinkedBlockingQueue<Object>(threadNum);
		linkerFilteingRead = new Thread(new LinkerFilteringRead(fastq1, fastq2, queue, END));
		linkerFilteingProccess = new Thread[threadNum];
		for (int i = 0; i < linkerFilteingProccess.length; i++) {
			linkerFilteingProccess[i] = new Thread(new LinkerFilteringProcess_Long(queue, END, lfp, arrayOfBufferedWriter));
		}
	}
	
	public BigFileProcess(String samFile1, String samFile2, int mappingScoreCutoff, int threadNum, BufferedWriter outBedpe, MappingStatistics mappingStatistics) {
		queue = new LinkedBlockingQueue<Object>(threadNum);
		mappingRead = new Thread(new MappingRead(samFile1, samFile2, queue, END));
		mappingProcess = new Thread[threadNum];
		for (int i = 0; i < mappingProcess.length; i++) {
			mappingProcess[i] = new Thread(new MappingProcess(mappingScoreCutoff, outBedpe, mappingStatistics, queue, END));
		}
	}
	
	public void start() {
        linkerFilteingRead.start();
        for (int i = 0; i < linkerFilteingProccess.length; i++) {
            linkerFilteingProccess[i].start();
        }
    }
	
	public void join() {
		try {
            linkerFilteingRead.join();
            for (int i = 0; i < linkerFilteingProccess.length; i++) {
                linkerFilteingProccess[i].join();
            }
        } catch (InterruptedException e) {
        	e.printStackTrace();
        }
    }
	
	public void start_mapping() {
		mappingRead.start();
        for (int i = 0; i < mappingProcess.length; i++) {
        	mappingProcess[i].start();
        }
    }
	
	public void join_mapping() {
		try {
			mappingRead.join();
            for (int i = 0; i < mappingProcess.length; i++) {
            	mappingProcess[i].join();
            }
        } catch (InterruptedException e) {
        	e.printStackTrace();
        }
    }
}