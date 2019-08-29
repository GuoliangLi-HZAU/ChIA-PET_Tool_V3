package multhread;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import LGL.align.LocalAlignment;
import LGL.chiapet.LinkerFiltering_FastQ_PET_longread;
import LGL.data.FastQMulthread;
import LGL.util.SeqUtil;

public class LinkerFilteringProcess_Long implements Runnable {
	
	private BlockingQueue<Object> queue;
	private Object END;
	private LinkerFiltering_FastQ_PET_longread lfp;
	private BufferedWriter[] arrayOfBufferedWriter;
	
	public LinkerFilteringProcess_Long(BlockingQueue<Object> queue, Object END, LinkerFiltering_FastQ_PET_longread lfp, 
			BufferedWriter[] arrayOfBufferedWriter) {
		this.queue = queue;
		this.END = END;
		this.lfp = lfp;
		this.arrayOfBufferedWriter = arrayOfBufferedWriter;
	}
	
	public void run() {
		try {
			while (true) {
				Object obj = queue.take();
				if (obj == END) {
					queue.put(obj);
					break;
				} else {
					@SuppressWarnings("unchecked")
					ArrayList<FastQMulthread> fastqList = (ArrayList<FastQMulthread>)obj;
					process(fastqList);
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void process(ArrayList<FastQMulthread> fastqList) throws IOException {
		for (int i = 0; i < fastqList.size(); i += 2) {
			processOnePET(fastqList.get(i), fastqList.get(i + 1));
			synchronized(lfp) {
				lfp.setnPETs(lfp.getnPETs() + 1);
			}
		}
	}

	public void processOnePET(FastQMulthread fastQ1, FastQMulthread fastQ2) throws IOException {
	    int[] arrayOfInt1 = processOneSequence(fastQ1.getFastq()[1], 0);
	    int[] arrayOfInt2 = processOneSequence(fastQ2.getFastq()[1], 1);
	    int i = arrayOfInt1[1] - arrayOfInt1[0];
	    int j = arrayOfInt2[1] - arrayOfInt2[0];
	    int index = 0;
	    if ((arrayOfInt1[3] >= lfp.getminimum_linker_alignment_score()) 
	    		&& (arrayOfInt2[3] >= lfp.getminimum_linker_alignment_score()) 
	    		&& (i >= lfp.getminimum_tag_length()) 
	    		&& (j >= lfp.getminimum_tag_length()) 
	    		&& (i <= lfp.getmaximum_tag_length()) 
	    		&& (j <= lfp.getmaximum_tag_length())) {
	    	index = (arrayOfInt1[2] * lfp.getnLinkers() + arrayOfInt2[2]) * 8;
	    	synchronized(lfp) {
        		lfp.setlinkerCompositionDistribution(arrayOfInt1[2], arrayOfInt2[2], lfp.getlinkerCompositionDistribution(arrayOfInt1[2], arrayOfInt2[2]) +
        				1);
        	}
	    	synchronized(arrayOfBufferedWriter) {
		    	if ((i <= 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 1);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 1);
		    		}
		    	} else if ((i <= 55) && (j > 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 2);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 2);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 3);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 3);
		    		}
		    	} else if ((i > 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 4);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 4);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 5);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 5);
		    		}
		    	} else {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 6);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 6);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 7);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 7);
		    		}
			    }
	    	}
	    } else if ((arrayOfInt1[3] >= lfp.getminimum_linker_alignment_score()) && (i >= lfp.getminimum_tag_length()) && 
	    		(i <= lfp.getmaximum_tag_length())) {
	    	if (arrayOfInt1[2] == 0) {
	    		index = (arrayOfInt1[2] * lfp.getnLinkers() + 1) * 8;
	    	} else {
	    		index = (arrayOfInt1[2] * lfp.getnLinkers()) * 8;
	    	}
	    	synchronized(lfp) {
        		lfp.setlinkerCompositionDistribution(arrayOfInt1[2], lfp.getnLinkers(), lfp.getlinkerCompositionDistribution(arrayOfInt1[2], 
        				lfp.getnLinkers()) + 1);
        	}
	    	synchronized(arrayOfBufferedWriter) {
		    	if ((i <= 55) && (fastQ2.getFastq()[1].length() <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write3(fastQ2, index + 1);
		    		} else {
		    			write4(fastQ2, index + 1);
		    		}
		    	} else if ((i <= 55) && (fastQ2.getFastq()[1].length() > 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 2);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 2);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write3(fastQ2, index + 3);
		    		} else {
		    			write4(fastQ2, index + 3);
		    		}
		    	} else if ((i > 55) && (fastQ2.getFastq()[1].length() <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 4);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 4);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write3(fastQ2, index + 5);
		    		} else {
		    			write4(fastQ2, index + 5);
		    		}
		     	} else {
		     		if (lfp.getflip_head_tag() == 1) {
		     			write1(fastQ1, arrayOfInt1, index + 6);
		    		} else {
		    			write2(fastQ1, arrayOfInt1, index + 6);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write3(fastQ2, index + 7);
		    		} else {
		    			write4(fastQ2, index + 7);
		    		}
		     	}
	    	}
	    } else if ((arrayOfInt2[3] >= lfp.getminimum_linker_alignment_score()) && (j >= lfp.getminimum_tag_length()) && 
	    		(j <= lfp.getmaximum_tag_length())) {
	    	if (arrayOfInt2[2] == 0) {
	    		index = (lfp.getnLinkers() + arrayOfInt2[2]) * 8;
	    	} else {
	    		index = (arrayOfInt2[2]) * 8;
	    	}
	    	synchronized(lfp) {
        		lfp.setlinkerCompositionDistribution(lfp.getnLinkers(), arrayOfInt2[2], lfp.getlinkerCompositionDistribution(lfp.getnLinkers(), 
        				arrayOfInt2[2]) + 1);
        	}
	    	synchronized(arrayOfBufferedWriter) {
		    	if ((fastQ1.getFastq()[1].length() <= 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write3(fastQ1, index);
		    		} else {
		    			write4(fastQ1, index);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 1);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 1);
		    		}
		    	} else if ((fastQ1.getFastq()[1].length() <= 55) && (j > 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write3(fastQ1, index + 2);
		    		} else {
		    			write4(fastQ1, index + 2);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 3);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 3);
		    		}
		    	} else if ((fastQ1.getFastq()[1].length() > 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write3(fastQ1, index + 4);
		    		} else {
		    			write4(fastQ1, index + 4);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 5);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 5);
		    		}
		    	} else {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write3(fastQ1, index + 6);
		    		} else {
		    			write4(fastQ1, index + 6);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 7);
		    		} else {
		    			write2(fastQ2, arrayOfInt2, index + 7);
		    		}
		    	}
	    	}
	    } else {
	    	synchronized(lfp) {
        		lfp.setnAmbiguousLinkerComposition(lfp.getnAmbiguousLinkerComposition() + 1);
        	}
	    	if (lfp.getoutput_data_with_ambiguous_linker_info() == 1) {
	    		index = lfp.getnLinkers() * lfp.getnLinkers() * 8;
	    		synchronized(arrayOfBufferedWriter) {
	    			write4(fastQ1, index);
	    			write4(fastQ2, index + 1);
	    		}
	    	}
	    }
	    if (lfp.getdebug_level() >= 2) {
	    	synchronized(lfp) {
		    	lfp.getdebug_output().write(fastQ1.getFastq()[1].substring(arrayOfInt1[0], arrayOfInt1[1]) + "\t" + 
	    	fastQ1.getFastq()[1].substring(arrayOfInt1[1]) + "\t" + arrayOfInt1[2] + "\t" + arrayOfInt1[3] + "\t" + arrayOfInt1[4] + "\t" + 
		    			arrayOfInt1[5] + "\t" + arrayOfInt1[6] + "\t");
		    	lfp.getdebug_output().write(fastQ2.getFastq()[1].substring(arrayOfInt2[0], arrayOfInt2[1]) + "\t" + 
	    	fastQ2.getFastq()[1].substring(arrayOfInt2[1]) + "\t" + arrayOfInt2[2] + "\t" + arrayOfInt2[3] + "\t" + arrayOfInt2[4] + "\t" + 
		    			arrayOfInt2[5] + "\t" + arrayOfInt2[6] + "\t");
		    	lfp.getdebug_output().write(fastQ1.getFastq()[0] + "\t");
		    	lfp.getdebug_output().write(fastQ2.getFastq()[0]);
		    	lfp.getdebug_output().newLine();
	    	}
	    }
    }
	
	public int[] processOneSequence(String seq, int iRead) throws IOException {
		
		LocalAlignment localAligner = new LocalAlignment(lfp.getmaxLinkerLength(), lfp.getmaxLinkerLength());
		
	    int bestScore = -1;
	    int secondBestScore = -1;
	    int bestLinkerIndex = -1;
	    int minJ = -1;
	    int minI = -1;
	    for (int i = 0; i < lfp.getlinkers().length; i++) {
	    	localAligner.align(lfp.getlinkers()[i], seq);
	    	int score = localAligner.getMaxScore();
	    	if (bestScore < score) {
	    		secondBestScore = bestScore;
	    		bestScore = score;
		        bestLinkerIndex = i;
		        minJ = localAligner.getMinJ();
		        minI = localAligner.getMinI();	        
	    	} else if (secondBestScore < score) {
	    		secondBestScore = score;
	    	}
	    }
	    if ((bestScore >= 0) && (bestScore <= lfp.getmaxLinkerLength())) {
	    	synchronized(lfp) {
				lfp.setscoreDistribution(iRead, bestScore, lfp.getscoreDistribution(iRead, bestScore) + 1);
			}
	    }
	    int secondBestScoreDiff = bestScore - secondBestScore;
	    if ((secondBestScoreDiff >= 0) && (secondBestScoreDiff <= 2 * lfp.getmaxLinkerLength())) {
	    	synchronized(lfp) {
        		lfp.setsecondBestScoreDiffDistribution(iRead, secondBestScoreDiff, lfp.getsecondBestScoreDiffDistribution(iRead, secondBestScoreDiff) + 1);
			}
	    }
	    synchronized(lfp) {
	        if (lfp.getmaxSecondBestScoreDiff() < secondBestScoreDiff) {
        		lfp.setmaxSecondBestScoreDiff(secondBestScoreDiff);
			}
	    }
	    int tag_Start = 0;
	    int tag_End = minJ - minI;
	    if (tag_End < 0) {
	    	tag_End = 0;
	    }
	    if (tag_End < 0) {
	    	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, 0, lfp.gettagLengthDistribution(iRead, 0) + 1);
			}
	    } else if (tag_End >= lfp.getmaximum_tag_length()) {
	    	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, lfp.getmaximum_tag_length() - 1, lfp.gettagLengthDistribution(iRead, lfp.getmaximum_tag_length() - 1) +
        				1);
			}
	    } else {
	    	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, tag_End, lfp.gettagLengthDistribution(iRead, tag_End) + 1);
			}	    	
	    }
	    synchronized(lfp) {
	    	if (lfp.getmaxRealTagLength() < tag_End) {
        		lfp.setmaxRealTagLength(tag_End);
			}
	    }
	    int[] arrayOfInt = new int[8];
	    arrayOfInt[0] = tag_Start;
	    arrayOfInt[1] = tag_End;
	    arrayOfInt[2] = bestLinkerIndex;
	    arrayOfInt[3] = bestScore;
	    arrayOfInt[4] = secondBestScoreDiff;
	    arrayOfInt[5] = minJ;
	    arrayOfInt[6] = minI;
	    return arrayOfInt;
    }
	
	public void write1(FastQMulthread fastQ, int[] arrayOfInt, int index) {
		try {
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[0]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(SeqUtil.revComplement(fastQ.getFastq()[1].substring(arrayOfInt[0], arrayOfInt[1])));
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[2]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(new StringBuffer(fastQ.getFastq()[3].substring(arrayOfInt[0], arrayOfInt[1])).reverse().toString());
			arrayOfBufferedWriter[index].newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void write2(FastQMulthread fastQ, int[] arrayOfInt, int index) {
		try {
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[0]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[1].substring(arrayOfInt[0], arrayOfInt[1]));
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[2]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[3].substring(arrayOfInt[0], arrayOfInt[1]));
			arrayOfBufferedWriter[index].newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void write3(FastQMulthread fastQ, int index) {
		try {
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[0]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(SeqUtil.revComplement(fastQ.getFastq()[1]));
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[2]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(new StringBuffer(fastQ.getFastq()[3]).reverse().toString());
			arrayOfBufferedWriter[index].newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void write4(FastQMulthread fastQ, int index) {
		try {
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[0]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[1]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[2]);
			arrayOfBufferedWriter[index].newLine();
			arrayOfBufferedWriter[index].write(fastQ.getFastq()[3]);
			arrayOfBufferedWriter[index].newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
