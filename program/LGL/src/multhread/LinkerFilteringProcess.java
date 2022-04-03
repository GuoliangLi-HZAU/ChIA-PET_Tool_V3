package multhread;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import LGL.align.LocalAlignment;
import LGL.chiapet.LinkerFiltering_FastQ_PET;
import LGL.data.FastQMulthread;

public class LinkerFilteringProcess implements Runnable {

	private BlockingQueue<Object> queue;
	private Object END;
	private LinkerFiltering_FastQ_PET lfp;
	private BufferedWriter[] fileOut;
	
	public LinkerFilteringProcess(BlockingQueue<Object> queue, Object END, LinkerFiltering_FastQ_PET lfp, BufferedWriter[] fileOut) {
		this.queue = queue;
		this.END = END;
		this.lfp = lfp;
		this.fileOut = fileOut;
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
        int[] results_1 = processOneSequence(fastQ1.getFastq()[1], 0);
        int[] results_2 = processOneSequence(fastQ2.getFastq()[1], 1);
        int tag_length_1 = results_1[1] - results_1[0];
        int tag_length_2 = results_2[1] - results_2[0];
        if ((results_1[3] >= lfp.getminimum_linker_alignment_score())
                && ((results_2[3]) >= lfp.getminimum_linker_alignment_score())
                && (tag_length_1 >= lfp.getminimum_tag_length())
                && (tag_length_2 >= lfp.getminimum_tag_length())
                && (tag_length_1 <= lfp.getmaximum_tag_length())
                && (tag_length_2 <= lfp.getmaximum_tag_length())
                && (results_1[4] >= lfp.getminSecondBestScoreDiff())
                && (results_2[4] >= lfp.getminSecondBestScoreDiff())
                && (results_1[7] == 1)
                && (results_2[7] == 1)) {
        	int index = (results_1[2] * lfp.getnLinkers() + results_2[2]) * 2;
        	synchronized(lfp) {
        		lfp.setlinkerCompositionDistribution(results_1[2], results_2[2], lfp.getlinkerCompositionDistribution(results_1[2], results_2[2]) + 1);
        		// update linker composition distribution
        	}
        	synchronized(fileOut) {
	            if (lfp.getflip_head_tag() == 1) {// output head
	            	fileOut[index].write(fastQ1.getFastq()[0]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(LGL.util.SeqUtil.revComplement(fastQ1.getFastq()[1].substring(results_1[0], results_1[1])));
	                fileOut[index].newLine();
	                fileOut[index].write(fastQ1.getFastq()[2]);
	                fileOut[index].newLine();
	                fileOut[index].write(new StringBuilder(fastQ1.getFastq()[3].substring(results_1[0], results_1[1])).reverse().toString());
	                fileOut[index].newLine();
	            } else {
	            	fileOut[index].write(fastQ1.getFastq()[0]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[1].substring(results_1[0], results_1[1]));
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[2]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[3].substring(results_1[0], results_1[1]));
	            	fileOut[index].newLine();
	            }
	            if (lfp.getflip_tail_tag() == 1) {// output tail
	            	fileOut[index + 1].write(fastQ2.getFastq()[0]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(LGL.util.SeqUtil.revComplement(fastQ2.getFastq()[1].substring(results_2[0], results_2[1])));
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[2]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(new StringBuilder(fastQ2.getFastq()[3].substring(results_2[0], results_2[1])).reverse().toString());
	            	fileOut[index + 1].newLine();
	            } else {
	            	fileOut[index + 1].write(fastQ2.getFastq()[0]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[1].substring(results_2[0], results_2[1]));
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[2]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[3].substring(results_2[0], results_2[1]));
	            	fileOut[index + 1].newLine();
	            }
        	}
        } else {// PETs with ambiguous linkers
        	synchronized(lfp) {
        		lfp.setnAmbiguousLinkerComposition(lfp.getnAmbiguousLinkerComposition() + 1);
        	}
            if (lfp.getoutput_data_with_ambiguous_linker_info() == 1) {
            	int index = lfp.getnLinkers() * lfp.getnLinkers() * 2;
            	synchronized(fileOut) {
	            	fileOut[index].write(fastQ1.getFastq()[0]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[1]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[2]);
	            	fileOut[index].newLine();
	            	fileOut[index].write(fastQ1.getFastq()[3]);
	            	fileOut[index].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[0]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[1]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[2]);
	            	fileOut[index + 1].newLine();
	            	fileOut[index + 1].write(fastQ2.getFastq()[3]);
	            	fileOut[index + 1].newLine();
            	}
            }
        }       
        if (lfp.getdebug_level() >= 2) {// output the sequences and linker alignment for debugging purpose
        	synchronized(lfp) {
	            lfp.getdebug_output().write(fastQ1.getFastq()[1].substring(results_1[0], results_1[1]) + "\t" + fastQ1.getFastq()[1].substring(results_1[1]) +"\t" +
        	results_1[2] + "\t" + results_1[3] + "\t" + results_1[4] + "\t" + results_1[5] + "\t" + results_1[6] + "\t");
	            lfp.getdebug_output().write(fastQ2.getFastq()[1].substring(results_2[0], results_2[1]) + "\t" + fastQ2.getFastq()[1].substring(results_2[1]) +"\t" +
        	results_2[2] + "\t" + results_2[3] + "\t" + results_2[4] + "\t" + results_2[5] + "\t" + results_2[6] + "\t");
	            lfp.getdebug_output().write(fastQ1.getFastq()[0] + "\t");
	            lfp.getdebug_output().write(fastQ2.getFastq()[0]);
	            lfp.getdebug_output().newLine();
        	}
        }
    }
	
	// iRead: 0 for read1, 1 for read2
	public int[] processOneSequence(String seq, int iRead) throws IOException {
		LocalAlignment localAligner = new LocalAlignment(lfp.getmaxLinkerLength(), lfp.getmaxLinkerLength());
		
        // align different linkers
        int bestScore = -1;
        int secondBestScore = -1;
        int bestLinkerIndex = -1;
        int minI = -1;// index in linker
        int minJ = -1;// index in sequence
        int barcodeStatus = -1;
        for (int i = 0; i < lfp.getlinkers().length; i++) {
        	localAligner.align(lfp.getlinkers()[i], seq, lfp.getminimum_tag_length()-1);
            int score = localAligner.getMaxScore();
            if (bestScore < score) {
                secondBestScore = bestScore;
                bestScore = score;
                bestLinkerIndex = i;
                minI = localAligner.getMinI();// index in linker
                minJ = localAligner.getMinJ();// index in sequence
                barcodeStatus = testBarcode(localAligner, lfp.getbarCodeStart(i), lfp.getbarCodeLength(i));
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
        //===vvv 2013-08-20: the original implementation set the maximum output tag length as 20
        //int tag_Start = minJ - minI - 20;
        //if (tag_Start < 0) {
        //    tag_Start = 0;
        //}
        //===^^^ 2013-08-20: set the output tag length as the available reads
        int tag_Start = 0;
        int tag_End = minJ - minI;
        if(lfp.AutoLinker.equalsIgnoreCase("true") && minI>3) { //can clip shoter than 4
	    	tag_End = minJ;
	    }
        if (tag_End < 0) {
            tag_End = 0;
        }
        if (tag_End < 0) {
        	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, 0, lfp.gettagLengthDistribution(iRead, 0) + 1);
			}
        } else if (tag_End >= lfp.getmaximum_tag_length()) {
        	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, lfp.getmaximum_tag_length() - 1, lfp.gettagLengthDistribution(iRead, lfp.getmaximum_tag_length() - 1) + 1);
			}
        } else {
        	synchronized(lfp) {
        		lfp.settagLengthDistribution(iRead, tag_End, lfp.gettagLengthDistribution(iRead, tag_End) + 1);
			}
        }
        // record the real tag length, for output the tag length distribution purpose
        synchronized(lfp) {
            if (lfp.getmaxRealTagLength() < tag_End) {
        		lfp.setmaxRealTagLength(tag_End);
			}
        }
        // sequence index starts with 0
        int[] results = new int[8];
        results[0] = tag_Start;
        results[1] = tag_End;
        results[2] = bestLinkerIndex;
        results[3] = bestScore;
        results[4] = secondBestScoreDiff;
        results[5] = minJ;
        results[6] = minI;
        results[7] = barcodeStatus;
        return results;
    }
	
	// 1: barcode test passed; 0: barcode test failed
	public int testBarcode(LocalAlignment a, int offset, int len) {
        if ((offset < 0) || (len < 0)) {
            // either offset or len is less than 0, barcode test is not applicable
            return 1;
        }
        // ****** assumptions ******
        // 1) the aligned strings are from the beginning of the original sequences, and to the end of the best 
        //local alignment
        // 2) the aligned str1 is the linker sequence
        // 3) the insertions and deletions are represented with '-'
        // 4) matched alignment status is represented with '|'
        String alignedStr1 = a.getAlignedStr1();
        String alignedStatus = a.getAlignedStatus();
        //StringBuilder trimmedAlignedStatus = new StringBuilder();
        StringBuffer trimmedAlignedStatus = new StringBuffer();
        for (int i = 0; i < alignedStr1.length(); i++) {
            if (alignedStr1.charAt(i) != '-') {
                trimmedAlignedStatus.append(alignedStatus.charAt(i));
            }
        }
        if (lfp.getdebug_level() >= 2) {
        	synchronized(lfp) {
	            if (lfp.getiOutput() < lfp.getnOutput()) {
	                System.out.println("originalStr1: " + a.getStr1());
	                System.out.println("alignedStr1:  " + alignedStr1);
	                System.out.println("status:       " + alignedStatus);
	                System.out.println("alignedStr2:  " + a.getAlignedStr2());
	                System.out.println("originalStr2: " + a.getStr2());
	                System.out.println("trimmedAlignedStatus: " + trimmedAlignedStatus);
	                lfp.setiOutput(lfp.getiOutput() + 1);
	            }
        	}
        }
        int barcodePassed = 1;// 1: barcode test passed; 0: barcode test failed
        if (trimmedAlignedStatus.length() < offset + len - 1) {
            // the barcode is not completely covered
            barcodePassed = 0;
        } else {
            int i = offset - 1;
            for (int j = 0; j < len; j++) {
                if (trimmedAlignedStatus.charAt(i) != '|') {
                    barcodePassed = 0;
                    break;
                }
                i++;
            }
        }
        return barcodePassed;
	}
}