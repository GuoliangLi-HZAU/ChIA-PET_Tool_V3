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

	public int alignsegpart(String seq1, String seq2, int start1, int end1, int start2, int end2) throws IOException {
		int newscore = -1;
		if(end1-start1<13 || end2-start2<13) {
			return newscore;
		}
		String three1 = seq1.substring(start1, end1);
		String five2 = SeqUtil.revComplement(seq2.substring(start2, end2));
		newscore = getlocalAlignmaxScore(three1, five2);
		//System.out.println(three1+" "+five2+" "+newscore);
		return newscore;
	}
	public void processOnePET(FastQMulthread fastQ1, FastQMulthread fastQ2) throws IOException {
	    int[] arrayOfInt1 = processOneSequence(fastQ1.getFastq()[1], 0);
	    int[] arrayOfInt2 = processOneSequence(fastQ2.getFastq()[1], 1);
	    int i = arrayOfInt1[1] - arrayOfInt1[0];
	    int j = arrayOfInt2[1] - arrayOfInt2[0];
	    
	    int index = 0;
	    // read1 and read2 align score > set minimum score 14, then process more than 1 linker case
	    int Noverlap1 = 1;
    	int Noverlap2 = 1;
    	int newscore = -1;

	    if((arrayOfInt1[3] >= lfp.getminimum_linker_alignment_score()) 
	    		&& (arrayOfInt2[3] >= lfp.getminimum_linker_alignment_score()) &&
	    		lfp.search_all_linker ) {
	    	//System.out.println("search_all_linker mode!!");
	    	// more than 1 linker in one part-read
	    	
	    	//System.out.println("1111AAA " + arrayOfInt1[7] + " " + arrayOfInt1[8]);
	    	/*
	    	 * Read1:
	    	 * [maxJ, secondminJ-secondminI]    [secondmaxJ, end] Or [maxJ, end]
	    	 * (7, 9-10) (11, end]  Or (7, end]
	    	 * 
	    	 * Read2:
	    	 * [0, minJ-minI] [maxJ, secondminJ] [secondmaxJ, end] Or [0, minJ-minI] [maxJ, end]
	    	 * [0, 5-6) (7,9) (11,end]  Or [0, 5-6) (7,end]
	    	 * 
	    	 * */
	    	int[] index1 = new int[4];
	    	int[] index2 = new int[6];
	    	if(arrayOfInt1[4] == 0) {
	    		index1[0] = arrayOfInt1[7]+1;
	    		index1[1] = arrayOfInt1[9]-arrayOfInt1[10];
	    		index1[2] = arrayOfInt1[11]+1;
	    		index1[3] = fastQ1.getFastq()[1].length();
	    	}else {
	    		index1[0] = arrayOfInt1[7]+1;
	    		index1[1] = fastQ1.getFastq()[1].length();
	    		index1[2] = 0;
	    		index1[3] = 1;
	    	}
	    	if(arrayOfInt2[4] == 0) {
	    		index2[0] = 0;
	    		index2[1] = arrayOfInt2[5]-arrayOfInt2[6];
	    		index2[2] = arrayOfInt2[7]+1;
	    		index2[3] = arrayOfInt2[9];
	    		index2[4] = arrayOfInt2[11]+1;
	    		index2[5] = fastQ2.getFastq()[1].length();
	    	}else {
	    		index2[0] = 0;
	    		index2[1] = arrayOfInt2[5]-arrayOfInt2[6];
	    		index2[2] = arrayOfInt2[7]+1;
	    		index2[3] = fastQ2.getFastq()[1].length();
	    		index2[4] = 0;
	    		index2[5] = 1;
	    	}
	    	
	    	for(int nl = 0; nl<4; nl+=2) {
	    		for(int nk = 0; nk < 6; nk+=2) {
	    			int minlen = (index1[nl+1]-index1[nl]<index2[nk+1]-index2[nk])?(index1[nl+1]-index1[nl]):(index2[nk+1]-index2[nk]);
	    			//System.out.println("sss " + index1[nl] + " " + index1[nl+1] + " "+ index2[nk] + " " + index2[nk+1]);
	    			newscore = alignsegpart(fastQ1.getFastq()[1], fastQ2.getFastq()[1], index1[nl], index1[nl+1], index2[nk], index2[nk+1]);
	        	    if( (float)newscore/minlen > 0.8 ) {
	        	    	Noverlap1 = nl/2+1;
	    		    	Noverlap2 = nk/2;
	    		    	//System.out.println("overlapscore " + Noverlap1 + " " + Noverlap2 + " " + newscore);
	    		    	break;
	        	    }
	    		}
	    	}
    		/* Noverlap1 5'->3'
    		 * ----linker1-----linker2----
    		 * Noverlap2 5'->3'
    		 * ----linker1-----linker2----
    		 * So, (Noverlap1+Noverlap2)%2==1 is normal case
    		 * */
	    }
	    
	    /*System.out.println(fastQ1.getFastq()[1] + " [] " + fastQ2.getFastq()[1] + " " + Nlinker1 + " " +
			       Nlinker2 + " " + chuanleORdajie);
	    
			    for(int nx = 0; nx < 14; nx++) {
			    	System.out.print(arrayOfInt1[nx] + " ");
			    }
			    System.out.println();
			    for(int nx = 0; nx < 14; nx++) {
			    	System.out.print(arrayOfInt2[nx] + " ");
			    }
			    System.out.println();
			  */
	    //System.out.println("---score--- " + arrayOfInt1[3] + " " + arrayOfInt2[3] + " " + i + " " + j);
	    
	    //System.out.println(lfp.getminimum_linker_alignment_score() + " " + lfp.getminimum_tag_length() + " " + lfp.getmaximum_tag_length());
	    if ((arrayOfInt1[3] >= lfp.getminimum_linker_alignment_score()) 
	    		&& (arrayOfInt2[3] >= lfp.getminimum_linker_alignment_score()) 
	    		&& (i >= lfp.getminimum_tag_length()) 
	    		&& (j >= lfp.getminimum_tag_length()) 
	    		&& (i <= lfp.getmaximum_tag_length()) 
	    		&& (j <= lfp.getmaximum_tag_length())
	    		//&& i>15 && j>15 // make sure the length of read is longer than 15
	    		&& Noverlap1 > lfp.printallreads() //only print overlap part, we can know +-
	    		)
	    {
		    //index, array1 array2, 0 0 0 1_1, 0 1 8 1_2, 1 0 16 2_1, 1 1 24 2_2
	    	index = (arrayOfInt1[2] * lfp.getnLinkers() + arrayOfInt2[2]) * 8;
	    	//System.out.println(arrayOfInt1[2] + " " + lfp.getnLinkers() + " " + arrayOfInt2[2] + " " + index);
	    	if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==1) {
	    		if(index<8 || index >= 24) {
	    			index = (arrayOfInt1[2] + 1) * 8;
	    		}
	    	}
	    	synchronized(lfp) {
        		lfp.setlinkerCompositionDistribution(arrayOfInt1[2], arrayOfInt2[2], lfp.getlinkerCompositionDistribution(arrayOfInt1[2], arrayOfInt2[2]) +
        				1);
        	}
	    	synchronized(arrayOfBufferedWriter) {
		    	if ((i <= 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ1, arrayOfInt1, index);
		    			}else
		    			    write2(fastQ1, arrayOfInt1, index);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 1);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ2, arrayOfInt2, index+1);
		    			}else
		    			    write2(fastQ2, arrayOfInt2, index + 1);
		    		}
		    	} else if ((i <= 55) && (j > 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 2);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ1, arrayOfInt1, index+2);
		    			}else
		    			    write2(fastQ1, arrayOfInt1, index + 2);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 3);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ2, arrayOfInt2, index+3);
		    			}else
		    			    write2(fastQ2, arrayOfInt2, index + 3);
		    		}
		    	} else if ((i > 55) && (j <= 55)) {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 4);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ1, arrayOfInt1, index+4);
		    			}else
		    			    write2(fastQ1, arrayOfInt1, index + 4);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 5);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ2, arrayOfInt2, index+5);
		    			}else
		    			    write2(fastQ2, arrayOfInt2, index + 5);
		    		}
		    	} else {
		    		if (lfp.getflip_head_tag() == 1) {
		    			write1(fastQ1, arrayOfInt1, index + 6);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ1, arrayOfInt1, index+6);
		    			}else
		    			    write2(fastQ1, arrayOfInt1, index + 6);
		    		}
		    		if (lfp.getflip_tail_tag() == 1) {
		    			write1(fastQ2, arrayOfInt2, index + 7);
		    		} else {
		    			if(lfp.search_all_linker && Noverlap1>0 && (Noverlap1+Noverlap2)%2==0) {
		    				write1(fastQ2, arrayOfInt2, index+7);
		    			}else
		    			    write2(fastQ2, arrayOfInt2, index + 7);
		    		}
			    }
	    	}
	    } else if ((arrayOfInt1[3] >= lfp.getminimum_linker_alignment_score()) && (i >= lfp.getminimum_tag_length()) && 
	    		(i <= lfp.getmaximum_tag_length())
	    		//&& i>15
	    		) {
	    	// force turned to 1_2 2_1, only cause long reads mode, here is not good??
	    	if(lfp.MAP2Linker.equalsIgnoreCase("true")) {
	    		if (arrayOfInt1[2] == 0) {
		    		index = (arrayOfInt1[2] * lfp.getnLinkers()) * 8;
		    	} else {
		    		index = (arrayOfInt1[2] * lfp.getnLinkers() + 1) * 8;
		    	}
	    	}else {
		    	if (arrayOfInt1[2] == 0) {
		    		index = (arrayOfInt1[2] * lfp.getnLinkers() + 1) * 8;
		    	} else {
		    		index = (arrayOfInt1[2] * lfp.getnLinkers()) * 8;
		    	}
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
	    		(j <= lfp.getmaximum_tag_length())
	    		//&& j>15
	    		) {
	    	// force turned to 1_2 2_1, only cause long reads mode, here is not good??
	    	if(lfp.MAP2Linker.equalsIgnoreCase("true")) {
	    		if (arrayOfInt2[2] == 0) {
	    			index = (arrayOfInt2[2]) * 8;
		    	} else {
		    		index = (lfp.getnLinkers() + arrayOfInt2[2]) * 8;
		    	}
	    	}else {
	    		if (arrayOfInt2[2] == 0) {
		    		index = (lfp.getnLinkers() + arrayOfInt2[2]) * 8;
		    	} else {
		    		index = (arrayOfInt2[2]) * 8;
		    	}
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
	
	// get max local alignment score of seq1 and seq2
	public int getlocalAlignmaxScore(String seq1, String seq2) throws IOException {
		
		LocalAlignment localAligner = new LocalAlignment(seq1.length(), seq2.length());
		
	    localAligner.align(seq1, seq2, 0);
	    int score = localAligner.getMaxScore();
	    
	    return score;
	}
	
	public int[] processOneSequence(String seq, int iRead) throws IOException {
		
		LocalAlignment localAligner = new LocalAlignment(lfp.getmaxLinkerLength(), lfp.getmaxLinkerLength());
		
	    int bestScore = -1;
	    int secondBestScore = -1;
	    int bestLinkerIndex = -1;
	    int minJ = -1;
	    int minI = -1;
	    int maxJ = -1;
	    int maxI = -1;
	    int secondminJ = -1;
	    int secondminI = -1;
	    int secondmaxJ = -1;
	    int secondmaxI = -1;
	    for (int i = 0; i < lfp.getlinkers().length; i++) {
	    	localAligner.align(lfp.getlinkers()[i], seq, lfp.getminimum_tag_length()-1); // 17
	    	int score = localAligner.getMaxScore();
	    	if (bestScore < score || (bestScore == score && maxJ > localAligner.getMaxJ() ) ) {
	    		secondBestScore = bestScore;
	    		secondminJ = minJ;
	    		secondminI = minI;
	    		secondmaxJ = maxJ;
	    		secondmaxI = maxI;
	    		bestScore = score;
		        bestLinkerIndex = i;
		        minJ = localAligner.getMinJ();
		        minI = localAligner.getMinI();
		        maxJ = localAligner.getMaxJ();
		        maxI = localAligner.getMaxI();
	    	}else if(bestScore == score && maxJ < localAligner.getMaxJ()) { 
	    		secondBestScore = bestScore;
	    		secondminJ = localAligner.getMinJ();
	    		secondminI = localAligner.getMinI();
	    		secondmaxJ = localAligner.getMaxJ();
	    		secondmaxI = localAligner.getMaxI();
	    		//bestScore = score;
		        //bestLinkerIndex = i;
	    	}else if (secondBestScore < score) {
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
	    //System.out.println("---chengshengceshiyixia----- " + minI);
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
	    int[] arrayOfInt = new int[14];
	    arrayOfInt[0] = tag_Start;
	    arrayOfInt[1] = tag_End;
	    arrayOfInt[2] = bestLinkerIndex;
	    arrayOfInt[3] = bestScore;
	    arrayOfInt[4] = secondBestScoreDiff;
	    arrayOfInt[5] = minJ;
	    arrayOfInt[6] = minI;
	    arrayOfInt[7] = maxJ;
	    arrayOfInt[8] = maxI;
	    arrayOfInt[9] = secondminJ;
	    arrayOfInt[10] = secondminI;
	    arrayOfInt[11] = secondmaxJ;
	    arrayOfInt[12] = secondmaxI;
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
