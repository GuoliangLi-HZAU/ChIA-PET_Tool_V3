package multhread;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;

import LGL.util.MappingStatistics;

public class MappingProcess implements Runnable {
	
	private int mappingScoreCutoff;
	private BufferedWriter outBedpe;
	private MappingStatistics mappingStatistics;
	private BlockingQueue<Object> queue;
	private Object END;

	public MappingProcess(int mappingScoreCutoff, BufferedWriter outBedpe, MappingStatistics mappingStatistics, BlockingQueue<Object> queue, Object END) {
		this.mappingScoreCutoff = mappingScoreCutoff;
		this.outBedpe = outBedpe;
		this.mappingStatistics = mappingStatistics;
		this.queue = queue;
		this.END = END;
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
					ArrayList<String> list = (ArrayList<String>)obj;
					process(list);
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
	
	public void process(ArrayList<String> list) {
		int[][] statistics = new int[3][3];
		Arrays.fill(statistics[0], 0);
        Arrays.fill(statistics[1], 0);
        Arrays.fill(statistics[2], 0);
		for (int i = 0; i < list.size(); i = i + 2) {
			String[] strs1 = list.get(i).split("\t");
			String[] strs2 = list.get(i + 1).split("\t");
			if (!strs1[0].equals(strs2[0])) {
				System.out.println("Error: read names of PET ends don't match");
				System.out.println(strs1[0]+"<------>"+strs2[0]);
				continue;
			}
			//System.out.println(strs1[0]);
			
			if (classify(strs1) == 0 && classify(strs2) == 0) {
				statistics[0][0]++;
			} else if (classify(strs1) == 1 && classify(strs2) == 0) {
				statistics[1][0]++;
			} else if (classify(strs1) == 2 && classify(strs2) == 0) {
				statistics[2][0]++;
			} else if (classify(strs1) == 0 && classify(strs2) == 1) {
				statistics[0][1]++;
			} else if (classify(strs1) == 1 && classify(strs2) == 1) {
				statistics[1][1]++;
				buildBedpe(strs1, strs2);
			} else if (classify(strs1) == 2 && classify(strs2) == 1) {
				statistics[2][1]++;
			} else if (classify(strs1) == 0 && classify(strs2) == 2) {
				statistics[0][2]++;
			} else if (classify(strs1) == 1 && classify(strs2) == 2) {
				statistics[1][2]++;
			} else {
				statistics[2][2]++;
			}
		}
		try {
			synchronized(outBedpe) {
				outBedpe.flush();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		synchronized(mappingStatistics) {
			mappingStatistics.addStatistics(statistics);
		}
	}
	
	public int classify(String[] strs) {
        int category_temp = 0;
        if (strs[2].equals("*")) {
            category_temp = 0;
        } else if (Integer.valueOf(strs[4]) >= mappingScoreCutoff) {
            category_temp = 1;
        } else {
            category_temp = 2;
        }
        return category_temp;
    }
    
    public void buildBedpe(String[] strs1, String[] strs2) {
    	StringBuilder read1 = new StringBuilder();
    	read1.append(strs1[2]);
    	read1.append("\t");
    	read1.append(Integer.valueOf(strs1[3]) - 1);
    	read1.append("\t");
    	read1.append(Integer.valueOf(strs1[3]) - 1 + strs1[9].length());
    	read1.append("\t");
    	StringBuilder read2 = new StringBuilder();
    	read2.append(strs2[2]);
    	read2.append("\t");
    	read2.append(Integer.valueOf(strs2[3]) - 1);
    	read2.append("\t");
    	read2.append(Integer.valueOf(strs2[3]) - 1 + strs2[9].length());
    	read2.append("\t");
    	StringBuilder bedpe1 = new StringBuilder();
    	bedpe1.append(read2);
    	bedpe1.append(read1);
    	bedpe1.append(strs1[0]);
    	bedpe1.append("\t");
    	bedpe1.append(strs2[4]);
    	bedpe1.append("\t");
    	bedpe1.append(getStrand(strs2[1]));
    	bedpe1.append("\t");
    	bedpe1.append(getStrand(strs1[1]));
    	StringBuilder bedpe2 = new StringBuilder();
    	bedpe2.append(read1);
    	bedpe2.append(read2);
    	bedpe2.append(strs1[0]);
    	bedpe2.append("\t");
    	bedpe2.append(strs1[4]);
    	bedpe2.append("\t");
    	bedpe2.append(getStrand(strs1[1]));
    	bedpe2.append("\t");
    	bedpe2.append(getStrand(strs2[1]));
		try {
			synchronized(outBedpe) {
				if (strs1[2].equals(strs2[2]) && Integer.valueOf(strs1[3]) > Integer.valueOf(strs2[3])) {
					outBedpe.write(bedpe1.toString());
					outBedpe.newLine();
				} else if (strs1[2].compareTo(strs2[2]) > 0) {
					outBedpe.write(bedpe1.toString());
					outBedpe.newLine();
				} else {
					outBedpe.write(bedpe2.toString());
					outBedpe.newLine();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    public String getStrand(String s) {
    	int x = Integer.valueOf(s);
		String strand = "+";
		if ((x & 16) != 0) {
			strand = "-";
		}
		return strand;
	}
}
