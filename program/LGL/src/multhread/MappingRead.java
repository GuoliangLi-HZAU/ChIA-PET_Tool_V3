package multhread;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

public class MappingRead implements Runnable {

	private String samFile1;
	private String samFile2;
	private BlockingQueue<Object> queue;
	private Object END;
	
	public MappingRead(String samFile1, String samFile2, BlockingQueue<Object> queue, Object END) {
		this.samFile1 = samFile1;
		this.samFile2 = samFile2;
		this.queue = queue;
		this.END = END;
	}
	
	public void run() {
		try {
			BufferedReader reader1 = new BufferedReader(new FileReader(samFile1));
			BufferedReader reader2 = new BufferedReader(new FileReader(samFile2));
			String line1 = "";
			String line2 = "";
			ArrayList<String> list = new ArrayList<String>();
			int lineNum = 0;
			int mark = 0;
			while ((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null) {
				if (!line1.startsWith("@")) {
					list.add(line1);
					list.add(line2);
					lineNum++;
					mark = 1;
					if (lineNum % 10000 == 0) {
						queue.put(list);
						list = new ArrayList<String>();
						mark = 2;
					}
				}
			}
			if (mark == 1) {// put >10000 into queue
				queue.put(list);
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
