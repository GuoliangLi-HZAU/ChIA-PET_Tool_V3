package process;

import java.util.ArrayList;

@SuppressWarnings("rawtypes")
public class SortK implements Comparable {
	
	private String content;// one row
    private ArrayList<String> columns = new ArrayList<String>();// the column to sort
    
    public SortK() {
    	
    }
    
    public void setContent(String content) {
    	this.content = content;
    }
    
    public void setColumns(String str) {
    	columns.add(str);
    }
    
    public String getContent() {
    	return content;
    }
    
    public String getColumns(int i) {
    	return columns.get(i);
    }
    
	public int compareTo(Object o) {
		int result = 0;
		for (int i = 0; i < columns.size(); i++) {
			result = columns.get(i).compareTo(((SortK)o).getColumns(i));
			if (result != 0) {
				break;
			}
		}
		return result;
	}
}
