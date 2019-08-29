package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public class Shell {
	
    public void runShell(String path) {
    	File f = new File(path);
        if (f.exists()) {
	        try {
				String command = "sh " + path;// notice the blank behind sh
				Process process = Runtime.getRuntime().exec(command);
				BufferedReader out = new BufferedReader(new InputStreamReader(process.getInputStream()));
				BufferedReader error = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				String message = null;
				while ((message = out.readLine()) != null) {
					System.out.println(message);
				}
				while ((message = error.readLine()) != null) {
					System.out.println(message);
				}
				process.waitFor();
				out.close();
				error.close();
	        } catch (IOException e) {
					e.printStackTrace();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
        } else {
        	System.out.println("Error: "+path+" doesn't exist");
        }
    }
}
