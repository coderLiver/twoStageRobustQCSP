package ccg_first;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class Input {
	int numBay;
	int numQC;
	double traverseTime;
	double[] processTime;
	double uncertain;
	int[] initLocation;
	double budget;
	
	public void dataInput(String pathname) {
		try {
			File file = new File(pathname);
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String line = null;
			
			if((line = bufferedReader.readLine()) != null) {
				numBay = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				numQC = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				traverseTime = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1); //É¾³ýÊ×Î²×Ö·û
				String[] strings = line.split(",");
				processTime = new double[strings.length];
				for (int i = 0; i < strings.length; i++) {
					processTime[i] = Integer.parseInt(strings[i]);
				}
			}
			
			if((line = bufferedReader.readLine()) != null) {
				uncertain = Double.parseDouble(line);
			}
			
			if((line = bufferedReader.readLine()) != null) {
				budget = Double.parseDouble(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1);
				String[] strings = line.split(",");
				initLocation = new int[strings.length];
				for (int i = 0; i < strings.length; i++) {
					initLocation[i] = Integer.parseInt(strings[i]);
				}
			}
			
			bufferedReader.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
