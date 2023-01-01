package ccg_unidirectional;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class StaticMethod {
	public static boolean intListEqual(int[] a, int[] b) {
		if (a.length != b.length) {
			return false;
		}else {
			for (int i = 0; i < b.length; i++) {
				if (a[i] != b[i]) {
					return false;
				}
			}
			return true;
		}
	}
	
	// 创建文件或在指定文件中写入内容，路径默认在./result of benchmark/下，可以创建一级文件夹
	// txtPath示例：./result of benchmark
	// txtName要以.txt结尾
	public static void txtFileCreatandWrite(String txtPath, String txtName, String writeContent) {
		File file = new File(txtPath);
		if(!file.exists()){//如果文件夹不存在
			file.mkdir();//创建文件夹
		}
		
		try {
			FileWriter writer = new FileWriter(txtPath + "/" + txtName, true);
			writer.write(writeContent);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// 重载，相比于前面的方法，该版本只用于创建新的空文件，不用于写入
	public static void txtFileCreatandWrite(String txtPath, String txtName) {
		File file = new File(txtPath);
		if(!file.exists()){//如果文件夹不存在
			file.mkdir();//创建文件夹
		}
		
		try {
			FileWriter writer = new FileWriter(txtPath + "/" + txtName);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// 清空文件夹中的所有文件
	public static boolean clearFolder(String txtPath) {
		File file = new File(txtPath);
		
		if(!file.exists()){//判断是否待删除目录是否存在
			System.err.println("The dir are not exists!");
			return false;
		}
		
		String[] content = file.list();//取得当前目录下所有文件和文件夹
		for(String name : content){
			File temp = new File(txtPath, name);
			
			if(temp.isDirectory()){//判断是否是目录
				clearFolder(temp.getAbsolutePath());//递归调用，删除目录里的内容
				temp.delete();//删除空目录
			}else{
				if(!temp.delete()){//直接删除文件
					System.err.println("Failed to delete " + name);
					System.exit(0);
				}
			}
		}
		return true;
	}
}
