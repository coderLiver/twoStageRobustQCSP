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
	
	// �����ļ�����ָ���ļ���д�����ݣ�·��Ĭ����./result of benchmark/�£����Դ���һ���ļ���
	// txtPathʾ����./result of benchmark
	// txtNameҪ��.txt��β
	public static void txtFileCreatandWrite(String txtPath, String txtName, String writeContent) {
		File file = new File(txtPath);
		if(!file.exists()){//����ļ��в�����
			file.mkdir();//�����ļ���
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
	
	// ���أ������ǰ��ķ������ð汾ֻ���ڴ����µĿ��ļ���������д��
	public static void txtFileCreatandWrite(String txtPath, String txtName) {
		File file = new File(txtPath);
		if(!file.exists()){//����ļ��в�����
			file.mkdir();//�����ļ���
		}
		
		try {
			FileWriter writer = new FileWriter(txtPath + "/" + txtName);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// ����ļ����е������ļ�
	public static boolean clearFolder(String txtPath) {
		File file = new File(txtPath);
		
		if(!file.exists()){//�ж��Ƿ��ɾ��Ŀ¼�Ƿ����
			System.err.println("The dir are not exists!");
			return false;
		}
		
		String[] content = file.list();//ȡ�õ�ǰĿ¼�������ļ����ļ���
		for(String name : content){
			File temp = new File(txtPath, name);
			
			if(temp.isDirectory()){//�ж��Ƿ���Ŀ¼
				clearFolder(temp.getAbsolutePath());//�ݹ���ã�ɾ��Ŀ¼�������
				temp.delete();//ɾ����Ŀ¼
			}else{
				if(!temp.delete()){//ֱ��ɾ���ļ�
					System.err.println("Failed to delete " + name);
					System.exit(0);
				}
			}
		}
		return true;
	}
}
