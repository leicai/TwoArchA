package jmetal.experiments;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.fastHypervolume.FastHypervolume;
//import jmetal.qualityIndicator.hypeHypervolume.HypEHypervolume;
import jmetal.util.JMException;
import jmetal.util.comparators.DominanceComparator;

public class Config {
	public static String Result_path="D:\\cygwin64\\home\\cailei\\Experiment_Data\\";
//	public static String Result_path="D:\\cygwin64\\home\\cailei\\Work\\workspace\\jmetal\\";
	public static Problem getProblem(String name, int objs) throws JMException {
		Problem problem = null;
		if (name.equals("DTLZ1")) {
			Object[] params = { "Real", objs + 4, objs };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.equals("DTLZ7")) {
			Object[] params = { "Real", objs + 19, objs };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.startsWith("DTLZ")) {
			Object[] params = { "Real", objs + 9, objs };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.startsWith("WFG")) {
			Object[] params = { "Real", objs - 1, 25 - objs, objs };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else {
			System.out.println("Error: function type " + name + " invalid");
			System.exit(-1);
		}
		return problem;
	}
	
	public static Problem setProblem(String name, int obj) throws JMException, ClassNotFoundException {
		Problem problem = null;
		if (name.equals("DTLZ1") || name.equals("SDTLZ1")) {
			Object[] params = { "Real", obj + 4, obj };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.equals("DTLZ7") || name.equals("NDTLZ7")) {
			Object[] params = { "Real", obj + 19, obj };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.startsWith("DTLZ") || name.equals("SDTLZ2")
				|| name.equals("CDTLZ2")) {
			Object[] params = { "Real", obj + 9, obj };
			problem = (new ProblemFactory()).getProblem(name, params);
		} else if (name.startsWith("WFG") ||name.startsWith("NWFG")) {
			Object[] params = { "Real", obj - 1, 25 - obj, obj };
			problem = (new ProblemFactory()).getProblem(name, params);
		}

		else {
			System.out.println("Error: function type " + name + " invalid");
			System.exit(-1);
		}
		return problem;
	}

	
	public static int getPopulationSize(int obj){
		if (obj == 3){
			return 91;
		}
		else if (obj == 4){
			return 166;
		}
		else if (obj == 5){
			return 105;
		}
		else if (obj == 6){
			return 182;
		}
		else if (obj == 7){
			return 238;
		}
		else if (obj == 8){
			return 240;
		}
		else if (obj == 9){
			return 210;
		}
		else if (obj == 10){
			return 276;
		}
		else if (obj == 13){
			return 182;
		}
		else if (obj == 15){
			return 136;
		}
		else 
			return 200;
	}
	
	public static int[] setDivs(int obj) {
		int divs[] = new int[2];

		if (obj == 2) {
			divs[0] = 99;
			divs[1] = 0;
		} else if (obj == 3) {
			divs[0] = 12;
			divs[1] = 0;
		} else if (obj == 4) {
			divs[0] = 8;
			divs[1] = 0;
		} else if (obj == 5) {
			divs[0] = 6;
			divs[1] = 0;
		} else if (obj == 6){
			divs[0] = 4;
			divs[1] = 3;
		} else if (obj == 7){
			divs[0] = 4;
			divs[1] = 2;
		} else if (obj == 8) {
			divs[0] = 3;
			divs[1] = 3;
		} else if (obj == 9) {
			divs[0] = 3;
			divs[1] = 2;
		} else if (obj == 10) {
			divs[0] = 3;
			divs[1] = 2;
		} else if (obj == 13) {
			divs[0] = 2;
			divs[1] = 2;
		} else if (obj == 15) {
			divs[0] = 2;
			divs[1] = 1;
		} else if (obj == 20) {
			divs[0] = 2;
			divs[1] = 1;
		}
		
		else {
			divs[0] = 90;
			divs[1] = 0;
		}
		return divs;
	}
	public static void filterSolutions(SolutionSet pop, Solution ref) {
		int i = 0;
		while (i < pop.size()) {
			Solution sol = pop.get(i);
			int res = new DominanceComparator().compare(sol, ref);
			if (res != -1)
				pop.remove(i);
			else
				i++;
		}
	}
	
	public static void printObject(Solution sol){
		for(int i=0;i<sol.getNumberOfObjectives();i++){
			System.out.print(sol.getObjective(i)+" ");
		}
		System.out.println();
	}
	

	
	public static void printVector(double[] a) {
		for(int i=0;i<a.length;i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println();
	}
	
	public static void printVector(int[] a) {
		for(int i=0;i<a.length;i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println();
	}
	

	



	public static Solution getReferencePoint(String name, int objs) {

		Solution ref = new Solution(objs);
		if (name.equals("DTLZ1")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 1);
		} else if (name.equals("DTLZ2") || name.equals("DTLZ4")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 2);
		} else if (name.equals("DTLZ3")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 3);
		} else if (name.equals("DTLZ5")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 4);
		} else if (name.equals("DTLZ6")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 11);
		} else if (name.equals("DTLZ7")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 21);
		} else if (name.startsWith("WFG")) {
			for (int i = 0; i < objs; i++)
				ref.setObjective(i, 2 * (i + 1) + 1);
			// ref.setObjective(i, 1);
		} else if (name.startsWith("NWFG")){
			for(int i=0;i<objs;i++){
				ref.setObjective(i, 2);
			}
		}
		else {

			System.out.println("Error: function type " + name + " invalid");
			System.exit(-1);

		}
		return ref;
	}
	
	public static void printFitness(SolutionSet population) {
		for(int i=0;i<population.size();i++) {
			Solution sol = population.get(i);
			System.out.print(sol.getFitness()+" ");
		}
		System.out.println();
	}


	public static String collectData(String path) throws IOException{
		File file=new File(path);
	    InputStreamReader osw = new InputStreamReader(new FileInputStream(file));
	    BufferedReader br      = new BufferedReader(osw);
	    String line="";
	    while((line=br.readLine())!=null){
	    	String[] strline =line.split(" ");
	    	if(strline[0].equals("Average:")){
	    		return strline[1];
	    	}
	    }
		return "NULL";
	}
	
	public static double[] getNadirPoint(String name, int obj) {
		double nadir[] = new double[obj];
		if (name.equals("DTLZ1")) {
			for (int j = 0; j < obj; j++)
				nadir[j] = 0.5;
		} else if (name.equals("DTLZ2") || name.equals("DTLZ3")
				|| name.equals("DTLZ4")) {
			for (int j = 0; j < obj; j++)
				nadir[j] = 1.0;
		} else if (name.equals("DTLZ7")) {
			for (int j = 0; j < obj - 1; j++)
				nadir[j] = 1.0;
			nadir[obj - 1] = 2 * obj;
		} else if (name.equals("SDTLZ1")) {
			if (obj == 3 || obj == 5) {
				for (int j = 0; j < obj; j++)
					nadir[j] = 0.5 * Math.pow(10, j);
			} else if (obj == 8) {
				for (int j = 0; j < obj; j++)
					nadir[j] = 0.5 * Math.pow(3, j);
			} else if (obj == 10) {
				for (int j = 0; j < obj; j++)
					nadir[j] = 0.5 * Math.pow(2, j);
			} else if (obj == 15) {
				for (int j = 0; j < obj; j++)
					nadir[j] = 0.5 * Math.pow(1.2, j);
			}

		} else if (name.equals("SDTLZ2")) {
			if (obj == 3 || obj == 5) {
				for (int j = 0; j < obj; j++)
					nadir[j] = Math.pow(10, j);
			} else if (obj == 8 || obj == 10) {
				for (int j = 0; j < obj; j++)
					nadir[j] = Math.pow(3, j);
			} else if (obj == 15) {
				for (int j = 0; j < obj; j++)
					nadir[j] = Math.pow(2, j);
			}
		} else if (name.startsWith("WFG")) {
			for (int j = 0; j < obj; j++)
				nadir[j] = 2 * (j + 1);
		} else if (name.startsWith("NWFG")) {
			for (int j = 0; j < obj; j++)
				nadir[j] = 1.0;
		} else if (name.startsWith("NDTLZ7")){
			for (int j = 0; j < obj; j++)
				nadir[j] = 1.0;
		}

		return nadir;
	}
	
	public static String getCenterFilePath(int num_subProblem,int num_obj) {
		return "centers/"+Integer.toString(num_subProblem)+"_"+Integer.toString(num_obj)+".txt";
	}
	
	public static double getRunTime(String filePath) throws IOException{
		FileInputStream fis = new FileInputStream(filePath);
		InputStreamReader isr = new InputStreamReader(fis);
		BufferedReader br = new BufferedReader(isr);
		
		String aux = br.readLine();
		String value=null;
		while (aux != null) {
			if (aux.startsWith("Time:")){
				value = aux.trim();
			}
			aux = br.readLine();
		}
		String[] parts = value.split(":");
		br.close();
		return Double.parseDouble(parts[1]);
	}
	
	public static List<Double> getMetricValueSet(String filePath, int times) throws IOException{
		List<Double> list = new ArrayList<Double>();
		FileInputStream fis = new FileInputStream(filePath);
		InputStreamReader isr = new InputStreamReader(fis);
		BufferedReader br = new BufferedReader(isr);
		
		String aux = br.readLine();
		int count = 0;
		while(aux!=null&count<times){
			if (!aux.startsWith("Average")){
				double val = Double.parseDouble(aux.trim());
				list.add(val);
				count++;
			}
			aux = br.readLine();
		}
		br.close();
		int size = list.size();
		
		int remain = times - size;
		int index = 0;
		
		while (remain > 0){
			double val = list.get(index);
			list.add(val);
			index = (index + 1) % size;
			remain--;
		}
		return list;
		
	}
	
	public static List<Double> getMetricValues(String path,String algName, String probName,int objs,int times,String metric) throws IOException{
		List<Double> list = new ArrayList<Double>();
		String filePath=null;
		if(metric.equals("HV")||metric.equals("HY")){
			 filePath = path+"/"+metric+"/"+algName+"/"+probName+"/"+objs+"-objective.txt";;
		}else if(metric.equals("HV1")){
			filePath = path+"/"+algName+"/"+probName+"/"+objs+"/HV1.txt";
		}
		else if(metric.equals("GD")){
			 filePath = path+"/"+algName+"/"+probName+"/"+objs+"/GD.txt";
		}
		
		FileInputStream fis = new FileInputStream(filePath);
		InputStreamReader isr = new InputStreamReader(fis);
		BufferedReader br = new BufferedReader(isr);
		
		String aux = br.readLine();
		
		int count = 0;
		while (aux != null & count < times) {
			if (!aux.startsWith("Average")){
				double val = Double.parseDouble(aux.trim());
				list.add(val);
				count++;
			}
			aux = br.readLine();
		}
		br.close();
		
		
		int size = list.size();
		
		int remain = times - size;
		int index = 0;
		
		while (remain > 0){
			double val = list.get(index);
			list.add(val);
			index = (index + 1) % size;
			remain--;
		}
		return list;
	}
	
	public static void printLatexTable(String[][] table) throws FileNotFoundException{
		int rows = table.length;
		int cols = table[0].length;
		
		for(int i=0;i<rows;i++){
			System.out.print("&");
			for(int j=0;j<cols;j++){
				if(j!=cols-1)
					System.out.print(table[i][j]+"&");
				else
					System.out.println(table[i][j]+"\\\\");
			}// for j
		}// for i
	}
	
	public static void printLatexTableWithHead(String[][] table,int[] heads) throws FileNotFoundException{
		int rows = table.length;
		int cols = table[0].length;
		
		for(int i=0;i<rows;i++){
			System.out.print("&");
			System.out.print(heads[i]);
			System.out.print("&");
			for(int j=0;j<cols;j++){
				if(j!=cols-1)
					System.out.print(table[i][j]+"&");
				else
					System.out.println(table[i][j]+"\\\\");
			}// for j
		}// for i
	}
	
	public static double[] getWorstPointFormFile(String path) throws IOException{
		SolutionSet population = readPopulationFromFile(path);
		return getWorstPoint(population);
	}
	
	public static double[] getWorstPoint(SolutionSet population){
		int objs = population.get(0).getNumberOfObjectives();
		double worstPoints[]=new double[objs];
		for(int i=0;i<objs;i++){
			worstPoints[i]=Double.MIN_VALUE;
		}
		for(int i=0;i<population.size();i++){
			for(int j=0;j<objs;j++){
				double value = population.get(i).getObjective(j);
				worstPoints[j] = worstPoints[j]>value?worstPoints[j]:value;
			}// for j
		}// for i
		return worstPoints;
	}
	
	public static void writeTableToLatexFile(String[][] table,String filepath) throws FileNotFoundException{
		//For latex table
		File file = new File(filepath);
		
		if(!file.exists()){
			String newpath=file.getParent();
			File newfile=new File(newpath);
			newfile.mkdirs();
		}
		
		PrintWriter out = new PrintWriter(new OutputStreamWriter(
				new FileOutputStream(file, false)));
		
		int rows=table.length; 			// acquire rows
		int cols=table[0].length;		// acquire columns 
		for(int i=0;i<rows;i++){
			out.print("&");
			for(int j=0;j<cols;j++){
				if(j!=cols-1)
					out.print(table[i][j]+"&");
				else
					out.println(table[i][j]+"\\\\");
			}
			
		}
		out.close();

	}
	
	public static void writeTableToFile(String[][] table,String filepath) throws FileNotFoundException{
		//For latex table
		File file = new File(filepath);
		
		if(!file.exists()){
			String newpath=file.getParent();
			File newfile=new File(newpath);
			newfile.mkdirs();
		}
		
		PrintWriter out = new PrintWriter(new OutputStreamWriter(
				new FileOutputStream(file, false)));
		
		int rows=table.length; 			// acquire rows
		int cols=table[0].length;		// acquire columns 
		for(int i=0;i<rows;i++){
			for(int j=0;j<cols;j++){
				if(j!=cols-1)
					out.print(table[i][j]+" ");
				else
					out.println(table[i][j]);
			}
			
		}
		out.close();

	}
	
	public static void writeTableToFile(double[][] table,String filepath) throws FileNotFoundException{
		//For latex table
		File file = new File(filepath);
		
		if(!file.exists()){
			String newpath=file.getParent();
			File newfile=new File(newpath);
			newfile.mkdirs();
		}
		
		PrintWriter out = new PrintWriter(new OutputStreamWriter(
				new FileOutputStream(file, false)));
		
		int rows=table.length; 			// acquire rows
		int cols=table[0].length;		// acquire columns 
		for(int i=0;i<rows;i++){
			for(int j=0;j<cols;j++){
				if(j!=cols-1)
					out.print(table[i][j]+" ");
				else
					out.println(table[i][j]);
			}
			
		}
		out.close();

	}
	
	public static double[][] readMatrixFromFile(String path,int rows, int cols)throws IOException{
		double[][] lambda = new double[rows][cols];
		File filename = new File(path); 
		InputStreamReader osw = new InputStreamReader(new FileInputStream(filename));
		BufferedReader br      = new BufferedReader(osw);
		String line="";
		int index_row = 0;
		while((line=br.readLine())!=null){
			String[] strSolutions =line.split(" ");
			double[] row = new double[cols];
			for(int i=0;i<strSolutions.length;i++){
				double p = Double.parseDouble(strSolutions[i]);
//				sol.setObjective(i,p);
				lambda[index_row][i] = p;
			}
			index_row++;

		}
		return lambda;
        
	}
	
	public static double[] readVectorFromFile(String path,int cols)throws IOException{

		double[] vector = new double[cols];
		File filename = new File(path); 
		InputStreamReader osw = new InputStreamReader(new FileInputStream(filename));
		BufferedReader br      = new BufferedReader(osw);
		String line="";

		if((line=br.readLine())!=null){
			String[] strSolutions =line.split(" ");
			double[] row = new double[cols];
			for(int i=0;i<strSolutions.length;i++){
				double p = Double.parseDouble(strSolutions[i]);
//				sol.setObjective(i,p);
				vector[i] = p;
			}
		}


		return vector;
        
	}
	public static SolutionSet readPopulationFromFile(String path,int capacity)throws IOException{
		// Read solutions form txt file.
		SolutionSet slt =new SolutionSet();
		slt.setCapacity(capacity);
		File filename = new File(path); 
	    InputStreamReader osw = new InputStreamReader(new FileInputStream(filename));
	    BufferedReader br      = new BufferedReader(osw);
	    String line="";

	    while((line=br.readLine())!=null){
		    String[] strSolutions =line.split(" ");
		    if(!strSolutions[0].equals("Time:")){// The last line of file is not need to load
		    	Solution sol =new Solution(strSolutions.length);
			    for(int i=0;i<strSolutions.length;i++){
			    	double p = Double.parseDouble(strSolutions[i]);
			    	sol.setObjective(i,p);
			    }
			    slt.add(sol);
		    }
	    }
	    return slt;
	}
	
	public static SolutionSet readPopulationFromFile(String path)throws IOException{
		// Read solutions form txt file.
		SolutionSet slt =new SolutionSet();
		File filename = new File(path); 
	    InputStreamReader osw = new InputStreamReader(new FileInputStream(filename));
	    BufferedReader br      = new BufferedReader(osw);
	    String line="";

	    while((line=br.readLine())!=null){
		    String[] strSolutions =line.split(" ");
		    if(!strSolutions[0].equals("Time:")){// The last line of file is not need to load
		    	Solution sol =new Solution(strSolutions.length);
			    for(int i=0;i<strSolutions.length;i++){
			    	double p = Double.parseDouble(strSolutions[i]);
			    	sol.setObjective(i,p);
			    }
			    slt.add(sol);
		    }
	    }
	    return slt;
	}
	
	public static void printMatrix(double[][] matrix){
		for(int i=0;i<matrix.length;i++){
			for (int j=0;j<matrix[i].length;j++){
				System.out.print(matrix[i][j]+" ");
			}
			System.out.println();
		}
	}
	

}
