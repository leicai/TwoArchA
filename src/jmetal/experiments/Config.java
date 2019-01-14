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

	public static double getVolumeForReference(Solution ref) {
		double v = 1.0;
		for (int j = 0; j < ref.getNumberOfObjectives(); j++)
			v *= ref.getObjective(j);
		return v;
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
}
