package jmetal.metaheuristics.TwoArchA;

import java.io.FileNotFoundException;
import java.util.HashMap;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.experiments.Config;

import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.fastHypervolume.FastHypervolume;
import jmetal.util.JMException;
import jmetal.experiments.Config;

public class TwoArchA_Main {

	/**
	 * @param args
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	public static void main(String[] args) throws JMException, ClassNotFoundException {
		// TODO Auto-generated method stub
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use

		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters
		String probName = "WFG9";
		int obj = 10;

		problem = Config.setProblem(probName, obj);

		algorithm = new TwoArchA(problem);

		algorithm.setInputParameter("maxEvaluations", 100000);
		algorithm.setInputParameter("T", 20);
		algorithm.setInputParameter("delta", 0.9);
		algorithm.setInputParameter("nr", 2);

		int[] divs = Config.setDivs(obj);
		algorithm.setInputParameter("div1", divs[0]);
		algorithm.setInputParameter("div2", divs[1]);

		// Mutation and Crossover for Real codification
		parameters = new HashMap();
		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 20.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);


		long initTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		long estimatedTime = System.currentTimeMillis() - initTime;

		FastHypervolume fastHypervolume = new FastHypervolume();
		Solution referPoint = Config.getReferencePoint(probName, obj);
		double area = Config.getVolumeForReference(referPoint);
		Config.filterSolutions(population, referPoint);
		double meValue = fastHypervolume.computeHypervolume(population, referPoint) / area;
		System.out.println(meValue);
	}

}