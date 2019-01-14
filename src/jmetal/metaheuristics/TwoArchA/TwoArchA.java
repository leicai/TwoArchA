package jmetal.metaheuristics.TwoArchA;

//MOEAD.java
//
//Author:
//   Antonio J. Nebro <antonio@lcc.uma.es>
//   Juan J. Durillo <durillo@lcc.uma.es>
//
//Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

import jmetal.core.*;
import jmetal.metaheuristics.moead.Utils;
import jmetal.util.JMException;

import jmetal.util.PseudoRandom;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Comparator;
import java.util.StringTokenizer;
import java.util.Vector;

public class TwoArchA extends Algorithm {


	private int populationSize_;
	/**
	 * Stores the population
	 */
	private SolutionSet CA_;
	private SolutionSet DA_;
	/**
	 * Z vector (ideal point)
	 */
	double[] z_;
	/**
	 * Lambda vectors
	 */
	// Vector<Vector<Double>> lambda_ ;
	double[][] lambda_;
	/**
	 * T: neighbour size
	 */
	int T_;
	/**
	 * Neighborhood
	 */
	int[][] neighborhood_;
	/**
	 * delta: probability that parent solutions are selected from neighbourhood
	 */
	double delta_;
	/**
	 * nr: maximal number of solutions replaced by each child solution
	 */
	int nr_;
	Solution[] indArray_;
	String functionType_;
	int generations_;
	/**
	 * Operators
	 */
	Operator crossover_;
	Operator mutation_;

	int div1_;
	int div2_;


	int evaluations_;
	
	
	public TwoArchA(Problem problem) {
		super(problem);

		functionType_ = "_TCHE";

	} // DMOEA

	public SolutionSet execute() throws JMException, ClassNotFoundException {

		
		int maxEvaluations;
		maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations"))
				.intValue();
		evaluations_ = 0;

		div1_ = ((Integer) this.getInputParameter("div1")).intValue();
		div2_ = ((Integer) this.getInputParameter("div2")).intValue();
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();
		
	//	System.out.println(lambda_.length);

		populationSize_ = vg.getVectors().length;

		CA_ = new SolutionSet(populationSize_);
		DA_ = new SolutionSet(populationSize_);
		
		indArray_ = new Solution[problem_.getNumberOfObjectives()];

		T_ = ((Integer) this.getInputParameter("T")).intValue();
		nr_ = ((Integer) this.getInputParameter("nr")).intValue();
		delta_ = ((Double) this.getInputParameter("delta")).doubleValue();

		neighborhood_ = new int[populationSize_][T_];

		z_ = new double[problem_.getNumberOfObjectives()];

		crossover_ = operators_.get("crossover"); 
		mutation_ = operators_.get("mutation"); 


		initNeighborhood();

		// STEP 1.2. Initialize population
		initPopulation();

		// STEP 1.3. Initialize z_
		initIdealPoint();

		// STEP 2. Update
		do {
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);

			for (int i = 0; i < populationSize_; i++) {
				int n = permutation[i]; // or int n = i;
				// int n = i ; // or int n = i;
				int type;
				double rnd = PseudoRandom.randDouble();

				// STEP 2.1. Mating selection based on probability
				if (rnd < delta_) // if (rnd < realb)
				{
					type = 1; // neighborhood
				} else {
					type = 2; // whole population
				}
				
				
				Vector<Integer> p = new Vector<Integer>();
				matingSelection(p, n, 2, type);

				// STEP 2.2. Reproduction
				Solution child;
				Solution[] parents = new Solution[2];

				parents[0] = CA_.get(p.get(0));
				parents[1] = DA_.get(p.get(1));

				Solution[] offSpring = (Solution[]) crossover_.execute(parents);

				child = offSpring[0];


				// Apply mutation
				mutation_.execute(child);

				// Evaluation
				problem_.evaluate(child);
				
				evaluations_++;

				updateReference(child);
				
				updateCA(child, n, type);
				updateDA(child, n, type);
	
			} // for
			generations_++;

		} while (evaluations_ < maxEvaluations);

		return DA_;
	}



	public void initNeighborhood() {
		double[] x = new double[populationSize_];
		int[] idx = new int[populationSize_];


		for (int i = 0; i < populationSize_; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize_; j++) {
				x[j] = Utils.distVector(lambda_[i], lambda_[j]);
				// x[j] = dist_vector(population[i].namda,population[j].namda);
				idx[j] = j;
				// System.out.println("x["+j+"]: "+x[j]+
				// ". idx["+j+"]: "+idx[j]) ;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize_, T_);
			// minfastsort(x,idx,population.size(),niche);

			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	/**
* 
*/
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);

			evaluations_++;
			
			CA_.add(newSolution);
			DA_.add(new Solution(newSolution));
		} // for
	} // initPopulation

	/**
* 
*/
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			z_[i] = 1.0e+30;
			indArray_[i] = new Solution(problem_);
			problem_.evaluate(indArray_[i]);

		} // for

		for (int i = 0; i < populationSize_; i++) {
			updateReference(CA_.get(i));
		} // for
		for (int i = 0; i < populationSize_; i++) {
			updateReference(DA_.get(i));
		} // for
	} // initIdealPoint

	/**
* 
*/
	public void matingSelection(Vector<Integer> list, int cid, int size,
			int type) {
		// list : the set of the indexes of selected mating parents
		// cid : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			// if (flag) list.push_back(p);
			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection
	
	
	/**
	 * 
	 * @param individual
	 */
	void updateReference(Solution individual) {
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < z_[n]) {
				z_[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference

	/**
	 * @param individual
	 * @param id
	 * @param type
	 */
	
	void updateDA(Solution indiv, int id, int type){
		// indiv: child solution
		// id: the id of current subproblem
		// type: update solutions in - neighborhood (1) or whole population
		// (otherwise)
		int size;
		int time;
		
		time = 0;
		if (type == 1) {
			size = neighborhood_[id].length;
		} else {
			size = DA_.size();
		}
		int[] perm = new int[size];

		Utils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[id][perm[i]];
			} else {
				k = perm[i]; // calculate the values of objective function
								// regarding the current subproblem
			}
			double[] f1, f2;

//			f1 = fitnessFunction(DA_.get(k), lambda_[k]);
//			f2 = fitnessFunction(indiv, lambda_[k]);
			f1 = getDistances(DA_.get(k), lambda_[k]);
			f2 = getDistances(indiv,lambda_[k]);

			if (f2[1] < f1[1]) {
				DA_.replace(k, new Solution(indiv));
				// population[k].indiv = indiv;
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed
			// 'limit'
			if (time >= nr_) {
				return;
			}
		}
	}

	void updateCA(Solution indiv, int id, int type) {
		// indiv: child solution
		// id: the id of current subproblem
		// type: update solutions in - neighborhood (1) or whole population
		// (otherwise)
		int size;
		int time;

		time = 0;

		if (type == 1) {
			size = neighborhood_[id].length;
		} else {
			size = CA_.size();
		}
		int[] perm = new int[size];

		Utils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[id][perm[i]];
			} else {
				k = perm[i]; // calculate the values of objective function
								// regarding the current subproblem
			}
			double f1, f2;

			f1 = fitnessFunction(CA_.get(k), lambda_[k]);
			f2 = fitnessFunction(indiv, lambda_[k]);

			if (f2 < f1) {
				CA_.replace(k, new Solution(indiv));
				// population[k].indiv = indiv;
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed
			// 'limit'
			if (time >= nr_) {
				return;
			}
		}
	} // updateProblem


	double fitnessFunction(Solution individual, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE")) {
			double maxFun = -1.0e+30;

			double sum = 0;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(individual.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}

				sum += feval;
			} // for

			fitness = maxFun;
		} // if
		else if (functionType_.equals("_PBI")) {

			double d1, d2, nl;
			// double theta = 1.0 / Math.tan(alpha_);

			double theta = 5.0;
			// System.out.println(theta);

			d1 = d2 = nl = 0.0;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - z_[i]) * lambda[i];

				nl += (lambda[i] * lambda[i]);
			}
			nl = Math.sqrt(nl);
			d1 = Math.abs(d1) / nl;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {

				d2 += ((individual.getObjective(i) - z_[i]) - d1
						* (lambda[i] / nl))
						* ((individual.getObjective(i) - z_[i]) - d1
								* (lambda[i] / nl));
			}
			d2 = Math.sqrt(d2);

			fitness = (d1 + theta * d2);
			// fitness = d1 + 1e6 * d2;
			// fitness = d2;

		} else if (functionType_.equals("_VADS")) {
			double inter = 0.0;
			double dx = 0.0;
			double dy = 0.0;
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				dx += (individual.getObjective(i) * individual.getObjective(i));
				dy += (lambda[i] * lambda[i]);
				inter += (individual.getObjective(i) * lambda[i]);
			}
			dx = Math.sqrt(dx);
			dy = Math.sqrt(dy);

			// double cos = inter / (dx * dy);

			int q = 100;
			double val = (q + 1) * Math.log(dx) - q * Math.log(inter / dy);

			fitness = Math.exp(val);

		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation} // MOEAD
	
	double[] getDistances(Solution individual, double[] lambda) {
		double d1, d2, nl;

		d1 = d2 = nl = 0.0;

		for (int i = 0; i < individual.getNumberOfObjectives(); i++) {
			d1 += (individual.getObjective(i) - z_[i]) * lambda[i];

			nl += (lambda[i] * lambda[i]);
		}
		nl = Math.sqrt(nl);
		d1 = Math.abs(d1) / nl;

		for (int i = 0; i < individual.getNumberOfObjectives(); i++) {

			d2 += ((individual.getObjective(i) - z_[i]) - d1 * (lambda[i] / nl))
					* ((individual.getObjective(i) - z_[i]) - d1
							* (lambda[i] / nl));
		}
		d2 = Math.sqrt(d2);

		double ds[] = new double[2];
		ds[0] = d1;
		ds[1] = d2;
		return ds;

	}

} 
