

/* count.cpp -- Master program

 Copyright 2002, 2003 Raymond Hemmecke, Ruriko Yoshida
 Copyright 2006, 2007 Matthias Koeppe

 This file is part of LattE.

 LattE is free software; you can redistribute it and/or modify it
 under the terms of the version 2 of the GNU General Public License
 as published by the Free Software Foundation.

 LattE is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with LattE; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */



#include "count.h"


/**
 * the polynomial series sometimes is in the form 0 + num*t
 * because not all functions in the mainCountDriver return prober polynomials
 * encoding the lattice point count.
 */
void CountAnswerContainer::checkPolynomial()
{
	if ( seriesExpansion.length() == 2
			&& seriesExpansion[0] == 0)
	{
		assert(numLaticePoints == 0 || numLaticePoints == seriesExpansion[1]);
		numLaticePoints = seriesExpansion[1];
		seriesExpansion.kill();
	}//if
}//checkPolynomial




/* ----------------------------------------------------------------- */
CountAnswerContainer mainCountDriver(int argc, char *argv[])
{
	set_program_name(argv[0]);
	CountAnswerContainer countAnswerContainer;

	int i;
	unsigned int flags = 0, print_flag = 0, output_cone = 0;
	char printfile[127], Save_Tri[127], Load_Tri[127], Print[127],
		removeFiles[127];
	bool approx;
	bool ehrhart_polynomial, ehrhart_series, ehrhart_taylor;
	bool triangulation_specified = false;
	double sampling_factor = 1.0;
	long int num_samples = -1;
	ReadPolyhedronData read_polyhedron_data;

	struct BarvinokParameters *params = new BarvinokParameters;

	latte_banner(cerr);

	if (argc < 2)
	{
	  cerr << "usage: " << shell_quote(argv[0]) << " [OPTIONS...] " << "INPUTFILE" << endl;
	  cerr << "Type `" << shell_quote(argv[0]) << " --help' \n"
	       << "for a list of options and input specifications." << endl;
	  THROW_LATTE_MSG(LattException::ue_BadCommandLineOptionCount, "too few options used");
	}

	//setbuf(stdout,0);

	cerr << "Invocation: ";
	for (i = 0; i < argc; i++)
	{
	  cerr << shell_quote(argv[i]) << " ";
	}
	cerr << endl;

	strcpy(removeFiles, "yes");
	strcpy(printfile, "no");
	strcpy(Save_Tri, "no");
	strcpy(Load_Tri, "no");
	strcpy(Print, "no");
	approx = false;
	ehrhart_polynomial = false;
	params->substitution = BarvinokParameters::PolynomialSubstitution;
	params->decomposition = BarvinokParameters::DualDecomposition;
	params->max_determinant = 1;

	for (i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "nodecom", 3) == 0 || strncmp(argv[i],
				"--nodecomposition", 5) == 0 || strncmp(argv[i],
				"--no-decomposition", 7) == 0)
			params->max_determinant = 0;
		else if (strncmp(argv[i], "uni", 3) == 0)
			strcpy(read_polyhedron_data.assumeUnimodularCones, "yes");
		//else if(strncmp(argv[i],"simp",4)==0) {strcpy(printfile,"yes"); flags |= PRINT;}
		else if (strncmp(argv[i], "file", 4) == 0)
			strcpy(read_polyhedron_data.Memory_Save, "no");
		//else if(strncmp(argv[i],"single",6)==0) strcpy(Singlecone,"yes");
		//else if(strncmp(argv[i],"ehrhartsimp",3)==0) strcpy(rationalCone,"yes");
		else if (strncmp(argv[i], "memsave", 7) == 0)
			strcpy(read_polyhedron_data.Memory_Save, "yes");
		else if (strncmp(argv[i], "printcones", 3) == 0)
			strcpy(Print, "yes");
		//else if(strncmp(argv[i],"hull",3)==0) strcpy (inthull, "yes");
		else if (strncmp(argv[i], "rem", 3) == 0)
		{
			strcpy(removeFiles, "no");
			strcpy(read_polyhedron_data.Memory_Save, "no");
		} else if (strncmp(argv[i], "trisave", 7) == 0)
		{
			strcpy(Save_Tri, "yes");
			flags |= SAVE;
		} else if (strncmp(argv[i], "triload", 7) == 0)
		{
			strcpy(Load_Tri, "yes");
			flags |= LOAD;
		} else if (strncmp(argv[i], "--exponential", 5) == 0)
			params->substitution = BarvinokParameters::ExponentialSubstitution;
		else if (strncmp(argv[i], "--polynomial", 6) == 0)
			params->substitution = BarvinokParameters::PolynomialSubstitution;
		else if (strncmp(argv[i], "--maxdet=", 9) == 0)
			params->max_determinant = atoi(argv[i] + 9);
		else if (strncmp(argv[i], "--irrational-all-primal", 14) == 0
				|| strncmp(argv[i], "--all-primal", 5) == 0)
			params->decomposition
					= BarvinokParameters::IrrationalAllPrimalDecomposition;
		else if (strncmp(argv[i], "--irrational-primal", 5) == 0)
			params->decomposition
					= BarvinokParameters::IrrationalPrimalDecomposition;
		else if (strcmp(argv[i], "--dual") == 0) // Don't use strncmp to
			// avoid clash with --dualization=...
			params->decomposition = BarvinokParameters::DualDecomposition;
		else if (strncmp(argv[i], "--count-lattice-points", 7) == 0)
		{
			// Default.
		} else if (strncmp(argv[i], "--multivariate-generating-function", 7)
				== 0)
		{
			params->substitution = BarvinokParameters::NoSubstitution;
		} else if (strncmp(argv[i], "--ehrhart-polynomial", 11) == 0)
		{
			ehrhart_polynomial = true;
			params->substitution = BarvinokParameters::ExponentialSubstitution;
		} else if (strncmp(argv[i], "--ehrhart-series", 11) == 0)
		{
			ehrhart_series = true;
			strcpy(read_polyhedron_data.dualApproach, "yes");
			strcpy(printfile, "yes");
			flags |= PRINT;
		} else if (strncmp(argv[i], "--simplified-ehrhart-series", 14) == 0)
		{
			ehrhart_series = true;
			strcpy(read_polyhedron_data.dualApproach, "yes");
			strcpy(read_polyhedron_data.rationalCone, "yes");
		} else if (strncmp(argv[i], "--ehrhart-taylor=", 17) == 0)
		{
			strcpy(read_polyhedron_data.taylor, "yes");
			read_polyhedron_data.degree = atoi(argv[i] + 17);
			strcpy(read_polyhedron_data.dualApproach, "yes");
		} else if (strncmp(argv[i], "--avoid-singularities", 7) == 0)
		{
			params->shortvector = BarvinokParameters::SubspaceAvoidingLLL;
		} else if (parse_standard_triangulation_option(argv[i], params))
		{
			if (strncmp(argv[i], "--triangulation=", 16) == 0)
				triangulation_specified = true;
		} else if (parse_standard_dualization_option(argv[i], params))
		{
		} else if (parse_standard_smith_option(argv[i], params))
		{
		} else if (strncmp(argv[i], "--approximate", 7) == 0)
			approx = true;
		else if (strncmp(argv[i], "--sampling-factor=", 18) == 0)
			sampling_factor = atof(argv[i] + 18);
		else if (strncmp(argv[i], "--num-samples=", 14) == 0)
			num_samples = atol(argv[i] + 14);
		else if (strncmp(argv[i], "--random-seed=", 14) == 0)
		{
			unsigned int seed = atoi(argv[i] + 14);
			seed_random_generator(seed);
		} else if (strcmp(argv[i], "--help") == 0)
		{
			read_polyhedron_data.show_options(cerr);
			cerr << "Options that control what to compute:" << endl
					<< "  --count-lattice-points                   Compute the number of lattice points"
					<< endl
					<< "                                           (default)"
					<< endl
					<< "  --multivariate-generating-function       Compute the multivariate generating function of"
					<< endl
					<< "                                           the set of lattice points of the polyhedron"
					<< endl
					<< "  --ehrhart-polynomial                     Compute an Ehrhart polynomial of an integral polytope"
					<< endl
					<< "  --ehrhart-series                         Compute the unsimplified Ehrhart series"
					<< endl
					<< "                                           as a univariate rational function"
					<< endl
					<< "  --simplified-ehrhart-series              Compute the simplified Ehrhart series"
					<< endl
					<< "                                           as a univariate rational function"
					<< endl
					<< "  --ehrhart-taylor=N                       Compute the first N terms of the Ehrhart series"
					<< endl;
			cerr << "Options for the Barvinok algorithm:" << endl
					<< "  --dual                                   Triangulate and signed-decompose in the dual space"
					<< endl
					<< "                                           (traditional method, default)"
					<< endl
					<< "  --irrational-primal                      Triangulate in the dual space, signed-decompose"
					<< endl
					<< "                                           in the primal space using irrationalization"
					<< endl
					<< "  --irrational-all-primal                  Triangulate and signed-decompose in the primal space"
					<< endl
					<< "                                           using irrationalization"
					<< endl
					<< "  --maxdet=N                               Decompose down to an index (determinant) of N"
					<< endl
					<< "                                           instead of index 1 (unimodular cones)"
					<< endl
					<< "  --no-decomposition                       Do not signed-decompose simplicial cones"
					<< endl;
			cerr << "Options for specialization:" << endl
					<< "  --polynomial                             Use polynomial substitution for specialization"
					<< endl
					<< "                                           (traditional method, default)"
					<< endl
					<< "  --exponential                            Use exponential substitution for specialization"
					<< endl
					<< "                                           (recommended for maxdet > 1)"
					<< endl;
			cerr << "Algorithmic options for subproblems:" << endl;
			show_standard_smith_option(cerr);
			show_standard_dualization_option(cerr);
			show_standard_triangulation_options(cerr);
			exit(0);//THROW_LATTE(LattException::ue_HelpMenuDisplayed);
		} else if (read_polyhedron_data.parse_option(argv[i]))
		{
		} else
		{
			cerr << "Unknown command/option " << argv[i] << endl;
			THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, argv[i]);
		}
	}//for i.

	if (read_polyhedron_data.expect_filename)
	{
		cerr << "Filename missing" << endl;
		THROW_LATTE(LattException::ue_FileNameMissing);
	}

	if (params->shortvector == BarvinokParameters::SubspaceAvoidingLLL)
	{
		if (params->decomposition
				== BarvinokParameters::IrrationalAllPrimalDecomposition)
		{
			/* Triangulation will be done in the primal space, so all
			 triangulation methods are fine. */
		} else
		{
			/* Triangulation will be done in the dual space, so we must
			 avoid using facets whose normal vectors lie in the
			 subspace. */
			if (triangulation_specified)
			{
				if (params->triangulation
						!= BarvinokParameters::SubspaceAvoidingBoundaryTriangulation
						&& params->triangulation
								!= BarvinokParameters::SubspaceAvoidingSpecialTriangulation)
				{
					cerr
							<< "Warning: The requested triangulation method is not guaranteed to work with --avoid-singularities."
							<< endl;
				}
			} else
			{
				/* Not specified, so choose one that will work. */
				cerr << "Choosing SubspaceAvoidingBoundaryTriangulation method"
						<< endl;
				params->triangulation
						= BarvinokParameters::SubspaceAvoidingBoundaryTriangulation;
			}
		}
	}

	if (approx)
	{
		params->substitution = BarvinokParameters::ExponentialSubstitution;
		if (params->decomposition == BarvinokParameters::DualDecomposition)
		{
			cerr
					<< "Exponential approximation not implemented for dual decomposition; switching to irrational primal decomposition."
					<< endl;
			params->decomposition
					= BarvinokParameters::IrrationalPrimalDecomposition;
		}
	}

	if (read_polyhedron_data.minimize[0] == 'y')
		strcpy(read_polyhedron_data.maximum, "yes");
	if (read_polyhedron_data.grobner[0] == 'y')
		strcpy(read_polyhedron_data.equationsPresent, "yes");
	if (read_polyhedron_data.maximum[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (printfile[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (read_polyhedron_data.rationalCone[0] == 'y')
		strcpy(read_polyhedron_data.Memory_Save, "no");
	if (printfile[0] == 'y')
		print_flag = 1;

	if (read_polyhedron_data.rationalCone[0] == 'y')
	{

		//HugInt digit(argv[1]);
		//conv(output_cone, digit.BigInt);
		// User can use only Mode one
		output_cone = 3;
	}

	if (output_cone > 3)
		output_cone = 0;
	flags |= (output_cone << 1);

	const char *fileName = read_polyhedron_data.filename.c_str();

	if (read_polyhedron_data.dualApproach[0] == 'y')
	{
		flags |= DUAL_APPROACH;
	}

	/* INPUT HANDLING. */

	if (read_polyhedron_data.grobner[0] == 'y')
	{
		CheckGrobner(fileName, read_polyhedron_data.cddstyle);
		countAnswerContainer.numLaticePoints = SolveGrobner(fileName, read_polyhedron_data.nonneg,
				read_polyhedron_data.dualApproach,
				read_polyhedron_data.grobner,
				read_polyhedron_data.equationsPresent,
				read_polyhedron_data.cddstyle);
		return countAnswerContainer;
	}

	Polyhedron *Poly = read_polyhedron_data.read_polyhedron(params);

	/* Compute the facet information from tight inequalities if
	 possible.  It is essential that this is done BEFORE translating
	 vertexes to the origin (in Ehrhart polynomial mode) -- otherwise
	 tightness information is wrong. */
	if (not Poly->dualized && Poly->cones != NULL
			&& read_polyhedron_data.matrix != NULL
			&& read_polyhedron_data.Vrepresentation[0] != 'y')
	{
		/* Fill in the facets of all cones; we determine them by
		 taking all inequalities tight at the respective vertex. */
		params->dualize_time.start();
		computeTightInequalitiesOfCones(Poly->cones,
				read_polyhedron_data.matrix, Poly->numOfVars);
		params->dualize_time.stop();
		cerr << params->dualize_time;
	}

	if (ehrhart_polynomial)
	{
		/* Translate all cones to the origin, saving the original vertex. */
		listCone *cone;
		for (cone = Poly->cones; cone; cone = cone->rest)
		{
			ZZ scale_factor;
			cone->vertex->ehrhart_vertex = scaleRationalVectorToInteger(
					cone->vertex->vertex, Poly->numOfVars, scale_factor);
			if (scale_factor != 1)
			{
				cerr
						<< "Computation of Ehrhart (quasi-)polynomials is only implemented "
						<< "for integral polytopes." << endl
						<< "Use `--ehrhart-series' or `--simplfied-ehrhart-series' for computing "
						<< "the Ehrhart series of rational polytopes." << endl;
				THROW_LATTE_MSG(LattException::pe_RationalPolytope, "Ehrhart (quasi-)polynomials is only implemented for integral polytopes");
			}
			delete cone->vertex->vertex;
			cone->vertex->vertex = new rationalVector(Poly->numOfVars);
		}
	}

	params->Flags = flags;
	params->File_Name = (char*) fileName;
	params->Number_of_Variables = Poly->numOfVars;

	switch (params->decomposition)
	{
	case BarvinokParameters::DualDecomposition:
	case BarvinokParameters::IrrationalPrimalDecomposition:
		if (not Poly->dualized)
		{
			if (read_polyhedron_data.Vrepresentation[0] != 'y')
			{
				/* Compute all inequalities tight at the respective vertex.
				 Then dualizeCones just needs to swap rays and facets. */
				computeTightInequalitiesOfCones(Poly->cones,
						read_polyhedron_data.matrix, Poly->numOfVars);
			}
			dualizeCones(Poly->cones, Poly->numOfVars, params);
			Poly->dualized = true;
		}
		break;
	case BarvinokParameters::IrrationalAllPrimalDecomposition:
		cerr << "Irrationalizing polyhedral cones... ";
		cerr.flush();
		if (Poly->dualized)
		{
			cerr << "(First dualizing back... ";
			cerr.flush();
			dualizeCones(Poly->cones, Poly->numOfVars, params);
			cerr << "done; sorry for the interruption.) ";
			cerr.flush();
		} else if (Poly->cones != NULL)
		{
			if (Poly->cones->facets == NULL)
			{
				cerr << "(First computing facets for them... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
				cerr << "done; sorry for the interruption.) ";
				cerr.flush();
			} else if (Poly->cones->rays == NULL)
			{
				/* Only facets computed, for instance by using the 4ti2
				 method of computing vertex cones.  So dualize twice to
				 compute the rays. */
				cerr << "(First computing their rays... ";
				cerr.flush();
				dualizeCones(Poly->cones, Poly->numOfVars, params);
				dualizeCones(Poly->cones, Poly->numOfVars, params); // just swaps
				cerr << "done; sorry for the interruption.) ";
				cerr.flush();
			}
		}
		params->irrationalize_time.start();
		{
			listCone *cone;
			for (cone = Poly->cones; cone; cone = cone->rest)
				assert(lengthListVector(cone->facets) >= Poly->numOfVars);
		}
		irrationalizeCones(Poly->cones, Poly->numOfVars);
		params->irrationalize_time.stop();
		cerr << params->irrationalize_time;
		break;
	default:
		cerr << "Unknown BarvinokParameters::decomposition" << endl;
		abort();
	}

	//start of main task.
	try
	{
		switch (params->substitution)
		{
		case BarvinokParameters::NoSubstitution:
		{

			//multivariate generating function is printed.

			string rat_filename = read_polyhedron_data.filename + ".rat";
			DelegatingSingleConeParameters *write_params =
					new DelegatingSingleConeParameters(*params);
			delete params;
			params = write_params;
			ConeConsumer *writing_consumer =
					new GeneratingFunctionWritingConeConsumer(rat_filename);

			if (Poly->projecting_up_transducer)
				writing_consumer = compose(Poly->projecting_up_transducer,
						writing_consumer);
			writing_consumer = compose(
					new PointsInParallelepipedComputingConeTransducer(
							write_params), writing_consumer);
			write_params->SetConsumer(writing_consumer);
			cerr << "Computing multivariate generating function" << endl;

			listCone *cone;
			for (cone = Poly->cones; cone != NULL; cone = cone->rest)
				barvinokDecomposition_Single(cone, write_params);
			cerr << "Multivariate generating function written to `"<< rat_filename.c_str() << "'" << endl;
			countAnswerContainer.multivariateGenFunctionFileName = rat_filename;
			break;
		}
		case BarvinokParameters::PolynomialSubstitution:

			if (ehrhart_polynomial)
			{
				cerr
						<< "Computation of Ehrhart polynomials is only implemented "
						<< "for the exponential substitution (--exp)." << endl;
				THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, "Must use --exp for Ehrhart polynomials");
			}
			if (Poly->unbounded)
			{
				cerr << "The polyhedron is unbounded." << endl;
				THROW_LATTE(LattException::pe_Unbounded);
			}
			if (read_polyhedron_data.assumeUnimodularCones[0] == 'n')
			{
				if (read_polyhedron_data.Memory_Save[0] == 'n')
				{
					//I'm not yet sure how this computation is used.
					//I think --ehrhart-series uses this.

					listCone *decomposed_cones = decomposeCones(Poly->cones,
							not Poly->dualized, *params);
					freeListCone(Poly->cones);
					Poly->cones = decomposed_cones;
					// 	cerr << "Decomposed cones: " << endl;
					//	printListCone(Poly->cones, Poly->numOfVars);
					/* Compute points in parallelepipeds */
					computePointsInParallelepipeds(Poly->cones,
							Poly->numOfVars, params);
				}
				// Iterator through simplicial cones, DFS
				else
				{
					//prints the Taylor Expansion
					Standard_Single_Cone_Parameters *standard_params =
							new Standard_Single_Cone_Parameters(*params);
					delete params;
					params = standard_params;
					countAnswerContainer.seriesExpansion = decomposeAndComputeResidue(Poly->cones,
							read_polyhedron_data.degree, false,
							*standard_params);

				}
			}
			break;
		case BarvinokParameters::ExponentialSubstitution:

			//find number_of_lattice_points or the ehrhart polynomial.
			if (Poly->unbounded)
			{
				cerr << "The polyhedron is unbounded." << endl;
				THROW_LATTE(LattException::pe_Unbounded);
			}
			if (read_polyhedron_data.dualApproach[0] == 'y')
			{
				cerr
						<< "Exponential substitution is not yet implemented for the homogenized version."
						<< endl;
				THROW_LATTE_MSG(LattException::ue_BadCommandLineOption, "Exponential substitution is not implemented for homogenized polytopes");
			} else
			{
				if (approx)
				{
#ifdef HAVE_EXPERIMENTS
					{
						Write_Exponential_Sample_Formula_Single_Cone_Parameters *write_param
						= new Write_Exponential_Sample_Formula_Single_Cone_Parameters
						(*params, "Exponential_Sample_Formula", sampling_factor,
								num_samples);
						delete params;
						params = write_param;
						decomposeAndWriteExponentialSampleFormula(Poly->cones, *write_param);
					}
#else
					cerr << "Approximation code is not compiled in, sorry."
							<< endl;
					THROW_LATTE(LattException::ue_BadCommandLineOption);
#endif
				} else if (ehrhart_polynomial)
				{
					Exponential_Ehrhart_Parameters *exp_param =
							new Exponential_Ehrhart_Parameters(*params);
					delete params;
					params = exp_param;
					mpq_vector ehrhart_coefficients =
							decomposeAndComputeEhrhartPolynomial(Poly->cones,
									*exp_param);
					cerr << endl << "Ehrhart polynomial: ";
					{
						unsigned int i;
						for (i = 0; i < ehrhart_coefficients.size(); i++)
						{
							if (ehrhart_coefficients[i] > 0)
								cout << " + " << ehrhart_coefficients[i]
										<< " * t^" << i;
							else if (ehrhart_coefficients[i] < 0)
								cout << " - " << abs(ehrhart_coefficients[i])
										<< " * t^" << i;
						}
					}
					countAnswerContainer.ehrhart_coefficients = ehrhart_coefficients;
					cout << endl << endl;
				} else
				{
					Exponential_Single_Cone_Parameters *exp_param =
							new Exponential_Single_Cone_Parameters(*params);
					delete params;
					params = exp_param;
					Integer number_of_lattice_points =
							decomposeAndComputeExponentialResidue(Poly->cones,
									*exp_param);
					params->deliver_number_of_lattice_points(
							number_of_lattice_points);
					countAnswerContainer.numLaticePoints = number_of_lattice_points;
				}
			}
			break;
		default:
			cerr << "Unknown BarvinokParameters::substitution" << endl;
			abort();
		}//end switch.

		if (read_polyhedron_data.grobner[0] == 'y')
		{

			Poly->cones = ProjectUp(Poly->cones,
					read_polyhedron_data.oldnumofvars, Poly->numOfVars,
					read_polyhedron_data.templistVec);
			Poly->numOfVars = read_polyhedron_data.oldnumofvars;

		}
		if (Print[0] == 'y')
			printListCone(Poly->cones, Poly->numOfVars);

		//   printListVector(IntegralHull(Poly->cones,  inequalities, equations, Poly->numOfVars), Poly->numOfVars);
		if (read_polyhedron_data.Memory_Save[0] == 'n')
		{

			if (read_polyhedron_data.dualApproach[0] == 'n')
			{
				cerr << "Creating generating function.\n";
				//printListVector(templistVec, oldnumofvars); cerr << ProjU << endl;
				if (read_polyhedron_data.equationsPresent[0] == 'y')
				{
					Poly->cones = ProjectUp2(Poly->cones,
							read_polyhedron_data.oldnumofvars, Poly->numOfVars,
							read_polyhedron_data.AA, read_polyhedron_data.bb);
					Poly->numOfVars = read_polyhedron_data.oldnumofvars;
				}

				createGeneratingFunctionAsMapleInput(fileName, Poly->cones,
						Poly->numOfVars);
			}
			//printListCone(cones, Poly->numOfVars);

			cerr << "Printing decomposed cones to `decomposed_cones'." << endl;
			printListConeToFile("decomposed_cones", Poly->cones,
					Poly->numOfVars);

			if (read_polyhedron_data.dualApproach[0] == 'n')
			{
				cerr << "Starting final computation.\n";
				ZZ number = Residue(Poly->cones, Poly->numOfVars);
				params->deliver_number_of_lattice_points(number);
				countAnswerContainer.numLaticePoints = number;
			}

			if (read_polyhedron_data.dualApproach[0] == 'y')
			{
				cerr << "Starting final computation.\n";
				//cerr << "output_cone: " << output_cone;
				switch (params->decomposition)
				{
				case BarvinokParameters::IrrationalPrimalDecomposition:
				case BarvinokParameters::IrrationalAllPrimalDecomposition:
				{
#ifdef HAVE_EXPERIMENTS
					ofstream out("func.rat");
					out << "HS := ";
					TrivialMonomialSubstitutionMapleOutput(out, Poly->cones, Poly->numOfVars);
					out << ";";
#else
					cerr
							<< "Trivial monomial substitution not compiled in, sorry."
							<< endl;
#endif
					break;
				}
				case BarvinokParameters::DualDecomposition:
					countAnswerContainer.seriesExpansion = ResidueFunction(Poly->cones, Poly->numOfVars, print_flag,
							read_polyhedron_data.degree, output_cone, params);
					// ResidueFunction consumes cones.
					Poly->cones = NULL;
					break;
				default:
					assert(0);
				}
				//  Else we have already computed the residue.
			}
		}
	} catch (NotIrrationalException)
	{
		cerr << "Bug: Irrationalization failed" << endl;
		THROW_LATTE(LattException::bug_Unknown);
	}; //end of try-catch.


	freeListVector(read_polyhedron_data.templistVec);
	freeListVector(read_polyhedron_data.matrix);
	delete Poly;

	if (read_polyhedron_data.rationalCone[0] == 'y')
	{
		string new_name = string(argv[argc - 1]) + ".rat";
		cerr << endl << "Rational function written to `" << new_name << "'" << endl << endl;
		rename_with_error_check("simplify.sum", new_name);
	}

	if (printfile[0] == 'y')
	{
		string new_name = string(argv[argc - 1]) + ".rat";
		cerr << endl << "Rational function written to `" << new_name << "'" << endl << endl;
		rename_with_error_check("func.rat", new_name);
	}
	if ((removeFiles[0] == 'y')
			&& (read_polyhedron_data.dualApproach[0] == 'n'))
	{

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".ext");
		// system_with_error_check(command);

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".cdd");
		// system_with_error_check(command);

		if (read_polyhedron_data.Memory_Save[0] == 'n')
		{
			// strcpy(command, "rm -f ");
			// strcat(command, fileName);
			// strcat(command, ".maple");
			// system_with_error_check(command);
		}

		// strcpy(command, "rm -f ");
		// strcat(command, fileName);
		// strcat(command, ".ead");
		// system_with_error_check(command);

	}

	//cerr << "Computation done. " << endl;

	params->total_time.stop();
	cerr << params->total_time;

	{
		// until we have a more sophisticated test script --mkoeppe
		ofstream totalTime("totalTime");
		totalTime << params->total_time.get_seconds() << " ("
				<< params->read_time.get_seconds() << "r" << ", "
				<< params->vertices_time.get_seconds() << "v" << ", "
				<< params->irrationalize_time.get_seconds() << "i" << ", "
				<< params->dualize_time.get_seconds() << "d" << ", "
				<< params->triangulate_time.get_seconds() << "t" << ", "
				<< params->decompose_time.get_seconds() << "b" << ")" << endl;
		ofstream stats("latte_stats");
		params->print_statistics(stats);
	}

	delete params;
	countAnswerContainer.checkPolynomial();
	return countAnswerContainer;
}
/* ----------------------------------------------------------------- */



