Include "parameters.dat";
// #############################################################################
Group {
	// SubDomains
	host         = Region[1];
	des         = Region[2];
	holes         = Region[3];
	electrode_left   = Region[4];
	electrode_right   = Region[5];
	gap   = Region[6];
	/* deswithgap   = Region[7]; */
	// Boundaries
	// Lines
	// Domains
	Omega_design   = Region[{gap}];
	/* Omega_design   = Region[{deswithgap}]; */

	Omega_gap   = Region[{holes,gap}];
	electrodes   = Region[{electrode_left,electrode_right}];
	Omega          = Region[{host, gap,des, holes, electrodes}];
	// Points
  PrintPoint	=  Region[7];

	// SurfNeumann    = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
}


// #############################################################################
Function{
	I3[]=TensorDiag[1,1,1];
	If (coupling_flag)
		eps_xx[] = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ];
		eps_yy[] = Complex[ScalarField[XYZ[], 0, 1]{2}  ,ScalarField[XYZ[], 0, 1 ]{3} ];
		eps_zz[] = Complex[ScalarField[XYZ[], 0, 1]{4}  ,ScalarField[XYZ[], 0, 1 ]{5} ];
		eps_xx_gap[] = Complex[ScalarField[XYZ[], 0, 1]{6}  ,ScalarField[XYZ[], 0, 1 ]{7} ];
		eps_yy_gap[] = Complex[ScalarField[XYZ[], 0, 1]{8}  ,ScalarField[XYZ[], 0, 1 ]{9} ];
		eps_zz_gap[] = Complex[ScalarField[XYZ[], 0, 1]{10}  ,ScalarField[XYZ[], 0, 1 ]{11} ];
		epsilon[gap]    = TensorDiag[eps_xx_gap[],eps_yy_gap[],eps_zz_gap[]];
		epsilon[des]    = TensorDiag[eps_xx[],eps_yy[],eps_zz[]];
		/* epsilon[des]    = Complex[eps_incl_re, eps_incl_im] * I3[]; */
		/* epsilon[des]    = TensorDiag[eps_xx[],eps_yy[],eps_zz[]]; */
	Else
		epsilon[gap]  = Complex[eps_des_re , eps_des_im] * I3[];
		epsilon[des]  = Complex[eps_des_re , eps_des_im] * I3[];

  EndIf
	epsilon[host]  = Complex[eps_host_re , eps_host_im] * I3[];
	epsilon[holes]  = Complex[eps_incl_re, eps_incl_im] * I3[];
	epsilon[electrodes]  = Complex[eps_electrode_re , eps_electrode_im] * I3[];
	/* epsilon[Omega_design]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1]; */


}
// #############################################################################
Constraint {

  { Name fixed_potential; Type Assign;
    Case {
      { Region electrode_left; Value Ebias*gap; }
      { Region electrode_right; Value 0; }
    }
  }
	{ Name fixed_potential2; Type Assign;
    Case {
      { Region electrode_left; Value Ebias*gap; }
      { Region electrode_right; Value 0; }
    }
  }
}
// #############################################################################
Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
  { Name JLin ;
    Case {
      { Region All ; Jacobian Lin ; }
    }
  }
}
// #############################################################################
Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints   4 ; }
          { GeoElement Line        ; NumberOfPoints  32 ; }
          { GeoElement Triangle    ; NumberOfPoints  16 ; } //1, 3, 4, 6, 7, 12, 13, 16
          { GeoElement Tetrahedron ; NumberOfPoints  29 ; }
					{ GeoElement Hexahedron  ; NumberOfPoints  34 ; }
          { GeoElement Prism       ; NumberOfPoints  51 ; }
        }
      }
    }
  }
}
// #############################################################################
FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Node;
        Support Omega; Entity NodesOf[All]; }
				If (el_order== 2)
				{ Name sn2; NameOfCoef un2; Function BF_Node_2E;
					Support Omega; Entity EdgesOf[All]; }
				EndIf
    }
		Constraint {
      { NameOfCoef un; EntityType NodesOf;   NameOfConstraint fixed_potential; }
			If (el_order== 2)
				{NameOfCoef un2; EntityType EdgesOf; NameOfConstraint fixed_potential2;}
			EndIf
    }
  }
}

// #############################################################################
Formulation {{Name estat; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hgrad;}}
		Equation { Galerkin {[epsilon[]*Dof{d u} , {d u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }

								}
            }
        }
// #############################################################################
Resolution {
  { Name estat;
    System {
      { Name S; NameOfFormulation estat; Type ComplexValue;}
    }
    Operation {
	Generate[S]; Solve[S]; SaveSolution[S];


    }
  }
}


        // Hinc[] : Complex[0,1] * 1/omega0 * 1/mu0 * Curl Einc[];
        // H_d    : Complex[0,1] * 1/omega0 * 1/mu0 * {Curl u};
// #############################################################################
PostProcessing {
    { Name postpro; NameOfFormulation estat; NameOfSystem S;
            Quantity {

							{ Name potential; Value { Local { [ {u} ]; In Omega; Jacobian JVol; } } }
							{ Name E; Value { Local { [ -{d u} ]; In Omega; Jacobian JVol; } } }
							{ Name Ex; Value { Local { [CompX[ -{d u} ]]; In Omega; Jacobian JVol; } } }
							{ Name Ey; Value { Local { [CompY[ -{d u} ]]; In Omega; Jacobian JVol; } } }
							{ Name Ez; Value { Local { [CompZ[ -{d u} ]]; In Omega; Jacobian JVol; } } }
							{ Name normE; Value { Local { [ Norm[{d u}] ]; In Omega; Jacobian JVol; } } }
							{ Name epsilon; Value { Local { [ CompXX[epsilon[]] ]; In Omega; Jacobian JVol; } } }
							{ Name int_polarization; Value { Integral { [ epsilon[]*(-{d u}) ]  ; Integration Int_1; In Omega_gap; Jacobian JVol; } } }
							{ Name int_Efield; Value { Integral { [(-{d u}) ]  ; Integration Int_1; In Omega_gap; Jacobian JVol; } } }

}}}

// #############################################################################
PostOperation {
	{ Name postop_fields_pos; NameOfPostProcessing postpro ;
			Operation {

		Print[ potential , OnElementsOf Omega, File "potential.pos"];
		Print[ E , OnElementsOf Omega, File "E.pos"];
		Print[ Ex , OnElementsOf Omega, File "Ex.pos"];
		Print[ Ey , OnElementsOf Omega, File "Ey.pos"];
		Print[ Ez , OnElementsOf Omega, File "Ez.pos"];
		Print[ normE , OnElementsOf Omega, File "normE.pos"];
		/* Print[ E0x , OnElementsOf Omega, File "E0x.pos"]; */
		/* Print[ E0y , OnElementsOf Omega, File "E0y.pos"]; */
		/* Print[ E0z , OnElementsOf Omega, File "E0z.pos"]; */
		/* Print[ E0 , OnElementsOf Omega, File "E0.pos"]; */
		/* Print[ Etotx , OnElementsOf Omega, File "Etotx.pos"]; */
		/* Print[ Etoty , OnElementsOf Omega, File "Etoty.pos"]; */
		/* Print[ Etotz , OnElementsOf Omega, File "Etotz.pos"]; */
		/* Print[ Edif , OnElementsOf Omega, File "E.pos"]; */
		/* Print[ Edifx , OnElementsOf Omega, File "Edifx.pos"]; */
		/* Print[ Edify , OnElementsOf Omega, File "Edify.pos"]; */
		/* Print[ Edifz , OnElementsOf Omega, File "Edifz.pos"]; */
	Print[ epsilon , OnElementsOf Omega, File "epsilon.pos"];
/* Print[ epsilon_annex , OnElementsOf Omega, File "epsilon_annex.pos"]; */
			}
	}

    /* { Name postopQ; NameOfPostProcessing postpro ;
        Operation {
	    Print[ normalized_losses1[Omega_source] , OnGlobal, File "Q.txt", Format Table ];
        }
    } */
		{ Name postop_epsilon; NameOfPostProcessing postpro ;
				Operation {
			Print [ epsilon , OnElementsOf des, File "epsilon.pos"];
				}
		}

		{ Name postop_mean; NameOfPostProcessing postpro ;
  	Operation {
      Print [ int_polarization[Omega_gap], OnElementsOf PrintPoint, Format SimpleTable, File "int_polarization.txt"];
      Print [ int_Efield[Omega_gap], OnElementsOf PrintPoint, Format SimpleTable, File "int_Efield.txt"];
      }
  	}

		{ Name postop_E; NameOfPostProcessing postpro ;
  	Operation {
			/* Print[E, OnElementsOf Omega_design, Depth 0, Format SimpleTable, File "E.txt" ]; */
		Print [ E, OnElementsOf des, Format ElementTable, File "Edes.txt"];
		Print [ E, OnElementsOf gap, Format ElementTable, File "Egap.txt"];
  		}
  	}

}
// #############################################################################
