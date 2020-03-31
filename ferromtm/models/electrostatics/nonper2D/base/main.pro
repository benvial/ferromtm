
Include "parameters.dat";

// #############################################################################
// #############################################################################

Group {
    // Domains
    /* pmlC           = Region[1]; */
  	/* pmlTB          = Region[2]; */
  	/* pmlLR          = Region[3]; */
  	host           = Region[9];
  	/* design         = Region[5]; */

    /* Omega_pml = Region[{pmlC,pmlTB,pmlLR}]; */
    design         = Region[10];
    If (inclusion_flag)
      incl          = Region[6];
      Omega_source    = Region[{incl,design}];
    Else
      Omega_source    = Region[{design}];
    EndIf
    /* Omega_nosource  = Region[{host,pmlC,pmlTB,pmlLR}]; */
    Omega           = Region[{design,incl,host}];

    Omega_gap   = Region[{design,incl}];

    Box_B = Region[100];
    Box_R = Region[200];
    Box_T = Region[300];
    Box_L = Region[400];

    // Boundaries
    SurfBlochLeft   = Region[500];
    SurfBlochRight  = Region[700];
    SurfBlochTop   = Region[600];
    SurfBlochBot  = Region[800];


    // Points
    PrintPoint	=  Region[10000];


}
// #############################################################################

Function{

    j[] = Complex[0.0, 1.0];
    // PML parameters

    /* sx[pmlC]         =   Complex[a_pml,-b_pml];
    sx[pmlTB]        =   1.;
    sx[pmlLR]        =   Complex[a_pml,-b_pml];

    sy[pmlC]  	  =   Complex[a_pml,-b_pml];
    sy[pmlTB] 	  =   Complex[a_pml,-b_pml];
    sy[pmlLR] 	  =   1.;
    sz[Omega_pml]      =   1.;

    dx = 2*h_pml  - domX_L +  domX_R;
    dy = 2*h_pml  - domY_B +  domY_T; */


    If (inclusion_flag)
      epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
      epsilonr[host]           = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
      If (coupling_flag)
        eps_xx[] = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ];
        eps_yy[] = Complex[ScalarField[XYZ[], 0, 1]{2}  ,ScalarField[XYZ[], 0, 1 ]{3} ];
        eps_zz[] = Complex[ScalarField[XYZ[], 0, 1]{4}  ,ScalarField[XYZ[], 0, 1 ]{5} ];
        epsilonr[design]    = TensorDiag[eps_xx[],eps_yy[],eps_zz[]];
      Else
        epsilonr[design]    = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];
      EndIf
    Else

    epsilonr[design]   = Complex[eps_des_re,eps_des_im] * TensorDiag[1,1,1];
    EndIf


    /* epsilonr[host]      = Complex[eps_host_re, eps_host_im] * TensorDiag[1,1,1]; */
    /* epsilonr[Omega_pml]   = eps_host_re*TensorDiag[sz[]*sy[]/sx[],sx[]*sz[]/sy[],sx[]*sy[]/sz[]]; */

    /* epsilonr[Omega_pml]      = Complex[eps_host_re, eps_host_im] * TensorDiag[1,1,1]; */


    /* epsilonr[host]      = TensorDiag[eps_xx[],eps_yy[],eps_zz[]]; */
    /* epsilonr[Omega_pml]      = TensorDiag[eps_xx[],eps_yy[],eps_zz[]]; */


    edir[] = Cos[theta] * Vector[1,0,0] + Sin[theta] * Vector[0,1,0];
    /* E_bias[Omega_source] = E_static * edir[]; */
    /* E_bias[Omega_nosource]  = 0 * edir[]; */

    lx = space2pml_L + space2pml_R + hx_des + 2*h_pml;

    E_bias[] = 0*E_static * edir[];
    /* E_bias[]  = Vector[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{6}, 0]; */
    /* xsi[Omega] = TensorDiag[1/CompYY[epsilonr[]], 1/CompXX[epsilonr[]],1]; */

}

// #############################################################################

/* Constraint {
} */

/* Constraint {

		{Name Bloch;
		    Case {
                    { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft;
                      Coefficient Complex[1.0,0.0]; Function Vector[$X-dx,$Y,$Z] ;
                    }
                    { Region SurfBlochTop; Type LinkCplx ; RegionRef SurfBlochBot;
                      Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y-dy,$Z] ;
                    }
			 }
		}
} */

Constraint {
  { Name Dirichlet; Type Assign;
    Case {
      { Region Box_R; Value 0.; }
      { Region Box_L ; Value E_static*hx_des; }
      /* { Region SurfBlochLeft; Value 1; } */
      /* { Region Box_L; Value -1.; } */
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
	  { GeoElement Point       ; NumberOfPoints  1 ; }
	  { GeoElement Line        ; NumberOfPoints  4 ; }
	  { GeoElement Triangle    ; NumberOfPoints  6 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  7 ; }
	}
      }
    }
  }
}

// #############################################################################

/* FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
   }
  }

} */

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      /* { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; } */
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      /* { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; } */
      /* { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; } */
      /* { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; } */
   }
  }

}

// #############################################################################

Formulation {
	    {Name electrostat_scalar; Type FemEquation;
    		Quantity {{ Name v; Type Local; NameOfSpace Hgrad;}}
		Equation {
          		Galerkin { [epsilonr[]*Dof{d v} , {d v}];
              In Omega; Jacobian JVol; Integration Int_1; }
              Galerkin { [ -epsilonr[] * E_bias[] , {d v} ];
              In  Omega; Jacobian JVol; Integration Int_1; }
                }
            }

        }
// #############################################################################


Resolution {
  { Name electrostat_scalar;
    System {
      { Name S; NameOfFormulation electrostat_scalar; Type ComplexValue;}
    }
    Operation {
      Generate[S] ;  Solve[S] ;SaveSolution[S] ;
    }
  }

}

// #############################################################################


PostProcessing {
    { Name postpro; NameOfFormulation electrostat_scalar; NameOfSystem S;
            Quantity {
            { Name epsilonr_xx; Value { Local { [CompXX[epsilonr[]]] ; In Omega; Jacobian JVol; } } }
            { Name epsilonr_yy; Value { Local { [CompYY[epsilonr[]]] ; In Omega; Jacobian JVol; } } }
            { Name v; Value { Local { [{v}] ; In Omega; Jacobian JVol; } } }
            { Name vx; Value { Local {[CompX[{d v}]] ; In Omega; Jacobian JVol; } } }
            { Name vy; Value { Local {[CompY[{d v}]] ; In Omega; Jacobian JVol; } } }
            { Name abs_dv; Value { Local {[Sqrt[(CompX[{d v}])^2 + (CompY[{d v}])^2]] ; In Omega; Jacobian JVol; } } }
            { Name epsilonr; Value { Local {[CompXX[epsilonr[]]] ; In Omega; Jacobian JVol; } } }
            { Name efield_tot; Value { Local { [ E_bias[] -{d v} ]  ; In Omega; Jacobian JVol; } } }
            { Name abs_Etot; Value { Local {[Sqrt[Abs[(CompX[E_bias[] - {d v}])]^2 + Abs[(CompY[E_bias[] - {d v}])]^2]] ; In Omega; Jacobian JVol; } } }
            { Name abs_Ptot; Value { Local {[Sqrt[Abs[(CompX[epsilonr[]*(E_bias[] - {d v})])]^2 + Abs[(CompY[epsilonr[]*(E_bias[] - {d v})])]^2]] ; In Omega; Jacobian JVol; } } }
            { Name polarization; Value { Local { [ epsilonr[]*(E_bias[] -{d v}) ]  ; In Omega; Jacobian JVol; } } }
            { Name int_polarization; Value { Integral { [ epsilonr[]*(E_bias[] -{d v}) ]  ; Integration Int_1; In Omega; Jacobian JVol; } } }
            { Name int_Efield; Value { Integral { [(E_bias[] -{d v}) ]  ; Integration Int_1; In Omega; Jacobian JVol; } } }

		     }
    }

}

// #############################################################################

PostOperation {

  	{ Name postop_E; NameOfPostProcessing postpro ;
  	Operation {
  		Print [ efield_tot, OnElementsOf design, Format ElementTable, File "Etot.txt"];
  		}
  	}

    { Name postop_mean; NameOfPostProcessing postpro ;
  	Operation {
      Print [ int_polarization[Omega_gap], OnElementsOf PrintPoint, Format SimpleTable, File "int_polarization.txt"];
      Print [ int_Efield[Omega_gap], OnElementsOf PrintPoint, Format SimpleTable, File "int_Efield.txt"];
      }
  	}


    { Name postop_fields_pos; NameOfPostProcessing postpro ;
        Operation {
              Print [ efield_tot , OnElementsOf Omega, File "efield_tot.pos", Name "efield_tot"];
        			Print [ v , OnElementsOf Omega, File "solution.pos", Name "v"];
              Print [ vx , OnElementsOf Omega, File "vx.pos", Name "vx"];
              Print [ vy , OnElementsOf Omega, File "vy.pos", Name "vy"];
              Print [ abs_dv , OnElementsOf Omega, File "abs_dv.pos", Name "abs_dv"];
              Print [ epsilonr , OnElementsOf Omega, File "epsilonr.pos", Name "epsilonr"];
              Print [ abs_Etot , OnElementsOf Omega, File "abs_Etot.pos", Name "abs_Etot"];
              Print [ abs_Ptot , OnElementsOf Omega, File "abs_Ptot.pos", Name "abs_Ptot"];
              Print [ polarization , OnElementsOf Omega, File "polarization.pos", Name "polarization"];

        	}
    }
    { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
            Print [ epsilonr_xx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
            { Niy-1, Nix-1} ,Format SimpleTable, File "epsilonr_xx.txt" ];
            Print [ epsilonr_yy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
            { Niy-1, Nix-1} ,Format SimpleTable, File "epsilonr_yy.txt" ];
            Print [ v , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
            { Niy-1, Nix-1} ,Format SimpleTable, File "v.txt" ];
            Print [ vx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
            { Niy-1, Nix-1} ,Format SimpleTable, File "vx.txt" ];
            Print [ vy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
            { Niy-1, Nix-1} ,Format SimpleTable, File "vy.txt" ];


        	}
    }

}


// #############################################################################
