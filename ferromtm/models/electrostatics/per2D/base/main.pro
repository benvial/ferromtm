Include "parameters.dat";


// #############################################################################
// #############################################################################

Group {
    // Domains
    host         = Region[1000];

    If (inclusion_flag)
      incl          = Region[2000];
      Omega          = Region[{host, incl}];
    Else
      Omega          = Region[{host}];
    EndIf

    // Boundaries
  	SurfBlochLeft   = Region[101];
  	SurfBlochRight  = Region[103];
  	SurfBlochTop   = Region[102];
  	SurfBlochBot  = Region[104];

  	// Points
  	PrintPoint      = Region[10000];

}
// #############################################################################

Function{

  If (inclusion_flag)
	epsilonr[incl]            = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
		If (coupling_flag)
	 	eps_xx[host] = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ];
		eps_yy[host] = Complex[ScalarField[XYZ[], 0, 1]{2}  ,ScalarField[XYZ[], 0, 1 ]{3} ];
		eps_zz[host] = Complex[ScalarField[XYZ[], 0, 1]{4}  ,ScalarField[XYZ[], 0, 1 ]{5} ];
		epsilonr[host]    = TensorDiag[eps_xx[],eps_yy[],eps_zz[]];
		Else
			epsilonr[host]           = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
		EndIf
	Else
		epsilonr[Omega]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1];
	EndIf

    /* f[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] ; */

    /* epsilonr[host] = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
    epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
    epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
    epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];


    epsilonr[incl] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]]; */

    edir[] = Cos[theta] * Vector[1,0,0] + Sin[theta] * Vector[0,1,0];
    E_bias[] = E_static * edir[];
    /* E_bias[]  = Vector[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{6}, 0]; */
    xsi[Omega] = TensorDiag[1/CompYY[epsilonr[]], 1/CompXX[epsilonr[]],1];
}

// #############################################################################

Constraint {

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

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      /*{ NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      /*{ NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
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
              In Omega; Jacobian JVol; Integration Int_1; }
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
            { Name abs_Etot; Value { Local {[Sqrt[(CompX[E_bias[] - {d v}])^2 + (CompY[E_bias[] - {d v}])^2]] ; In Omega; Jacobian JVol; } } }

		     }
    }

}

// #############################################################################

PostOperation {

  	{ Name postop_E; NameOfPostProcessing postpro ;
  	Operation {
  		Print [ efield_tot, OnElementsOf host, Format ElementTable, File "Etot.txt"];
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
