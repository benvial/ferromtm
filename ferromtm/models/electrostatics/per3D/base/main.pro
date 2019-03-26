Include "parameters.dat";
Group {
  // Domains
	host         = Region[1002];
	If (inclusion_flag)
		incl         = Region[1001];
		Omega          = Region[{host,incl}];
	Else
		Omega          = Region[{host}];
	EndIf

	// Boundaries
	SurfBlochXm   = Region[505];
	SurfBlochXp  = Region[506];
	SurfBlochYm   = Region[501];
	SurfBlochYp  = Region[502];
	SurfBlochZm   = Region[504];
	SurfBlochZp  = Region[503];

	PrintPoint      = Region[10000];
}

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

	E_bias[] = E_static*Vector[1,0,0] ;
        }



Constraint {

		{Name Bloch;
		    Case {
                    { Region SurfBlochXp; Type LinkCplx ; RegionRef SurfBlochXm;
                      Coefficient Complex[1.0,0.0]; Function Vector[$X-dx,$Y,$Z] ;
                    }
                    { Region SurfBlochYp; Type LinkCplx ; RegionRef SurfBlochYm;
                      Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y-dy,$Z] ;
                    }
                    { Region SurfBlochZp; Type LinkCplx ; RegionRef SurfBlochZm;
                      Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y,$Z-dz] ;
                    }
			 }
		}
}


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

Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints   4 ; }
          { GeoElement Line        ; NumberOfPoints  32 ; }
          { GeoElement Triangle    ; NumberOfPoints  16 ; } //1, 3, 4, 6, 7, 12, 13, 16
          { GeoElement Tetrahedron ; NumberOfPoints  29 ; }
          { GeoElement Prism       ; NumberOfPoints  51 ; }
	}
      }
    }
  }
}

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
//       { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
//       { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
//       { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
   }
  }
}

Formulation {
	    {Name electrostat; Type FemEquation;
    		Quantity {{ Name v; Type Local; NameOfSpace Hgrad;}}
		Equation {
		Galerkin { [epsilonr[]*Dof{d v} , {d v}];
                 		In Omega; Jacobian JVol; Integration Int_1; }
                 Galerkin { [ -epsilonr[] * E_bias[] , {d v} ];
                  		In Omega; Jacobian JVol; Integration Int_1; }
                }
            }


        }



Resolution {
  { Name electrostat;
    System {
      { Name S; NameOfFormulation electrostat; Type ComplexValue;}
    }
    Operation {
      Generate[S] ; Solve[S] ; SaveSolution[S] ;
    }
  }

}

PostProcessing {
    { Name postpro; NameOfFormulation electrostat; NameOfSystem S;
            Quantity {
               { Name potential; Value { Local { [ {v} ]  ; In Omega; Jacobian JVol; } } }
							 { Name efield; Value { Local { [ -{d v} ]  ; In Omega; Jacobian JVol; } } }
							 { Name efield_tot; Value { Local { [ E_bias[] -{d v} ]  ; In Omega; Jacobian JVol; } } }
							 /* { Name potential_x; Value { Local { [CompX[ {v} ]]  ; In Omega; Jacobian JVol; } } } */
							 { Name norm_efield; Value { Local { [-{d v} ]  ; In Omega; Jacobian JVol; } } }
							 { Name norm_efield_tot; Value { Local { [ Norm[E_bias[]-{d v}] ]  ; In Omega; Jacobian JVol; } } }
							{ Name epsxx; Value { Local { [ CompXX[epsilonr[] ] ] ; In Omega; Jacobian JVol; } } }
							 { Name epsyy; Value { Local { [ CompYY[epsilonr[] ]]  ; In Omega; Jacobian JVol; } } }
							 { Name epszz; Value { Local { [ CompZZ[epsilonr[] ] ] ; In Omega; Jacobian JVol; } } }

		     }
    }


}

PostOperation {

	{ Name postop_E; NameOfPostProcessing postpro ;
	Operation {
		Print [ efield_tot, OnElementsOf host, Format ElementTable, File "Etot.txt"];
		}
	}
{ Name postop_fields_txt; NameOfPostProcessing postpro ;
Operation {
	Print [ potential, OnElementsOf Omega, Format TimeTable, File "potential.txt"];
	Print [ efield, OnElementsOf Omega, Format TimeTable, File "efield.txt"];
	}
}
{ Name postop_fields_pos; NameOfPostProcessing postpro ;
Operation {
	Print [ potential, OnElementsOf Omega,  File "potential.pos"];
	Print [ efield, OnElementsOf Omega,  File "efield.pos"];
	Print [ efield_tot, OnElementsOf Omega,  File "efield_tot.pos"];
	Print [ norm_efield_tot, OnElementsOf Omega,  File "norm_efield_tot.pos"];
	Print [ epszz, OnElementsOf Omega,  File "epszz.pos"];
	Print [ epsxx, OnElementsOf Omega,  File "epsxx.pos"];
	Print [ epsyy, OnElementsOf Omega,  File "epsyy.pos"];
	}
}

}
