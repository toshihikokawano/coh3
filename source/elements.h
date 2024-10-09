//  Z - element name conversion table

#ifndef __CSTRING_H__
#define __CSTRING_H__
#include <cstring>
#endif



const int N_ELEMENTS = 141;

static std::string element_name[N_ELEMENTS]={"",
  "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
  "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca",
  "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Uue","Ubn",
  "Ubu","Ubb","Ubt","Ubq","Ubp","Ubh","Ubs","Ubo","Ube","Utn",
  "Utu","Utb","Utt","Utq","Utp","Uth","Uts","Uto","Ute","Uqn",
};

static inline int element_getZ(char *);

/**********************************************************/
/*      Convert Element Name into Z Number                */
/**********************************************************/
static inline int element_getZ(char *d)
{
  int z = 0;
  if( isalpha(d[0]) ){
    d[0] = toupper(d[0]);
    if(strlen(d) == 2) { 
        d[1] = tolower(d[1]);
    } else if(strlen(d) == 3) { 
        d[1] = tolower(d[1]);
        d[2] = tolower(d[2]);
    }
    for(int i=1 ; i<N_ELEMENTS ; i++){
      if(element_name[i] == (std::string)d){
        z = i;
        break;
      }
    }
  }
  else z = atoi(d);
  return(z);
}
