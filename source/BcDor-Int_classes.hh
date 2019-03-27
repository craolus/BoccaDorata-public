//
//
//  Class to store the 2-body matrix elements of the (effective) nuclear
// interaction.
//  The matrix elements stored are in p-n formalism.
//
//  The function 'Vpp(ia,ib,ic,id,J)' returns:
//
//      < phi_a, phi_b ; J | V | phi_c, phi_d ; J >
//
//
//   ia, ib, ic and id are the sp orbits of the model space associated
//  with the object (== the interaction).
//
//

#ifndef  _INT_CLASSES_H
#define  _INT_CLASSES_H

class VppInt_t {
  public:
   ModSpace_t *MdSp;  // pointer to the associated model space
   double   *Vpp_wrk;
   int      *iabJ_list;
   int      *icdJ_list;
   int      *Jipch_str, *Jipch_end;
// int *n_Jipch; # of sorted mtx el (by the values of J, ip and ch
   int      n_me_stored, n_me_alloc;
   //
   // the following were added to have different interactions
   double **Int_list;
   int    n_ints, n_ints_alloc;
   double *StrEn;


   // functions:
          VppInt_t();
          VppInt_t(       ModSpace_t *);
          VppInt_t(char*, ModSpace_t *, int n_all_new_ints=1, int n_me_min=0);
         ~VppInt_t();
   void   set_NULL(void );
   void   init_NULL(void );
   void   Clear_all(void);
   void   Clear_all(ModSpace_t *);
   void   free_mem(void );
   void   Allocate_Vpp_table(int, int, int n_new_ints= -1);
   void   Allocate_single_Vpp(void );
   void   Deallocate_single_Vpp(int );
   int    write(char* );
   int    write_bin(char* );
   int    write_ascii(char* );
   int    read(char*, int n_all_new_ints=1, int n_me_min=0 );
   int    read_bin(char* );
   int    read_ascii(char*, int n_all_new_ints=1, int n_me_min=0 );
   int    read_Oslofmt(char*, int n_all_new_ints=1, int n_me_min=0 );

  
  // manipulation
         void scale_chrg(double , int );
         void scaleall_chrg(double , int );
  inline void scale(double C){if (NULL != Vpp_wrk) for (int i=0; i<n_me_stored; ++i) Vpp_wrk[i] *= C;}
  inline void scaleall(double C){for (int i=0; i<n_ints; ++i) for(int j=0; j<n_me_stored; ++j) Int_list[i][j] *= C;}
  inline void set_zero(void ){if (NULL != Vpp_wrk) for (int i=0; i<n_me_stored; ++i) Vpp_wrk[i] = 0.0;}
  inline void set_zeroall(void ){for (int i=0; i<n_ints; ++i) for(int j=0; j<n_me_stored; ++j) Int_list[i][j] = 0.0;}


  inline void select(int i) {if (i>=0 && i<n_ints) Vpp_wrk=Int_list[i];}
  inline int get_wrk(void ) {for(int i=0; i<n_ints; ++i) if (Vpp_wrk == Int_list[i]) return i; return -100;}
   //int    add_file(double, char* );
   int    add_file(double, char*, int n_xe=-1, int ist=0);
   void   prepare_for_reading(void );
   double get_ph(int, int, int, int, int, int*);
   double    get(int, int, int, int, int, int*);
   //int     add(int, int, int, int, int, double);
   int    add(int, int, int, int, int, double*, int n_xe=-1, int ist=0);
   int    read_line_Vpp(char* , int*, int*, int*, int*, int*, double* );
   int    check_for_double_entries(void );
   int    seek_missing_Vpp(int iadd_zeroes=0);
   int    comp_Vpp_file(char* );
   void   interpolate_en(int, int, int, int fxd_e = 0);

   int    Build_2body_me( double (*)(int, int, int, int, int) );

   int    Suppress_Tmix(void );
  };


#endif
