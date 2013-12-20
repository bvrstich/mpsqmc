#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

ostream& operator<<(ostream& os, const Operator& theOp){

   os << "---------------------------------------------------------------------------------------" << endl;
   os << "Operator with d = " << theOp.gPhys_d() << endl;
   os << "Am I the identity? : " << theOp.AmIOpI() << endl;
   os << "Am I zero? : " << theOp.AmIOp0() << endl;
   os << "Contents : " << endl;

   for (int row=0; row<theOp.gPhys_d(); row++){
      os << "\t\t";
      for (int col=0; col<theOp.gPhys_d(); col++){
         os << theOp(row,col) << "\t";
      }
      os << endl;
   }
   os << "---------------------------------------------------------------------------------------" << endl;

   return os;
   
}

int Operator::gPhys_d() const {

   return phys_d; 

}
