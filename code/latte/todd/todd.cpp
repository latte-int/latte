#include "todd/todd.h"
#include "todd/gmp_pow.h"

mpq_class todd(int dim, int k, const mpz_vector &x) {
switch (dim) {
case 1:
switch (k) {
case 0:
return -1;
case 1:
return mpq_class(1,2)*x[0];
default:
abort();
}
case 2:
switch (k) {
case 0:
return 1;
case 1:
return -mpq_class(1,2)*x[0]-mpq_class(1,2)*x[1];
case 2:
return mpq_class(1,12)*pow(x[1],2)+(mpq_class(1,4)*x[1]+mpq_class(1,12)*x[0])*x[0];
default:
abort();
}
case 3:
switch (k) {
case 0:
return -1;
case 1:
return mpq_class(1,2)*x[0]+mpq_class(1,2)*x[1]+mpq_class(1,2)*x[2];
case 2:
return -mpq_class(1,12)*pow(x[2],2)+(-mpq_class(1,4)*x[2]-mpq_class(1,12)*x[1])*x[1]+(-mpq_class(1,4)*x[2]-mpq_class(1,4)*x[1]-mpq_class(1,12)*x[0])*x[0];
case 3:
return (mpq_class(1,24)*pow(x[2],2)+mpq_class(1,24)*x[1]*x[2])*x[1]+(mpq_class(1,24)*pow(x[2],2)+(mpq_class(1,8)*x[2]+mpq_class(1,24)*x[1])*x[1]+(mpq_class(1,24)*x[1]+mpq_class(1,24)*x[2])*x[0])*x[0];
default:
abort();
}
case 4:
switch (k) {
case 0:
return 1;
case 1:
return -mpq_class(1,2)*x[0]-mpq_class(1,2)*x[1]-mpq_class(1,2)*x[2]-mpq_class(1,2)*x[3];
case 2:
return mpq_class(1,12)*pow(x[3],2)+(mpq_class(1,4)*x[3]+mpq_class(1,12)*x[2])*x[2]+(mpq_class(1,4)*x[2]+mpq_class(1,4)*x[3]+mpq_class(1,12)*x[1])*x[1]+(mpq_class(1,4)*x[1]+mpq_class(1,4)*x[2]+mpq_class(1,4)*x[3]+mpq_class(1,12)*x[0])*x[0];
case 3:
return (-mpq_class(1,24)*pow(x[3],2)-mpq_class(1,24)*x[3]*x[2])*x[2]+(-mpq_class(1,24)*pow(x[3],2)+(-mpq_class(1,8)*x[3]-mpq_class(1,24)*x[2])*x[2]+(-mpq_class(1,24)*x[2]-mpq_class(1,24)*x[3])*x[1])*x[1]+(-mpq_class(1,24)*pow(x[3],2)+(-mpq_class(1,8)*x[3]-mpq_class(1,24)*x[2])*x[2]+(-mpq_class(1,8)*x[2]-mpq_class(1,8)*x[3]-mpq_class(1,24)*x[1])*x[1]+(-mpq_class(1,24)*x[2]-mpq_class(1,24)*x[1]-mpq_class(1,24)*x[3])*x[0])*x[0];
case 4:
return -mpq_class(1,720)*pow(x[3],4)+(mpq_class(1,144)*pow(x[3],2)-mpq_class(1,720)*pow(x[2],2))*pow(x[2],2)+((mpq_class(1,48)*pow(x[3],2)+mpq_class(1,48)*x[3]*x[2])*x[2]+(mpq_class(1,144)*pow(x[3],2)+(mpq_class(1,48)*x[3]+mpq_class(1,144)*x[2])*x[2]-mpq_class(1,720)*pow(x[1],2))*x[1])*x[1]+((mpq_class(1,48)*pow(x[3],2)+mpq_class(1,48)*x[3]*x[2])*x[2]+(mpq_class(1,48)*pow(x[3],2)+(mpq_class(1,16)*x[3]+mpq_class(1,48)*x[2])*x[2]+(mpq_class(1,48)*x[2]+mpq_class(1,48)*x[3])*x[1])*x[1]+(mpq_class(1,144)*pow(x[3],2)+(mpq_class(1,48)*x[3]+mpq_class(1,144)*x[2])*x[2]+(mpq_class(1,48)*x[2]+mpq_class(1,48)*x[3]+mpq_class(1,144)*x[1])*x[1]-mpq_class(1,720)*pow(x[0],2))*x[0])*x[0];
default:
abort();
}
case 5:
switch (k) {
case 0:
return -1;
case 1:
return mpq_class(1,2)*x[0]+mpq_class(1,2)*x[1]+mpq_class(1,2)*x[2]+mpq_class(1,2)*x[3]+mpq_class(1,2)*x[4];
case 2:
return -mpq_class(1,12)*pow(x[4],2)+(-mpq_class(1,4)*x[4]-mpq_class(1,12)*x[3])*x[3]+(-mpq_class(1,4)*x[3]-mpq_class(1,4)*x[4]-mpq_class(1,12)*x[2])*x[2]+(-mpq_class(1,4)*x[3]-mpq_class(1,4)*x[2]-mpq_class(1,4)*x[4]-mpq_class(1,12)*x[1])*x[1]+(-mpq_class(1,4)*x[1]-mpq_class(1,4)*x[3]-mpq_class(1,4)*x[4]-mpq_class(1,4)*x[2]-mpq_class(1,12)*x[0])*x[0];
case 3:
return (mpq_class(1,24)*pow(x[4],2)+mpq_class(1,24)*x[4]*x[3])*x[3]+(mpq_class(1,24)*pow(x[4],2)+(mpq_class(1,8)*x[4]+mpq_class(1,24)*x[3])*x[3]+(mpq_class(1,24)*x[4]+mpq_class(1,24)*x[3])*x[2])*x[2]+(mpq_class(1,24)*pow(x[4],2)+(mpq_class(1,8)*x[4]+mpq_class(1,24)*x[3])*x[3]+(mpq_class(1,8)*x[4]+mpq_class(1,8)*x[3]+mpq_class(1,24)*x[2])*x[2]+(mpq_class(1,24)*x[2]+mpq_class(1,24)*x[4]+mpq_class(1,24)*x[3])*x[1])*x[1]+(mpq_class(1,24)*pow(x[4],2)+(mpq_class(1,8)*x[4]+mpq_class(1,24)*x[3])*x[3]+(mpq_class(1,8)*x[4]+mpq_class(1,8)*x[3]+mpq_class(1,24)*x[2])*x[2]+(mpq_class(1,8)*x[3]+mpq_class(1,8)*x[2]+mpq_class(1,8)*x[4]+mpq_class(1,24)*x[1])*x[1]+(mpq_class(1,24)*x[3]+mpq_class(1,24)*x[2]+mpq_class(1,24)*x[4]+mpq_class(1,24)*x[1])*x[0])*x[0];
case 4:
return mpq_class(1,720)*pow(x[4],4)+(-mpq_class(1,144)*pow(x[4],2)+mpq_class(1,720)*pow(x[3],2))*pow(x[3],2)+((-mpq_class(1,48)*pow(x[4],2)-mpq_class(1,48)*x[4]*x[3])*x[3]+(-mpq_class(1,144)*pow(x[4],2)+(-mpq_class(1,48)*x[4]-mpq_class(1,144)*x[3])*x[3]+mpq_class(1,720)*pow(x[2],2))*x[2])*x[2]+((-mpq_class(1,48)*pow(x[4],2)-mpq_class(1,48)*x[4]*x[3])*x[3]+(-mpq_class(1,48)*pow(x[4],2)+(-mpq_class(1,16)*x[4]-mpq_class(1,48)*x[3])*x[3]+(-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[3])*x[2])*x[2]+(-mpq_class(1,144)*pow(x[4],2)+(-mpq_class(1,48)*x[4]-mpq_class(1,144)*x[3])*x[3]+(-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[3]-mpq_class(1,144)*x[2])*x[2]+mpq_class(1,720)*pow(x[1],2))*x[1])*x[1]+((-mpq_class(1,48)*pow(x[4],2)-mpq_class(1,48)*x[4]*x[3])*x[3]+(-mpq_class(1,48)*pow(x[4],2)+(-mpq_class(1,16)*x[4]-mpq_class(1,48)*x[3])*x[3]+(-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[3])*x[2])*x[2]+(-mpq_class(1,48)*pow(x[4],2)+(-mpq_class(1,16)*x[4]-mpq_class(1,48)*x[3])*x[3]+(-mpq_class(1,16)*x[3]-mpq_class(1,16)*x[4]-mpq_class(1,48)*x[2])*x[2]+(-mpq_class(1,48)*x[3]-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[2])*x[1])*x[1]+(-mpq_class(1,144)*pow(x[4],2)+(-mpq_class(1,48)*x[4]-mpq_class(1,144)*x[3])*x[3]+(-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[3]-mpq_class(1,144)*x[2])*x[2]+(-mpq_class(1,48)*x[3]-mpq_class(1,48)*x[4]-mpq_class(1,48)*x[2]-mpq_class(1,144)*x[1])*x[1]+mpq_class(1,720)*pow(x[0],2))*x[0])*x[0];
case 5:
return (-mpq_class(1,1440)*pow(x[4],4)-mpq_class(1,1440)*x[4]*pow(x[3],3))*x[3]+(-mpq_class(1,1440)*pow(x[4],4)+(mpq_class(1,288)*pow(x[4],2)-mpq_class(1,1440)*pow(x[3],2))*pow(x[3],2)+((mpq_class(1,288)*pow(x[4],2)+mpq_class(1,288)*x[4]*x[3])*x[3]+(-mpq_class(1,1440)*x[3]-mpq_class(1,1440)*x[4])*pow(x[2],2))*x[2])*x[2]+(-mpq_class(1,1440)*pow(x[4],4)+(mpq_class(1,288)*pow(x[4],2)-mpq_class(1,1440)*pow(x[3],2))*pow(x[3],2)+((mpq_class(1,96)*pow(x[4],2)+mpq_class(1,96)*x[4]*x[3])*x[3]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]-mpq_class(1,1440)*pow(x[2],2))*x[2])*x[2]+((mpq_class(1,288)*pow(x[4],2)+mpq_class(1,288)*x[4]*x[3])*x[3]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]+(mpq_class(1,288)*x[4]+mpq_class(1,288)*x[3])*x[2])*x[2]+(-mpq_class(1,1440)*x[4]-mpq_class(1,1440)*x[2]-mpq_class(1,1440)*x[3])*pow(x[1],2))*x[1])*x[1]+(-mpq_class(1,1440)*pow(x[4],4)+(mpq_class(1,288)*pow(x[4],2)-mpq_class(1,1440)*pow(x[3],2))*pow(x[3],2)+((mpq_class(1,96)*pow(x[4],2)+mpq_class(1,96)*x[4]*x[3])*x[3]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]-mpq_class(1,1440)*pow(x[2],2))*x[2])*x[2]+((mpq_class(1,96)*pow(x[4],2)+mpq_class(1,96)*x[4]*x[3])*x[3]+(mpq_class(1,96)*pow(x[4],2)+(mpq_class(1,32)*x[4]+mpq_class(1,96)*x[3])*x[3]+(mpq_class(1,96)*x[4]+mpq_class(1,96)*x[3])*x[2])*x[2]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]+(mpq_class(1,96)*x[4]+mpq_class(1,96)*x[3]+mpq_class(1,288)*x[2])*x[2]-mpq_class(1,1440)*pow(x[1],2))*x[1])*x[1]+((mpq_class(1,288)*pow(x[4],2)+mpq_class(1,288)*x[4]*x[3])*x[3]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]+(mpq_class(1,288)*x[4]+mpq_class(1,288)*x[3])*x[2])*x[2]+(mpq_class(1,288)*pow(x[4],2)+(mpq_class(1,96)*x[4]+mpq_class(1,288)*x[3])*x[3]+(mpq_class(1,96)*x[4]+mpq_class(1,96)*x[3]+mpq_class(1,288)*x[2])*x[2]+(mpq_class(1,288)*x[4]+mpq_class(1,288)*x[2]+mpq_class(1,288)*x[3])*x[1])*x[1]+(-mpq_class(1,1440)*x[4]-mpq_class(1,1440)*x[1]-mpq_class(1,1440)*x[3]-mpq_class(1,1440)*x[2])*pow(x[0],2))*x[0])*x[0];
default:
abort();
}
default:
abort();
}
}
