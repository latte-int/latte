#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <rctl.h>

int main(int argc, char **argv)
{
  int seconds;
  rctlblk_t *rblk;
  if (argc != 2) {
    fprintf(stderr, "usage: set-max-cpu-time SECONDS\n");
    exit(1);
  }
  seconds = atoi(argv[1]);
  printf("%d seconds\n", seconds);
  if ((rblk = malloc(rctlblk_size())) == NULL) {
    (void) perror("rblk malloc");
    exit(1);
  }
/*   if (getrctl("task.max-cpu-time", NULL, rblk, RCTL_FIRST) == -1) { */
/*     (void) perror("getrctl"); */
/*     exit(1); */
/*   } */
  rctlblk_set_local_action(rblk, RCTL_LOCAL_SIGNAL, SIGTERM);
  //rctlblk_set_local_flags(rblk, 0);
  rctlblk_set_privilege(rblk, RCPRIV_BASIC);
  rctlblk_set_value(rblk, seconds);
  if (setrctl("process.max-cpu-time", NULL, rblk, RCTL_INSERT) != 0){
    perror("setrctl");
    exit(1);
  }
  system("env LD_LIBRARY_PATH=/localapp/imosoft/sparc-sun-solaris2.7/alpha/gmp-4.1.4-gcc33/lib:$LD_LIBRARY_PATH  ../../code/latte/count ../../EXAMPLES/magic5x5");
  return 0;
}
