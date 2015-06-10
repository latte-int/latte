Install Extra Packages
- Ipopt: https://projects.coin-or.org/Ipopt
- I think it needed
  - readline: http://cnswww.cns.cwru.edu/php/chet/readline/rltop.html
  - ncurses: http://www.gnu.org/software/ncurses/
- Or use /homes/home03/b/bedutra/installDir

My Scip+Latte files as part of my scip distro
- my public build of scip+latte 
- the core latte files are in /homes/home03/b/bedutra/scipoptsuitelatte/scip-3.1.0
- the runBuild.sh file in the above folder is my workaround for building scip+latte
- These files are copied into the latte svn in /code/latte/nonlinearOptimization/scipFiles



Making scip With IPOPT
- What worked for me: make READLINE=false IPOPT=true -j10
- I could not get it to work with READLINE=true, even after playing with LDFLAGS/CFLAGS/etc.
- Then make will error but you get a message about scip-3.1.0/lib should have some sym links for ipopt. After making the sym links, running the same make command worked for me. 





