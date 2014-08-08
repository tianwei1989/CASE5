void den_step(PARA_DATA *para, REAL **var);

void vel_step(PARA_DATA *para, REAL **var);

void temp_step(PARA_DATA *para, REAL **var);

void FFD_solver(PARA_DATA *para, REAL **var);

void equ_solver(PARA_DATA *para, REAL **var, int Type, REAL *x);
