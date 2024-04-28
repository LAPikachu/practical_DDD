void lame(	char *tr,
		float *x,
		float *y,
		float *z,
		int *NMAX);

void perspective(float *x,
		float *y,
		float *z,
		float *xp,
		float *yp);

void segments(	float *x,
		float *y,
		int *XO,
		int *YO,
		int *X1,
		int *Y1);

void trace(	char *tr,
		int co,
		int *iseg,
		int *XO,
		int *YO,
		int *X1,
		int *Y1,
		int segnum,
		int *junct,
		int *nmax);

void traceps(	char *tr,
		int co,
		int *iseg,
		int *XO,
		int *YO,
		int *X1,
		int *Y1,
		int segnum,
		int *junct,
		int *nmax);

void look_seg(	int *iseg,
		int *ivec,
		int x_ptr,
		int y_ptr,
		int *XO,
		int *YO,
		int *X1,
		int *Y1,
		char *tr,
		int *junct);

void unsegout(	int *iseg,
		int *junct,
		int i);

void segout(	int *iseg,
		int *junct);

int pgdc3(	int a,
		int b,
		int c);

int mymax(	int a,
		int b);

int testevent(	float *xp,
		float *yp,
		float *zp,
		int *iseg,
		int *ivec,
		int *quit,
		int *junct,
		int *nmax,
		int *fini);

int modulo(int a, int b);

int pgdc( int aa, int bb);

void window_init(int argc, char **argv);

