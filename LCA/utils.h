#include "rmq-offline.h"

#define flog2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

struct Tuples
{
    INT a, p;
};

struct List{

  List * next;
  INT pos;
};

INT answer_rmqs ( INT * A, INT n, Query * Q, INT q, Query * Q_Prime,  INT * Af );
INT contract( INT * A, INT n, Query * Q, INT q, INT * AQ, INT * s, INT * Af, Query * Q_Prime );
INT recover ( INT * A, INT n, INT * AQ, INT s, INT * Af );
INT create ( INT * A, INT n, INT max, Query * Q, INT q, INT * l_0, List ** l_1, INT * AQ, INT * s, INT * Af, Query * Q_Prime);
INT marking( INT * A, INT max, Query * Q, INT q, INT * l_0, List ** l_1 );
