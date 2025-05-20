#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
using namespace std;

#ifdef _USE_64
typedef int64_t INT;
#else
typedef int INT;
#endif

/**
 * Query struct.
 *
 * @param L The left index of the query
 * @param R The right index of the query (R >= L)
 * @param O The index of the minimum
 */
struct Query
{
    INT L, R, O;
};

/**
 * This is the only library function.
 *
 * @param A A rewritable array of type INT
 * @param n The array size
 * @param Q An array of type Query
 * @param q The array size
 * @return The answers are returned in Q
 */
INT rmq_offline(INT *A, INT n, Query *Q, INT q);